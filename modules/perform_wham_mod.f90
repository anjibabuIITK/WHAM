MODULE wham
USE open_mpi
USE Read_Inputs
USE Write_Data


IMPLICIT NONE



INTERFACE perform_wham
  MODULE PROCEDURE perform_wham_2D
END INTERFACE perform_wham


INTEGER :: iter,nbin1,nbin2
!INTEGER:: i
REAL*8 :: cnvg, dummy_1, dummy_2, dummy_3
INTEGER ::  rank, gleng1_min,gleng1_max,gleng2,ngrid
!real*8, parameter :: au_to_kcal = 627.51
!real*8, parameter :: kj_to_kcal = 0.239006





CONTAINS
  SUBROUTINE perform_wham_2D(biased_prob)
   REAL*8, INTENT(IN) :: biased_prob(:,:,:) 
real*8 :: prob(nbin1,nbin2)
!   REAL*8 :: prob(:,:) 


 if (parent)  write(*,*) 'wham begins'
 CALL DistributeGrids_2d(ncv,grid0,grid,rank,gleng1_min,gleng1_max)

 gleng2=nint((grid(2,2)-grid(1,2))/grid(3,2))+1

 write(*,*)'new_grid', gleng1_min, gleng1_max,gleng2, rank, grid(1,1)
  iter=0
  scf_loop : do
       iter=iter+1
       call wham_scf_2d(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
       if (parent) write(*,*) 'iteration #', iter, 'convergence =',cnvg
       if (mod(iter,100) .eq. 0 ) then
       call write_fes(prob)
       endif
       if((cnvg.lt.toler).or.(iter.ge.20000))then
       if (parent) write(*,*)'** convergence achieved **'
       exit scf_loop
       end if
   end do scf_loop
  
  
  END SUBROUTINE perform_wham_2D
  

subroutine wham_scf_2d(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
!performs wham scf.
implicit none
integer:: i_s1, i_s2, i_umbr, nbin1, nbin2, umbr_n,ncv
integer:: nmax(umbr_n) 
real*8 :: umbr_k(umbr_n), umbr_mean(umbr_n),dum
real*8 :: num, den, dummy_v, dummy_s1, avg, del_s1, kt, dummy, cnvg
real*8 :: prob(nbin1,nbin2),a(umbr_n)
real*8 :: grid0(3,ncv), biased_prob(nbin1,nbin2,umbr_n),grid(3,ncv)
integer:: rank, gleng1_min,gleng1_max,gleng2,ngrid
real*8,allocatable :: dummy_a1(:),dummy_a(:) 

allocate(dummy_a(umbr_n))
allocate (dummy_a1(umbr_n) )

dummy_a = 0.0d0

!calculates probability at each grid_point.
do i_s1 =gleng1_min,gleng1_max !over US cv
do i_s2 =1,gleng2 !over MTD cv
   num = 0.0d0
   den = 0.0d0
   dummy_s1 = grid0(1,1)+dfloat(i_s1-1)*grid0(3,1)

   !calculates probability.
   do i_umbr=1,umbr_n
      del_s1=dummy_s1-umbr_mean(i_umbr)
      if ( del_s1 .gt. 3.14d0 ) del_s1 = del_s1 - 6.28d0
      if ( del_s1 .lt.-3.14d0 ) del_s1 = del_s1 + 6.28d0
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      num=num+dfloat(nmax(i_umbr))*biased_prob(i_s1,i_s2,i_umbr)
      den=den+dfloat(nmax(i_umbr))*a(i_umbr)*dummy_v
   enddo

   prob(i_s1,i_s2)=num/den
   if(prob(i_s1,i_s2).ne.prob(i_s1,i_s2)) prob(i_s1,i_s2)=1.0D-16 !remove NaN
   if(prob(i_s1,i_s2)+1.eq.prob(i_s1,i_s2)) prob(i_s1,i_s2)=1.0D-16 !remove infinity

!calculate a.
      dum=grid0(3,1)*grid0(3,2)
   do i_umbr=1,umbr_n
      del_s1=dummy_s1 - umbr_mean(i_umbr)
      if ( del_s1 .ge. 3.14d0 ) del_s1 = 6.28d0 - del_s1
      if ( del_s1 .le.-3.14d0 ) del_s1 = 6.28d0 + del_s1
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      dummy_a(i_umbr)=dummy_a(i_umbr) +dum*dummy_v*prob(i_s1,i_s2)
   enddo
enddo !end of MTD cv loop
enddo !end of US cv loop

!=======================================================================================

CALL GlobSumR(dummy_a,dummy_a1,umbr_n)
do i_umbr=1,umbr_n
dummy_a(i_umbr)=dummy_a1(i_umbr)
end do
!========================================================================================

!finds convergence and update a.
 avg = 0.0d0
 cnvg = 0.0d0
 do i_umbr=1,umbr_n
 dummy_a(i_umbr) = 1.0d0/dummy_a(i_umbr)
 cnvg = cnvg + dabs(dlog(dummy_a(i_umbr))-dlog(a(i_umbr)))
 a(i_umbr) = dummy_a(i_umbr)
 enddo
 cnvg = kt*cnvg
 end subroutine wham_scf_2d
!**************************************************************************************************!

SUBROUTINE DistributeGrids_2d(ncv,grid0,grid,rank,gleng1_min,gleng1_max)
!Distribute X grid over processors by mapping grid0 to grid 
IMPLICIT NONE
INTEGER :: ncv,rank
REAL*8 :: grid0(3,ncv), grid(3,ncv)
!
INTEGER :: i,ncpu,icpu,ngrids,ngrids_m,ngrids_y,ngrids_z,ngrids_o
INTEGER :: gleng1_min, gleng1_max

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

write(*,*)'NCPUS',ncpu
DO i=1,ncv
  grid(1:3,i)=grid0(1:3,i)
END DO
rank=0

IF(ncv.EQ.1)THEN
 ngrids_y=1
 ngrids_z=1
ELSE IF(ncv.EQ.2)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=1
ELSE IF(ncv.EQ.3)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=(nint((grid0(2,3)-grid0(1,3))/grid0(3,3))+1)
END IF

if(ncpu.eq.1) then 
gleng1_min=1
gleng1_max=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
end if

!Distribute X grids
if(icpu.eq.0)WRITE(*,'(3A12,3A16)') 'CPU','CV', 'GRID SIZE', 'GRID MIN', 'GRID MAX', 'GRID BIN'
CALL Sync_procs
IF(ncpu.GT.1)THEN
  ngrids=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  write(*,*)'NEW_GRID',ngrids,icpu,ncpu
  ngrids_o=ngrids
  ngrids=ngrids/ncpu
  IF(icpu.eq.ncpu-1)THEN
    ngrids_m=ngrids+mod(ngrids_o,ncpu)
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids_m-1)*grid0(3,1)
  ELSE
    ngrids_m=ngrids
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids-1)*grid0(3,1)
  END IF
  CALL Sync_procs
  WRITE(*,'(3I12,3F16.6)') icpu, 1, ngrids_m, grid(1,1), grid(2,1), grid(3,1)
  rank=ngrids_z*ngrids_y*ngrids*icpu
  gleng1_min=ngrids*icpu+1
  gleng1_max=ngrids*(icpu+1)
  if(icpu.eq.ncpu-1) gleng1_max=ngrids*(icpu+1)+mod(ngrids_o,ncpu)
END IF 
END SUBROUTINE DistributeGrids_2d
!****************************************************************************************!


END MODULE wham
