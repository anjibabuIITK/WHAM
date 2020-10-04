MODULE Read_Inputs
USE open_mpi

  IMPLICIT NONE
  INTEGER :: ncv, umbr_n
  REAL*8, ALLOCATABLE :: grid0(:,:),a(:),umbr_mean(:), umbr_k(:),grid(:,:)
  INTEGER, ALLOCATABLE :: nmax(:)
  REAL*8 :: kt,toler, dummy
  INTEGER :: i_umbr, i_s1, i_s2
  CHARACTER(LEN=50) :: cvfile,f1
  LOGICAL :: parent, periodic=.FALSE.
  REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
  REAL*8, PARAMETER :: kj_to_kcal = 0.239006

  CONTAINS
!!
!########################################################################################################################################
SUBROUTINE read_input
CALL MPI_Start
CALL Set_Parent(parent)

! Check for the input files
if(parent) CALL Check_files()

if(parent)then
  open (unit=1, file='input',status='old')
  !read number of cv, number of umbrella,kt in energy unit
  read(1,*)ncv, umbr_n, kt
end if
  kt=kb*kt

!broadcast ncv
CALL IBcast(ncv,1)
CALL IBcast(umbr_n,1)
CALL RBcast(kt,1)

!allocate grid0
ALLOCATE(grid(3,ncv))
ALLOCATE(grid0(3,ncv))
ALLOCATE(a(umbr_n))
ALLOCATE(umbr_mean(umbr_n))
ALLOCATE(umbr_k(umbr_n))
ALLOCATE(nmax(umbr_n))

a=1.d0

if (parent) then
 !read grid_min, grid_max, grid_bin
 do i=1,ncv
    read(1,*)grid0(1:3,i)
    write(*,'(i10,3f16.6)')i,grid0(1:3,i)
 end do
end if
!! Call subroutine according to dimension
   IF(ncv==2)CALL read_whaminput_2D

END SUBROUTINE read_input
!!#########################################################################################################################################
SUBROUTINE read_whaminput_2D
  IMPLICIT NONE
  INTEGER :: nbin1, nbin2
  REAL*8, ALLOCATABLE :: biased_prob(:,:,:),prob(:,:)
  

if (parent) then
nbin1=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1 
nbin2=nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1
if (ncv.ge.3)then
STOP '3 or more CVs not implemented'
end if
end if

!broadcast grids and bin info
CALL RBcast(grid0,3*ncv)
CALL IBcast(nbin1,1)
CALL IBcast(nbin2,1)
!! Difference in allocation
ALLOCATE(biased_prob(nbin1,nbin2,umbr_n))
ALLOCATE(prob(nbin1,nbin2))

if (parent) then
   open( unit =2, file= 'whaminput', status = 'old') 
   read(2,*)toler
   do i_umbr=1, umbr_n 
      write(*,*) 'umbrella simulation #', i_umbr
      !reads probability file name
      read(2,'(a)') cvfile
      write(*,*)'CVFILE',cvfile
      !reads force constant, r0, max points in that umbrella
      read(2,*) umbr_mean(i_umbr),umbr_k(i_umbr),nmax(i_umbr)
        umbr_k(i_umbr)=umbr_k(i_umbr)*kj_to_kcal !converted in kcal
       if(umbr_mean(i_umbr).gt. 3.14d0)umbr_mean(i_umbr)=umbr_mean(i_umbr)- 6.28d0
       if(umbr_mean(i_umbr).lt.-3.14d0)umbr_mean(i_umbr)=umbr_mean(i_umbr)+ 6.28d0
 
      f1 = cvfile 
      open( unit=3, file=f1, status='old' )
      do i_s1=1,nbin1 !US
         do i_s2=1,nbin2 !MTD
         read(3,*)dummy,dummy,biased_prob(i_s1,i_s2,i_umbr)
         end do
      end do
   enddo
end if

call RBcast(toler,1)
call RBcast(umbr_mean,umbr_n)
call RBcast(umbr_k,umbr_n)
call IBcast(nmax,umbr_n)
call RBcast(biased_prob,nbin1*nbin2*umbr_n)

END SUBROUTINE read_whaminput_2D
!!##############################################################################################################################################
SUBROUTINE Check_files()
implicit none
logical :: file_exists

! Looking for input file
INQUIRE(FILE="input", EXIST=file_exists)

    if(file_exists) then
       print*," input file : Found"
    else
       print*," input file : Not Found"
       stop
   endif

! Looking for whaminput file
INQUIRE(FILE="whaminput", EXIST=file_exists)

    if(file_exists) then
       print*," whaminput file : Found"
    else
       print*," whaminput file : Not Found"
       stop
   endif



END SUBROUTINE Check_files

!!##############################################################################################################################################

SUBROUTINE ApplyPeriodicity(x)
implicit none
real*8, intent(inout)::x
REAL*8,PARAMETER :: pi=4.d0*atan(1.d0)

! Applying Periodicity

   if (x .gt. pi ) x = x - 2.d0*pi
   if (x .lt.-pi ) x = x + 2.d0*pi

RETURN
END SUBROUTINE


!!##############################################################################################################################################





END MODULE Read_Inputs
