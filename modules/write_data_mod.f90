MODULE Write_Data
USE open_mpi
USE Read_Inputs

INTERFACE write_fes
  MODULE PROCEDURE print_pmf_2d
END SUBROUTINE write_fes

CONTAINS
  SUBROUTINE print_pmf_2d(prob)
  REAL*8, INTENT(IN) :: prob(:,:)
  LOGICAL :: parent
 integer:: i_s1, i_s2
 integer:: nbin1, nbin2, ncv
 real*8 :: grid0(3,ncv)
 real*8 :: s1, s2, kt, dum
 real*8, allocatable :: prob1(:,:)
 character (len=50):: f2

 allocate(prob1(nbin1,nbin2))
 prob1=0.0
 CALL GlobSumR(prob,prob1,nbin1*nbin2)
 if (parent)then
 f2= 'free_energy'
 open( unit =7 , file = f2, status =  'unknown' )
 do i_s1=1,nbin1 !US cv
  s1=DFLOAT(i_s1-1)*grid0(3,1)+grid0(1,1)
  do i_s2=1,nbin2 !MTD cv
     s2=DFLOAT(i_s2-1)*grid0(3,2)+grid0(1,2)
     dum= -kt*DLOG(prob1(i_s1,i_s2))
     write(7,'(4E16.8)') s1, s2, dum, prob(i_s1,i_s2)
  enddo
  write(7,*)
 enddo
 write(*,*) 'free energy written in ',f2
 close(7)
 end if
 endsubroutine print_pmf_2d



END MODULE Write_Data
