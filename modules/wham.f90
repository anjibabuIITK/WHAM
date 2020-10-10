PROGRAM main
USE open_mpi
USE Read_Inputs
USE wham
USE Write_Data


CALL read_input

print*, "Periodicty:",periodic   
!CALL perform_wham_2D(biased_prob)

!CALL write_fes(prob)

END PROGRAM main
