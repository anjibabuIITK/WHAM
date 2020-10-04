PROGRAM wham
USE open_mpi
USE Read_Inputs
USE wham
USE Write_Data


CALL read_input

CALL perform_wham(biased_prob)

CALL write_fes(prob)

END PROGRAM wham
