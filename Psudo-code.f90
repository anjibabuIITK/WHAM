! Module based Fortran program to perform wham for multiple dimensions
! 
!
PROGRAM wham
use openmpi
use ReadInputs
use wham
use WriteData

! Check input files exist or not
call CheckFiles()

! Read the Input files
Call ReadInput()

! Perform the WHAM
call perform_wham()

! write the final data to file
call Write_fes() 


END PROGRAM wham
