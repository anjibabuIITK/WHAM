clean: rm *.o *.exe *.mod

all:
	gfortran -c open_mpi_mod.f90
	gfortran -c read_inputs_mod.f90
	gfortran -c write_data_mod.f90
	gfortran -c perform_wham_mod.f90
	gfortran -c wham.f90
	gfortran -o wham.exe  open_mpi_mod.o read_inputs_mod.o write_data_mod.o perform_wham_mod.o wham.o


