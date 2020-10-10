
all:
	mpif90 -c open_mpi_mod.f90
	mpif90 -c read_inputs_mod.f90
	mpif90 -c write_data_mod.f90
	mpif90 -c perform_wham_mod.f90
	mpif90 -c wham.f90
	mpif90 -o wham.exe  open_mpi_mod.o read_inputs_mod.o write_data_mod.o perform_wham_mod.o wham.o


