include ../../make.inc

.SUFFIXES : .mod .o .F90

#
# Make sure the base compiler of mpif90 is the same as $(FC)
#
#mpi:    $(DQMCLIB) mpi.F90
#	mpif90  -funroll-all-loops -O3 $(PRG_FLAGS) -I$(HOME)/SRC \
#	-Bstatic -o testmpi mpi.F90 $(LIB) 
mpi_3d:    $(DQMCLIB) dqmc_3d.o mpi_3d.o
	$(FC) $(FLAGS) -o 3dmpi mpi_3d.o dqmc_3d.o $(LIB_MPI)
mpi_3d.o: dqmc_3d.F90
	$(FC) $(FLAGS) -c mpi_3d.F90
dqmc_3d.o: dqmc_3d.F90
	$(FC) $(FLAGS) -c dqmc_3d.F90
clean:
	rm -f 3dmpi;rm -f *.o;rm -f *.mod

