include ../../make.inc

.SUFFIXES : .mod .o .F90

#
# Make sure the base compiler of mpif90 is the same as $(FC)
#
#mpi:    $(DQMCLIB) mpi.F90
#	mpif90  -funroll-all-loops -O3 $(PRG_FLAGS) -I$(HOME)/SRC \
#	-Bstatic -o testmpi mpi.F90 $(LIB)
all: mpi_tdmxRT_3d
mpi_tdmxRT_3d: $(DQMCLIB) dqmc_3d.o tdmxRT_mpi.o
	$(FC) $(FLAGS) -o mpi_tdmxRT_3d tdmxRT_mpi.o dqmc_3d.o $(LIB_MPI)
dqmc_3d.o: dqmc_3d.F90
	$(FC) $(FLAGS) -c dqmc_3d.F90
tdmxRT_mpi.o: tdmxRT_mpi.F90
	$(FC) $(FLAGS) -c tdmxRT_mpi.F90
clean:
	rm -f *.o;rm -f mpi_tdmx_3d;rm -f fort.*
realclean:
	rm -f *.mod;rm -f *.o;rm -f mpi_tdmx_3d;rm -f fort.*

