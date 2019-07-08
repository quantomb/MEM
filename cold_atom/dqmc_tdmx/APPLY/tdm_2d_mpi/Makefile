include ../../make.inc

.SUFFIXES : .mod .o .F90

#
# Make sure the base compiler of mpif90 is the same as $(FC)
#
#mpi:    $(DQMCLIB) mpi.F90
#	mpif90  -funroll-all-loops -O3 $(PRG_FLAGS) -I$(HOME)/SRC \
#	-Bstatic -o testmpi mpi.F90 $(LIB)
all: mpi_tdmxRT 
mpi_tdmxRT: $(DQMCLIB) tdmxRT_mpi.o
	$(FC) $(FLAGS) -o mpi_tdmxRT tdmxRT_mpi.o $(LIB_MPI)
tdmxRT_mpi.o: tdmxRT_mpi.F90
	$(FC) $(FLAGS) -c tdmxRT_mpi.F90
clean:
	rm -f *.o;rm -f fort.*;rm -r mpi_tdmxRT
realclean:
	rm -f *.mod;rm -f *.o;rm -f mpi_tdmxRT;rm -f fort.*

