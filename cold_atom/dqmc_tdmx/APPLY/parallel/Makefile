include ../../make.inc

.SUFFIXES : .mod .o .F90

#
# Make sure the base compiler of mpif90 is the same as $(FC)
#
#mpi:    $(DQMCLIB) mpi.F90
#	mpif90  -funroll-all-loops -O3 $(PRG_FLAGS) -I$(HOME)/SRC \
#	-Bstatic -o testmpi mpi.F90 $(LIB) 
mpiout: $(DQMCLIB) mpi.o
	$(FC) $(FLAGS) -o testmpi mpi.o $(LIB_MPI)
mpi.o: mpi.F90
	$(FC) $(FLAGS) -c mpi.F90
clean:
	rm -f testmpi;rm -f *.o;rm -f *.mod

