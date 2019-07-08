include ../../make.inc

.SUFFIXES : .mod .o .F90
all: tdm tdmx tdmxRT

tdm:    $(DQMCLIB) tdm.o
	$(FC) $(FLAGS) -o tdm tdm.o $(LIB_MPI)
tdmx:    $(DQMCLIB) tdmx.o
	$(FC) $(FLAGS) -o tdmx tdmx.o $(LIB_MPI)
tdmxRT:    $(DQMCLIB) tdmxRT.o
	$(FC) $(FLAGS) -o tdmxRT tdmxRT.o $(LIB_MPI)
tdm.o:  tdm.F90
	$(FC) $(FLAGS) -c tdm.F90
tdmx.o:  tdmx.F90
	$(FC) $(FLAGS) -c tdmx.F90
tdmxRT.o:  tdmxRT.F90
	$(FC) $(FLAGS) -c tdmxRT.F90
clean:
	rm -f tdm;rm -f tdmx;rm -f fort.*;rm -f output;rm -f tdmxRT

