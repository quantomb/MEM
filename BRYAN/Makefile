# Makefile for bryan - library maintenance
# Author: M. Jarrell
# Please note that to compile this code you will likely have to change
# the mkl libary locations set by -L/opt/intel/mkl70... and to run this 
# code, the executable will have to know where the libraries are (intel 
# provides scripts such as ifortvars and mklvars to set these paths.

.SUFFIXES: .a .o .for .h
INTEL_F90=ifort 
#lapacklib=/home/zhaoxinxu/Codes/lib/lapack_dirtyrice.a
#blaslib=/home/zhaoxinxu/Codes/lib/blas_dirtyrice.a
blaslib=-lblas
lapacklib=-llapack
MKL=/opt/intel/Compiler/11.1/059/mkl/lib/em64t
MKLLIB=$(MKL)/libmkl_solver_lp64.a -Wl,--start-group $(MKL)/libmkl_intel_lp64.a $(MKL)/libmkl_intel_thread.a $(MKL)/libmkl_core.a $(MKL)/libmkl_blacs_lp64.a -Wl,--end-group -lpthread /opt/intel/Compiler/11.1/059/lib/intel64/libiomp5.a
#################################
#       Compiler Flags          #
#################################
# Flags for debugging
#FLAGS = -132 -check all
# Flags for pentium, athlon
#FLAGS = -132 -O3 -ipo
# Flags for recent Intel P4
FLAGS = -132 -ipo -xT 
# GNU Flags
#FLAGS=-ffixed-line-length-132
F90= ifort $(FLAGS)
#################################
#       Linker Flags            #
#################################
# Flags for debugging
#LFLAGS = -i-static -openmp
# Flags for pentium, athlon
#LFLAGS =$(MKL_LIBS) -ipo
# Flags for recent Intel P4
#LFLAGS =$(MKL_LIBS) -ipo -lsvml -xP
# GNU Flags
#LFLAGS =$(MKL_LIBS)
# Boudin
LFLAGS = -132 -ipo -xT
#link=$(INTEL_F90) $(LFLAGS)
link=$(INTEL_F90)
#################################
RM = rm -f

OBJS = module_global.o error.o init.o main.o\
	output.o  update.o  

bryan: $(OBJS)
#	$(link) -o bryan_linux *.o $(lapacklib) $(blaslib)
	$(link) -o bryan_linux *.o $(MKLLIB)

module_global.o: module_global.for
	$(F90) -c module_global.for

error.o: module_global.for error.for
	$(F90) -c error.for

init.o: module_global.for init.for
	$(F90) -c init.for

main.o: module_global.for main.for
	$(F90) -c main.for

output.o: module_global.for output.for
	$(F90) -c output.for

update.o: module_global.for update.for
	$(F90) -c update.for

clean:
	$(RM) *.o global.mod *.il *~

realclean:
	$(RM) *.o bryan_linux global.mod *.il *~

