# Makefile for bryan - library maintenance
# Author: M. Jarrell
# Please note that to compile this code you will likely have to change
# the mkl libary locations set by -L/opt/intel/mkl70... and to run this 
# code, the executable will have to know where the libraries are (intel 
# provides scripts such as ifortvars and mklvars to set these paths.

.SUFFIXES: .a .o .for .h

#################################
#       Compiler Flags          #
#################################
# Flags for debugging
#FLAGS = -132 -check all
# Flags for pentium, athlon
#FLAGS = -132 -O3 -ipo
# Flags for recent Intel P4
FLAGS = -132 -ipo -O3 -xP
# GNU Flags
#FLAGS=-ffixed-line-length-132
F90= $(INTEL_F90)  $(FLAGS) 
#################################
#       Linker Flags            #
#################################
# Flags for debugging
#LFLAGS =$(MKL_LIBS)
# Flags for pentium, athlon
#LFLAGS =$(MKL_LIBS) -ipo
# Flags for recent Intel P4
LFLAGS =$(MKL_LIBS) -ipo -lsvml -xP
# GNU Flags
#LFLAGS =$(MKL_LIBS)
link=$(INTEL_LINK) $(LFLAGS) 
#################################
RM = rm -f

OBJS = module_global.o error.o init.o main.o\
	output.o  update.o  

bryan: $(OBJS)
	$(link) -o bryan_linux *.o 

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

