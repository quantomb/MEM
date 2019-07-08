# Makefile for dqa - library maintenance
# Author: M. Jarrell

.SUFFIXES: .a .o .for .h

#################################
#       Compiler Flags          #
#################################
# Flags for debugging
#FLAGS = -132 -CB -fpe0 -g -c
# Flags for pentium, athlon
#FLAGS = -132 -O3 -ipo -c
# Flags for recent Intel P4
FLAGS = -132 -ipo -O3 -xP -c
F90= $(INTEL_F90)  $(FLAGS) 
#################################
#       Linker Flags            #
#################################
# Flags for debugging
#LFLAGS =$(MKL_LIBS) -g
# Flags for pentium, athlon
#LFLAGS =$(MKL_LIBS) -ipo
# Flags for recent Intel P4
LFLAGS =$(MKL_LIBS) -ipo -lsvml -xP
link=$(INTEL_LINK) $(LFLAGS) 
#################################

RM = rm -f

OBJS   = readmc_linux.o 
OBJSMS = readmc_linuxms.o

all:$(OBJS) $(OBJSMS)
	$(link) $(OBJS) -o readmc_linux  
	$(link) $(OBJSMS) -o readmc_linuxms  

readmc: $(OBJS)
	$(link) -o readmc_linux *.o 
	
readmcms: $(OBJSMS)
	$(link) -o readmc_linuxms *.o
	

readmc_linux.o: readmc_linux.for
	$(F90) readmc_linux.for
	
	
readmc_linuxms.o: readmc_linuxms.for
	$(F90) readmc_linuxms.for
	
clean:
	$(RM) *.o *.il *~

realclean:
	$(RM) *.o *.il readmc_linux readmc_linuxms *~
