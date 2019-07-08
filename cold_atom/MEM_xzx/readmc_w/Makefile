# Makefile for dqa - library maintenance
# Author: M. Jarrell

.SUFFIXES: .a .o .for .h

RM = rm -f
f77=ifort -xT 
flags=
blaslib=-lblas
MKL= /opt/intel/Compiler/11.1/059/mkl/lib/em64t
MKLLIB=$(MKL)/libmkl_solver_lp64.a -Wl,--start-group $(MKL)/libmkl_intel_lp64.a $(MKL)/libmkl_intel_thread.a $(MKL)/libmkl_core.a $(MKL)/libmkl_blacs_lp64.a -Wl,--end-group /opt/intel/Compiler/11.1/059/lib/intel64/libiomp5.a -lpthread
all: clean rdbins readmc readmcms default
readmc:
	$(f77) $(flags) -o runmc.out readmc_linux.f $(MKLLIB)
readmcms:
	$(f77) $(flags) -o runmcms.out readmc_linuxms.f $(MKLLIB)
rdbins:
	$(f77) $(flags) -o runbin.out rdbins.f
default:
	$(f77) $(flags) -o disp.out model_S1disper.f
clean:
	$(RM) *.out
