# Makefile for dqa - library maintenance
# Author: M. Jarrell

.SUFFIXES: .a .o .for .h

RM = rm -f
f90=ifort -xT 
flags= -132 -O3
cmflags= -132 -O3 -c
OBJS = readmc_linux.o
blaslib=-lblas
MKL=/opt/intel/Compiler/11.1/059/mkl/lib/em64t
MKLLIB=$(MKL)/libmkl_solver_lp64.a -Wl,--start-group $(MKL)/libmkl_intel_lp64.a $(MKL)/libmkl_intel_thread.a $(MKL)/libmkl_core.a $(MKL)/libmkl_blacs_lp64.a -Wl,--end-group /opt/intel/Compiler/11.1/059/lib/intel64/libiomp5.a -lpthread
all: rddcams readmc_linuxms readmc rdbin rdbin_dqmc disp sopt_dca
readmc:
	$(f90) $(flags) -o rdmc.out readmc_linux.f $(MKLLIB)
rdbin:
	$(f90) $(flags) -o rdbin.out rdbins_t.f
rdbin_dqmc: rdbins_dqmc.o
	$(f90) $(flags) -o rdbin_dqmc rdbins_dqmc.o
rdbins_dqmc.o: rdbins_dqmc.f
	$(f90) $(cmflags) -c rdbins_dqmc.f
disp:
	$(f90) $(flags) -o disp.out model_Default.f
rddcams: module_Global.o tables.o rddcams.o
	$(f90) -o rddcams module_Global.o tables.o rddcams.o $(MKLLIB)
readmc_linuxms: readmc_linuxms.o
	$(f90) -o readmc_linuxms readmc_linuxms.o $(MKLLIB)
tables.o: module_Global.f tables.f
	$(f90) $(cmflags) tables.f
module_Global.o: module_Global.f
	$(f90) $(cmflags) module_Global.f
rddcams.o: module_Global.f rddcams.f
	$(f90) $(cmflags) rddcams.f
readmc_linuxms.o: readmc_linuxms.f
	$(f90) $(cmflags) readmc_linuxms.f
sopt_dca: module_Global.o sopt_dca.o tables.o
	$(f90) $(flags) -o sopt_dca module_Global.o tables.o sopt_dca.o $(MKLLIB)
sopt_dca.o: module_Global.f sopt_dca.f
	$(f90) $(cmflags) sopt_dca.f
clean:
	$(RM) *.out; $(RM) rddcams; $(RM) readmc_linuxms;$(RM) sopt_dca
