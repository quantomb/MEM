#makefile for Linux boxes
RM = rm -f
DOBJ      = module_Global.o dsigcluster.o tables.o 
ROBJ      = module_Global.o rddca.o tables.o 
ROBJMS    = module_Global.o rddcams.o tables.o 
ROBJGW    = module_Global.o rddcagw.o tables.o
ROBJMSNB  = module_Global.o rddcams_lex.o tables.o 
SOBJ      = module_Global.o sopt_dca.o tables.o 
SPOBJ     = module_Global.o module_spectra.o tables.o spectra.o
SPSOBJ    = module_Global.o module_spectra.o tables.o spectra_star.o
EOBJ      = module_Global.o module_spectra.o tables.o energy_surface.o
NOBJ      = module_Global.o NCA2model.o
TWOOBJ    = module_Global.o module_spectra.o 2default.o tables.o
RDMCOBJ   = readmc_linuxms.o
RDMCSCOBJ = readmc_linuxms_btc.o
MFOBJ     = model_filter.o
#
#################################
#       Compiler Flags          #
#################################
INTEL_F90 = ifort
# Flags for debugging
#FLAGS = -132 -check all -c
# Flags for pentium, athlon, centrino
FLAGS = -132 -O3 -c
# Flags for recent Intel P4
#FLAGS = -132 -ipo -O3 -xP -c
F90= $(INTEL_F90)  $(FLAGS) 
#################################
#       Linker Flags            #
#################################
INTEL_LINK = ifort
# Flags for debugging
MKL= /opt/intel/Compiler/11.1/059/mkl/lib/em64t/
#MKL= /opt/intel/mkl/10.2.1.017/lib/em64t
MKL_LIBS=$(MKL)/libmkl_solver_lp64_sequential.a -Wl,--start-group $(MKL)/libmkl_intel_lp64.a $(MKL)/libmkl_sequential.a $(MKL)/libmkl_core.a -Wl,--end-group -lpthread
#LFLAGS =$(MKL_LIBS)
# Flags for pentium, athlon, centrino
#LFLAGS =-llapack -lblas
# Flags for recent Intel P4
# LFLAGS =$(MKL_LIBS)
link=$(INTEL_LINK) $(LFLAGS) 
#################################

all: $(ROBJ) $(ROBJMS) $(ROBJGW) $(ROBJMSNB) $(DOBJ) $(SOBJ) $(SPOBJ) $(SPSOBJ) $(NOBJ) $(EOBJ) $(TWOOBJ) $(RDMCOBJ) $(RDMCSCOBJ) $(MFOBJ)
	$(link) $(DOBJ) -o dsigclusteri $(MKL_LIBS);
	$(link) $(ROBJ) -o rddca $(MKL_LIBS);
	$(link) $(ROBJMS) -o rddcams $(MKL_LIBS);
	$(link) $(ROBJGW) -o rddcagw $(MKL_LIBS);
	$(link) $(ROBJMSNB) -o rddcams_lex $(MKL_LIBS);
	$(link) $(SOBJ) -o sopt_dca $(MKL_LIBS);
	$(link) $(SPOBJ) -o spectra $(MKL_LIBS);
	$(link) $(SPSOBJ) -o spectra_star $(MKL_LIBS);
	$(link) $(NOBJ) -o NCA2model $(MKL_LIBS);
	$(link) $(EOBJ) -o energy_surface $(MKL_LIBS);  
	$(link) $(TWOOBJ) -o 2default $(MKL_LIBS);
	$(link) $(RDMCOBJ) -o readmc_linuxms $(MKL_LIBS);
	$(link) $(RDMCSCOBJ) -o readmc_linuxms_btc $(MKL_LIBS);
	$(link) $(MFOBJ) -o model_filter $(MKL_LIBS)


model_filter: $(MFOBJ)
	$(link) $(MFOBJ) $(LFLAGS) -o model_filter 

readmc_linuxms: $(RDMCOBJ)
	$(link) $(RDMCOBJ) $(LFLAGS) -o readmc_linuxms 

readmc_linuxms_btc: $(RDMCSCOBJ)
	$(link) $(RDMCSCOBJ) $(LFLAGS) -o readmc_linuxms_btc 

dsigcluster: $(DOBJ)
	$(link) $(DOBJ) $(LFLAGS) -o dsigcluster  

rddca: $(ROBJ)
	$(link) $(ROBJ) -o rddca

rddcams: $(ROBJMS)
	$(link) $(ROBJMS) -o rddcams

rddcagw: $(ROBJGW)
	$(link) $(ROBJGW) -o rddcagw

rddcams_lex: $(ROBJMSNB)
	$(link) $(ROBJMSNB) -o rddcams_lex

sopt_dca: $(SOBJ)
	$(link) $(SOBJ) -o sopt_dca

spectra: $(SPOBJ)
	$(link) $(SPOBJ) -o spectra

spectra_star: $(SPSOBJ)
	$(link) $(SPSOBJ) -o spectra_star

2default: $(TWOOBJ)
	$(link) $(TWOOBJ) -o 2default

NCA2model: $(NOBJ)
	$(link) $(NOBJ) -o NCA2model

energy_surface: $(EOBJ)
	$(link) $(EOBJ) -o energy_surface

readmc_linuxms.o: readmc_linuxms.f 
	$(F90) readmc_linuxms.f

readmc_linuxms_btc.o: readmc_linuxms_btc.f 
	$(F90) readmc_linuxms_btc.f

model_filter.o: model_filter.f 
	$(F90) model_filter.f

module_Global.o: module_Global.f 
	$(F90) module_Global.f

module_spectra.o: module_spectra.f 
	$(F90) module_spectra.f

dsigcluster.o: module_Global.f dsigcluster.f
	$(F90) dsigcluster.f

tables.o: module_Global.f tables.f
	$(F90) tables.f

rddca.o: module_Global.f rddca.f
	$(F90) rddca.f

rddcams.o: module_Global.f rddcams.f
	$(F90) rddcams.f

rddcagw.o: module_Global.f rddcagw.f
	$(F90) rddcagw.f

rddcams_lex.o: module_Global.f rddcams_lex.f
	$(F90) rddcams_lex.f

sopt_dca.o: module_Global.f sopt_dca.f
	$(F90) sopt_dca.f

spectra.o: module_Global.f module_spectra.f spectra.f
	$(F90) spectra.f

spectra_star.o: module_Global.f module_spectra.f spectra_star.f
	$(F90) spectra_star.f

2default.o: module_Global.f module_spectra.f 2default.f
	$(F90) 2default.f

energy_surface.o: module_Global.f module_spectra.f energy_surface.f
	$(F90) energy_surface.f

NCA2model.o: module_Global.f NCA2model.f
	$(F90) NCA2model.f

realclean:
	$(RM) $(DOBJ) $(ROBJMS) $(ROBJGW) $(ROBJMSNB) $(ROBJ) $(SOBJ) $(SPOBJ) $(SPSOBJ) $(NOBJ) $(EOBJ) $(TWOOBJ) $(RDMCOBJ) $(RDMCSCOBJ) $(MFOBJ) global.mod module_spectra.mod dsigcluster sopt_dca rddca rddcams spectra spectra_star NCA2model energy_surface 2default readmc_linuxms model_filter rddcams_lex global.mod  module_spectra.mod *~

clean:
	$(RM) $(DOBJ) $(ROBJ) $(ROBJMS) $(ROBJGW) $(ROBJMSNB) $(SOBJ) $(SPOBJ) $(SPSOBJ) $(EOBJ) $(NOBJ) $(TWOOBJ) $(RDMCOBJ) $(RDMCSCOBJ) $(MFOBJ) global.mod  module_spectra.mod *~

