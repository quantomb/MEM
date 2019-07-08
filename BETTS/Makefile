#makefile for Linux boxes
RM = rm -f
DOBJ      = module_Global.o dsigcluster.o tables.o 
D2OBJ     = module_Global.o dsigcluster_sig.o tables.o
DMOBJ     = module_Global.o dsigcluster_M.o tables.o
ROBJ      = module_Global.o rddca.o tables.o 
ROBJSIG   = module_Global.o rddca_sig.o tables.o
ROBJMS    = module_Global.o rddcams.o tables.o
ROBJMSNB  = module_Global.o rddcams_lex.o tables.o 
SOBJ      = module_Global.o sopt_dca.o tables.o 
SPOBJ     = module_Global.o module_spectra.o tables.o spectra.o
PRINTSIG  = module_Global.o tables.o rddca_print_sig.o
INTER_SIGKWN = module_Global.o module_spectra.o tables.o interpolate_sigkwn.o
ETMJOBJ   = module_Global.o module_spectra.o tables.o eta_MJ.o
SPSOBJ    = module_Global.o module_spectra.o tables.o spectra_star.o
EOBJ      = module_Global.o module_spectra.o tables.o energy_surface.o
NOBJ      = module_Global.o NCA2model.o
TWOOBJ    = module_Global.o module_spectra.o 2default.o tables.o
TWOCHIPD    = module_Global.o module_spectra.o 2chipp_d.o tables.o
RDMCOBJ   = readmc_linuxms.o
RDMCSCOBJ = readmc_linuxms_btc.o
MFOBJ     = model_filter.o
#
#################################
#       Compiler Flags          #
#################################
INTEL_F90 = ifort
# Flags for debugging
# FLAGS = -132 -check all -traceback -c
# Flags for pentium, athlon, centrino
# FLAGS = -132 -O3 -c
FLAGS = -132 -fast -c
# Flags for recent Intel P4
# FLAGS = -132 -ipo -O3 -xP -c
#
#
F90= $(INTEL_F90)  $(FLAGS) 
F90T= $(INTEL_F90) -132  -c
#################################
#       Linker Flags            #
#################################
INTEL_LINK = ifort
# Flags for debugging
#LFLAGS =$(MKL_LIBS)
# Flags for pentium, athlon, centrino
#LFLAGS =$(MKL_LIBS)
LFLAGS=-llapack -lblas 
# Flags for recent Intel P4
#LFLAGS =$(MKL_LIBS) -ipo -lsvml -xP
link=$(INTEL_LINK) $(LFLAGS) 
#################################

all: $(ROBJ) $(ROBJMS) $(ROBJSIG) $(ROBJMSNB) $(DOBJ) $(D2OBJ) $(SOBJ) $(SPOBJ) $(ETOBJ) $(ETMJOBJ) $(SPSOBJ) $(NOBJ) $(EOBJ) $(TWOOBJ) $(RDMCOBJ) $(RDMCSCOBJ) $(MFOBJ) $(DMOBJ)
	$(link) $(DOBJ) -o dsigcluster;
	$(link) $(D2OBJ) -o dsigcluster_sig;
	$(link) $(DMOBJ) -o dsigcluster_M;
	$(link) $(ROBJ) -o rddca;
	$(link) $(ROBJSIG) -o rddca_sig;
	$(link) $(ROBJMS) -o rddcams;
	$(link) $(ROBJMSNB) -o rddcams_lex;
	$(link) $(SOBJ) -o sopt_dca;
	$(link) $(SPOBJ) -o spectra;
	$(link) $(ETMJOBJ) -o eta_MJ;
	$(link) $(SPSOBJ) -o spectra_star;
	$(link) $(NOBJ) -o NCA2model;
	$(link) $(EOBJ) -o energy_surface; 
	$(link) $(TWOOBJ) -o 2default;
	$(link) $(TWOCHIPD) -o 2chipp_d;
	$(link) $(RDMCOBJ) -o readmc_linuxms;
	$(link) $(RDMCSCOBJ) -o readmc_linuxms_btc;
	$(link) $(MFOBJ) -o model_filter


model_filter: $(MFOBJ)
	$(link) $(MFOBJ) $(LFLAGS) -o model_filter 

readmc_linuxms: $(RDMCOBJ)
	$(link) $(RDMCOBJ) $(LFLAGS) -o readmc_linuxms 

readmc_linuxms_btc: $(RDMCSCOBJ)
	$(link) $(RDMCSCOBJ) $(LFLAGS) -o readmc_linuxms_btc 

dsigcluster: $(DOBJ)
	$(link) $(DOBJ) $(LFLAGS) -o dsigcluster  

dsigcluster_sig: $(D2OBJ)
	$(link) $(D2OBJ) $(LFLAGS) -o dsigcluster_sig  

dsigcluster_M: $(DMOBJ)
	$(link) $(DMOBJ) $(LFLAGS) -o dsigcluster_M  

rddca: $(ROBJ)
	$(link) $(ROBJ) -o rddca

rddca_sig: $(ROBJSIG)
	$(link) $(ROBJSIG) -o rddca_sig

rddcams: $(ROBJMS)
	$(link) $(ROBJMS) -o rddcams

rddcams_lex: $(ROBJMSNB)
	$(link) $(ROBJMSNB) -o rddcams_lex

sopt_dca: $(SOBJ)
	$(link) $(SOBJ) -o sopt_dca

spectra: $(SPOBJ)
	$(link) $(SPOBJ) -o spectra

eta_MJ: $(ETMJOBJ)
	$(link) $(ETMJOBJ) -o eta_MJ

spectra_star: $(SPSOBJ)
	$(link) $(SPSOBJ) -o spectra_star

2default: $(TWOOBJ)
	$(link) $(TWOOBJ) -o 2default

2chipp_d: $(TWOCHIPD)
	$(link) $(TWOCHIPD) -o 2chipp_d

NCA2model: $(NOBJ)
	$(link) $(NOBJ) -o NCA2model

energy_surface: $(EOBJ)
	$(link) $(EOBJ) -o energy_surface

rddca_print_sig: $(PRINTSIG)
	$(link) $(PRINTSIG) -o rddca_print_sig

interpolate_sigkwn: $(INTER_SIGKWN)
	$(link) $(INTER_SIGKWN) -o interpolate_sigkwn

rddca_print_sig.o: rddca_print_sig.f
	$(F90) rddca_print_sig.f

interpolate_sigkwn.o: module_Global.f module_spectra.f interpolate_sigkwn.f
	$(F90) interpolate_sigkwn.f

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

dsigcluster_sig.o: module_Global.f dsigcluster_sig.f
	$(F90) dsigcluster_sig.f

dsigcluster_M.o: module_Global.f dsigcluster_M.f
	$(F90) dsigcluster_M.f

tables.o: module_Global.f tables.f
	$(F90T) tables.f

rddca.o: module_Global.f rddca.f
	$(F90) rddca.f

rddca_sig.o: module_Global.f rddca_sig.f
	$(F90) rddca_sig.f

rddcams.o: module_Global.f rddcams.f
	$(F90) rddcams.f

rddcams_lex.o: module_Global.f rddcams_lex.f
	$(F90) rddcams_lex.f

sopt_dca.o: module_Global.f sopt_dca.f
	$(F90) sopt_dca.f

spectra.o: module_Global.f module_spectra.f spectra.f
	$(F90) spectra.f

eta_MJ.o: module_Global.f module_spectra.f eta_MJ.f
	$(F90) eta_MJ.f

spectra_star.o: module_Global.f module_spectra.f spectra_star.f
	$(F90) spectra_star.f

2default.o: module_Global.f module_spectra.f 2default.f
	$(F90) 2default.f

2chipp_d.o: module_Global.f module_spectra.f 2chipp_d.f
	$(F90) 2chipp_d.f

energy_surface.o: module_Global.f module_spectra.f energy_surface.f
	$(F90) energy_surface.f

NCA2model.o: module_Global.f NCA2model.f
	$(F90) NCA2model.f

realclean:
	$(RM) $(DOBJ) $(D2OBJ) $(DMOBJ) $(ROBJMS) $(ROBJMSNB) $(ROBJ) $(SOBJ) $(SPOBJ) $(ETOBJ) $(ETMJOBJ) $(SPSOBJ) $(NOBJ) $(EOBJ) $(TWOOBJ) $(TWOCHIPD) $(RDMCOBJ) $(RDMCSCOBJ) $(MFOBJ) global.mod module_spectra.mod dsigcluster sopt_dca rddca rddcams spectra eta_MJ spectra_star NCA2model energy_surface 2default 2chipp_d readmc_linuxms model_filter rddcams_lex global.mod  module_spectra.mod *~

clean:
	$(RM) $(DOBJ) $(D2OBJ) $(DMOBJ) $(ROBJ) $(ROBJMS) $(ROBJMSNB) $(SOBJ) $(SPOBJ) $(ETOBJ) $(ETMJOBJ) $(SPSOBJ) $(EOBJ) $(NOBJ) $(TWOOBJ) $(TWOCHIPD) $(RDMCOBJ) $(RDMCSCOBJ) $(MFOBJ) global.mod  module_spectra.mod *~

