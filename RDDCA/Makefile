#makefile for Linux boxes
RM = rm -f
DOBJ   = module_Global.o dsigcluster.o tables.o 
ROBJ   = module_Global.o rddca.o tables.o 
RSCOBJ = module_Global.o rddca_sc.o tables.o 
SOBJ   = module_Global.o sopt_dca.o tables.o 
SPOBJ  = module_Global.o module_spectra.o tables.o spectra.o
EOBJ   = module_Global.o module_spectra.o tables.o energy_surface.o
NOBJ   = module_Global.o NCA2model.o
TWOOBJ = module_Global.o module_spectra.o 2default.o tables.o
# Set flags
  FLAGS =-c -132
  LFLAGS = 
  F77=$(INTEL_F90) $(FLAGS) 
  link= $(INTEL_F90)

all: $(ROBJ) $(DOBJ) $(SOBJ) $(SPOBJ) $(NOBJ) $(EOBJ) $(TWOOBJ) $(RSCOBJ)
	$(link) $(DOBJ) -o dsigcluster;
	$(link) $(ROBJ) -o rddca;
	$(link) $(RSCOBJ) -o rddca_sc;
	$(link) $(SOBJ) -o sopt_dca;
	$(link) $(SPOBJ) -o spectra;
	$(link) $(NOBJ) -o NCA2model;
	$(link) $(EOBJ) -o energy_surface; 
	$(link) $(TWOOBJ) -o 2default


dsigcluster: $(DOBJ)
	$(link) $(DOBJ) $(LFLAGS) -o  dsigcluster  

rddca: $(ROBJ)
	$(link) $(ROBJ) -o rddca

rddca_sc: $(RSCOBJ)
	$(link) $(RSCOBJ) -o rddca_sc

sopt_dca: $(SOBJ)
	$(link) $(SOBJ) -o sopt_dca

spectra: $(SPOBJ)
	$(link) $(SPOBJ) -o spectra

2default: $(TWOOBJ)
	$(link) $(TWOOBJ) -o 2default

NCA2model: $(NOBJ)
	$(link) $(NOBJ) -o NCA2model

energy_surface: $(EOBJ)
	$(link) $(EOBJ) -o energy_surface

module_Global.o: module_Global.f 
	$(F77) module_Global.f

module_spectra.o: module_spectra.f 
	$(F77) module_spectra.f

dsigcluster.o: module_Global.f dsigcluster.f
	$(F77) dsigcluster.f

tables.o: module_Global.f tables.f
	$(F77) tables.f

rddca.o: module_Global.f rddca.f
	$(F77) rddca.f

rddca_sc.o: module_Global.f rddca_sc.f
	$(F77) rddca_sc.f

sopt_dca.o: module_Global.f sopt_dca.f
	$(F77) sopt_dca.f

spectra.o: module_Global.f module_spectra.f spectra.f
	$(F77) spectra.f

2default.o: module_Global.f module_spectra.f 2default.f
	$(F77) 2default.f

energy_surface.o: module_Global.f module_spectra.f energy_surface.f
	$(F77) energy_surface.f

NCA2model.o: module_Global.f NCA2model.f
	$(F77) NCA2model.f

realclean:
	$(RM) $(DOBJ) $(ROBJ) $(SOBJ) $(RSCOBJ) $(SPOBJ) $(NOBJ) $(EOBJ) $(TWOOBJ) global.mod module_spectra.mod dsigcluster sopt_dca rddca rddca_sc spectra NCA2model energy_surface 2default global.mod  module_spectra.mod *~

clean:
	$(RM) $(DOBJ) $(ROBJ) $(SOBJ) $(RSCOBJ) $(SPOBJ) $(EOBJ) $(NOBJ) $(TWOOBJ) global.mod  module_spectra.mod *~

