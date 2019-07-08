############################################################################
#
#  Program:         QUEST V1.0
#
#  Module:          Makefile
#
#  Purpose:         Makefile for library
#
#  Modified:        8/20/2008
#
############################################################################

include ../make.inc

.SUFFIXES : .mod .o .F90

FLAGS = $(FC_FLAGS) 

%.o: 	%.F90
	$(FC) $(FLAGS) -c $<

SRC = 	dqmc_geom_param.F90 dqmc_latt.F90 dqmc_bonds.F90 dqmc_reclatt.F90 \
	dqmc_hamilt.F90 dqmc_symm.F90 \
	blas_mod.F90 lapack_mod.F90 dqmc_util.F90 dqmc_cfg.F90 \
	dqmc_wspace.F90 dqmc_struct.F90 dqmc_matb.F90 dqmc_ckb.F90 \
	dqmc_seqb.F90 dqmc_gfun.F90 dqmc_phy0.F90 \
	dqmc_phy2.F90 dqmc_hubbard.F90 dqmc_gtau.F90 dqmc_2dperl.F90 \
	dqmc_tdm1.F90 dqmc_tdm2.F90 dqmc_tdm.F90 dqmc_simqmc.F90 dqmc_tdmx.F90 \
	dqmc_tdmxRT.F90

OBJ = 	$(SRC:.F90=.o)

INC = 	$(SRC:.F90=.mod)

lib: 	$(OBJ)
	$(ARCH) $(ARFLAG) $(DQMCLIB) $(OBJ)
	$(RANLIB) $(DQMCLIB)

blas_mod.o: blas_mod.F90

lapack_mod.o: blas_mod.o lapack_mod.F90

dqmc_util.o: blas_mod.o lapack_mod.o dqmc_util.F90

dqmc_cfg.o: dqmc_util.o dqmc_cfg.F90

dqmc_wspace.o: lapack_mod.o dqmc_util.o	dqmc_wspace.F90

dqmc_struct.o: dqmc_util.F90

dqmc_matb.o: blas_mod.o lapack_mod.o dqmc_util.o dqmc_wspace.o dqmc_matb.F90

dqmc_ckb.o: blas_mod.o lapack_mod.o dqmc_util.o dqmc_wspace.o dqmc_ckb.F90

clean:
	rm -f *.o *.mod 
