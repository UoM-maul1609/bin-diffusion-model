BMM_DIR = bmm
MBD_DIR = mbd
DCC_DIR = dcc

BMM_SCE_DIR = $(BMM_DIR)/sce
BMM_OPT_DIR = $(BMM_DIR)/opt

.PHONY: bmm_code mbd_code dcc_code clean cleanall
CLEANDIRS = $(BMM_DIR) $(BMM_DIR)/osnf $(BMM_DIR)/sce $(BMM_DIR)/sce/osnf \
            $(MBD_DIR) $(MBD_DIR)/osnf $(DCC_DIR) $(DCC_DIR)/osnf ./

DEBUG = -g
MPI    =#-DMPI1
OPT    =-O3

# These three lines should be edited for your system. On systems
# that do not have separate Fortran and C libraries, set NETCDF_FOR and NETCDF_C
# to the same prefix and NETCDF_LIB to -lnetcdf as appropriate.
NETCDF_LIB=-lnetcdff

NETCDFLIB=-L ${NETCDF_FOR}/lib/ \
          -L ${NETCDF_C}/lib/
NETCDFMOD= ${NETCDF_FOR}/include/

FOR = gfortran -c
FOR2 = gfortran

AR = ar
RANLIB = ranlib
OBJ = o
FFLAGS = $(OPT) $(DEBUG) -o
FFLAGS2 = $(DEBUG) -O3 -o

BMM_INCLUDES = -I$(BMM_DIR) -I$(BMM_DIR)/osnf -I$(BMM_SCE_DIR) -I$(BMM_OPT_DIR)
BMM_LINK = $(BMM_DIR)/bin_microphysics_module.$(OBJ) \
           $(BMM_DIR)/b_micro_lib.a \
           $(BMM_DIR)/osnf/osnf_lib.a \
           $(BMM_SCE_DIR)/sce_micro_lib.a \
           $(BMM_SCE_DIR)/sce_module.$(OBJ) \
           $(BMM_OPT_DIR)/optics.a

main.exe: bmd_lib.a main.$(OBJ) bmm_code mbd_code dcc_code \
          bin_diffusion_model.$(OBJ)
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) \
		 bin_diffusion_model.$(OBJ) $(BMM_LINK) \
		 $(DCC_DIR)/diffusion_coefficients.$(OBJ) \
		 $(MBD_DIR)/diffusion.$(OBJ) $(MBD_DIR)/diff_lib.a \
		 -lm ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG) \
		 $(BMM_INCLUDES) -I$(DCC_DIR)

bmd_lib.a: bmm_code
	cp $(BMM_DIR)/osnf/osnf_lib.a bmd_lib.a

bin_diffusion_model.$(OBJ): bin_diffusion_model.f90 bmm_code dcc_code mbd_code
	$(FOR) bin_diffusion_model.f90 -I ${NETCDFMOD} $(BMM_INCLUDES) \
	     -I$(MBD_DIR) -I$(DCC_DIR) $(FFLAGS)bin_diffusion_model.$(OBJ)

main.$(OBJ): main.f90 bmm_code mbd_code bin_diffusion_model.$(OBJ)
	$(FOR) main.f90 -I ${NETCDFMOD} $(BMM_INCLUDES) -I$(MBD_DIR) \
	     $(FFLAGS)main.$(OBJ)

bmm_code:
	$(MAKE) -C $(BMM_DIR)

mbd_code:
	$(MAKE) -C $(MBD_DIR)

dcc_code:
	$(MAKE) -C $(DCC_DIR)

clean:
	rm -f *.exe *.o *.mod *~ *.a
	rm -rf *.dSYM

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean || true; \
	done
