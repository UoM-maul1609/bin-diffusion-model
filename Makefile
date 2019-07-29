BMM_DIR = bmm
MBD_DIR = mbd
DCC_DIR = dcc
OSNF_DIR = osnf

.PHONY: bmm_code cleanall
.PHONY: mbd_code cleanall
.PHONY: dcc_code cleanall
.PHONY: osnf_code cleanall
CLEANDIRS = $(BMM_DIR) $(BMM_DIR)/osnf $(MBD_DIR) $(MBD_DIR)/osnf \
            $(DCC_DIR) $(DCC_DIR)/osnf $(OSNF_DIR) ./

DEBUG = -fbounds-check -g
MPI    =#-DMPI1
OPT    =-O3

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
#NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
#NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
NETCDF_LIB=-lnetcdff 

NETCDFLIB=-L ${NETCDF_FOR}/lib/  \
          -L ${NETCDF_C}/lib/
NETCDFMOD= ${NETCDF_FOR}/include/


FOR = gfortran -c  
FOR2 = gfortran  

AR = ar 
RANLIB = ranlib 
OBJ = o
FFLAGS = $(OPT)  $(DEBUG) -o 
FFLAGS2 =  $(DEBUG) -O3 -o 


main.exe	:  bmd_lib.a  main.$(OBJ) bmm_code mbd_code dcc_code osnf_code \
        bin_diffusion_model.$(OBJ)
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ)  \
		 bin_diffusion_model.$(OBJ) $(BMM_DIR)/bin_microphysics_module.$(OBJ) \
		 $(DCC_DIR)/diffusion_coefficients.$(OBJ) \
		 $(MBD_DIR)/diffusion.$(OBJ) $(BMM_DIR)/b_micro_lib.a $(MBD_DIR)/diff_lib.a \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG) -I${BMM_DIR} -I${DCC_DIR} \
		 -I${OSNF_DIR} 
bmd_lib.a	:   osnf_code
	cp $(OSNF_DIR)/osnf_lib.a bmd_lib.a 
bin_diffusion_model.$(OBJ)	: bin_diffusion_model.f90 bmm_code dcc_code osnf_code
	$(FOR) bin_diffusion_model.f90 -I ${NETCDFMOD}  -I${BMM_DIR} -I${MBD_DIR} -I${DCC_DIR} \
	     -I${OSNF_DIR} $(FFLAGS)bin_diffusion_model.$(OBJ)
main.$(OBJ)   : main.f90 bmm_code mbd_code osnf_code bin_diffusion_model.$(OBJ) \
             $(BMM_DIR)/bin_microphysics_module.$(OBJ)
	$(FOR)  main.f90 -I ${NETCDFMOD} -I${BMM_DIR} -I${MBD_DIR} -I${OSNF_DIR} $(FFLAGS)main.$(OBJ) 
bmm_code:
	$(MAKE) -C $(BMM_DIR)

mbd_code:
	$(MAKE) -C $(MBD_DIR)

dcc_code:
	$(MAKE) -C $(DCC_DIR)

osnf_code:
	$(MAKE) -C $(OSNF_DIR)

clean: 
	rm *.exe  *.o *.mod *~ *.a;rm -R *.dSYM

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
	
	
