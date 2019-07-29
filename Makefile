OSNF_DIR = osnf

.PHONY: osnf_code clearnall
CLEANDIRS = $(OSNF_DIR) ./


DEBUG = -fbounds-check -g 
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
FFLAGS = $(OPT)  $(DEBUG)  -o 
FFLAGSOMP = -fopenmp-simd $(FFLAGS)
FFLAGS2 =  $(DEBUG) -O3 -o 


main.exe	:  main.$(OBJ) diffusion.$(OBJ)  \
			 diff_lib.a osnf_code 
	$(FOR2) $(FFLAGS)main.exe main.$(OBJ) diffusion.$(OBJ)  \
			 -lm diff_lib.a \
		 ${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
diff_lib.a   :  osnf_code 
	cp $(OSNF_DIR)/osnf_lib.a diff_lib.a 
diffusion.$(OBJ) : diffusion.f90 osnf_code diff_lib.a
	$(FOR) diffusion.f90 -I ${NETCDFMOD} -I${OSNF_DIR} $(FFLAGSOMP)diffusion.$(OBJ)
main.$(OBJ)   : main.f90 diffusion.$(OBJ)  
	$(FOR)  main.f90 -I ${NETCDFMOD} -I${OSNF_DIR} $(FFLAGS)main.$(OBJ) 

osnf_code:
	$(MAKE) -C $(OSNF_DIR)
clean :
	rm *.exe *.o *.mod *~ \
	*.a;rm -R *.dSYM
cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
