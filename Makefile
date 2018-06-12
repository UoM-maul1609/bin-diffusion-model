# this is an initial attempt at a makefile, needs to be cleaned up to follow the structure of other makefiles


# variables

DEBUG = -fbounds-check -g
OPT = -O3

FOR = gfortran -c
FOR2 = gfortran

OBJ = o
FFLAGS = $(OPT)  $(DEBUG) -o 
FFLAGS2 = $(DEBUG) -O3 -o


# dependancies

main.exe	: main.$(OBJ) diffusion_coefficients.$(OBJ)
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) diffusion_coefficients.$(OBJ)

diffusion_coefficients.$(OBJ)	: diffusion_coefficients.f90
	$(FOR) diffusion_coefficients.f90 \
		$(FFLAGS)diffusion_coefficients.$(OBJ)

main.$(OBJ)	 : main.f90 diffusion_coefficients.$(OBJ)
	$(FOR) main.f90 $(FFLAGS)main.$(OBJ)


clean: 
	rm *.exe *.o *.mod *~ *.a
