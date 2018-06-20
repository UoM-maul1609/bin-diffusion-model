# diffusion-coefficients

A repository containing ideal and published diffusion coefficient parameterisations for secondary organic aerosols. Currently unlicensed, implying that copyright is the owner's. Code can be used only with permission of the owner.

## files
The repository contains the following files:
- diffusion_coefficients.f90
- main.f90
- Makefile
- namelist.in
- README.md
- gitignore

## variables

### namelist.in variables
- *kp*: number of grid points, integer
- *molefrac*: water molefraction, real of length kp
- *t*: temperature, real
- *d_self*: self diffusion coefficients, real length 2
- *param*: parameterisation type, integer (see below for table of inputs)
- *compound*: organic component of aerosol, integer (see below for table of inputs)

Note that not all diffusion coefficient parameterisations require all of these inputs.

param               | compound
--------------------| --------
**1** constant      | -
**2** darken        | -
**3** vignes        | -
**4** Lienhard2014  | 1 citric acid
**5** Lienhard2015  | 1 levoglucaosan / 2 levoglucosan/NH4HSO4 / 3 raffinose / 4 3-MBTCA / 5 alpha-pinene / 6 sucrose / 7 citric acid / 8 shikimic acid
**6** Price2014     |
**7** Price2015     |
**8** Price2016     |
**9** Zobrist2011   |
**10** Shiraiwa2013 |

### output
- *d_coeff*: mutual diffusion coefficient

The main function prints the molefraction and corresponding diffusion coefficient to the terminal.

## diffusion_coefficients.f90 subroutines


## running the code
