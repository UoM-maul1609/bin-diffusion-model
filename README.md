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
- *ncomp*: number of components, integer (currently always 2)
- *molefrac*: water molefraction, real of length kp
- *t*: temperature (K), real
- *d_self*: self diffusion coefficients (m^2s^{-1}), real length 2
- *param*: parameterisation type, integer (see below for table of inputs)
- *compound*: organic component of aerosol, integer (see below for table of inputs)

Note that not all diffusion coefficient parameterisations require all of these inputs.

param               | compound
--------------------| --------
**1** constant      | -
**2** darken        | -
**3** vignes        | -
**4** Lienhard2014  | **1** citric acid
**5** Lienhard2015  | **1** levoglucaosan / **2** levoglucosan/NH4HSO4 / **3** raffinose / **4** 3-MBTCA / **5** alpha-pinene / **6** sucrose / **7** citric acid / **8** shikimic acid
**6** Price2014     | **1** sucrose / **2** levoglucosan / **3** MgSO4 / **4** raffinose
**7** Price2015     | **1** alpha-pinene
**8** Price2016     | **1** sucrose / **2** water
**9** Zobrist2011   | **1** sucrose
**10** Shiraiwa2013 | **1** alpha-pinene


### output
- *d_coeff*: mutual diffusion coefficient (m^2s^{-1}), real of length kp

The main function prints the molefraction and corresponding diffusion coefficient to the terminal.


## diffusion_coefficients.f90 subroutines

**diffusion_coefficient**

The main subroutine that calls all others to output the mutual diffusion coefficient

`diffusion_coefficient(kp, molefrac, t, d_self, param, compound, d_coeff)`


**Lienhard2014**

'citric acid' [233-280K]

`Lienhard2014(kp, t, molefrac, compound, d_coeff)`

http://pubs.rsc.org/en/content/articlehtml/2014/cp/c4cp01939c


**Lienhard2015**

'levoglucosan' [215-253K], 'levoglucosan/NH4HSO4' [200-263K], 'raffinose' [210-300K], '3-MBTCA' [240-300K], 'alpha-pinene' [190-300K], 'sucrose' [215-310K], 'citric acid' [223-290K], 'shikimic acid' [240-300K]

`Lienhard2015(kp, t, molefrac, compound, d_coeff)`

http://www.atmos-chem-phys-discuss.net/15/24473/2015/acpd-15-24473-2015.pdf [discussion]

http://www.atmos-chem-phys.net/15/13599/2015/acp-15-13599-2015.pdf [final paper]


**Price2014**

Price 2014 only finds the diffusion coefficients for water in sucrose, levoglucosan, MgSO4 and Raffinose at a temperature of 23.5+/-0.5degreeC or 298.65+/-0.5K. 'sucrose' [23.5degreeC],'levoglucosan' [23.5degreeC],'MgSO4' [23.5degreeC],'raffinose' [23.5degreeC]

`Price2014(molefrac, compound, d_coeff)`

http://www.atmos-chem-phys.net/14/3817/2014/


**Price2015**

Price 2015 uses the parameterisation given in Lienhard 2014 to find the diffusion coefficient for alpha-pinene. Diffusion is measured in a disk. alpha-pinene [240-280K]

`Price2015(kp, t, molefrac, compound, d_coeff)`

http://pubs.rsc.org/en/content/articlehtml/2015/sc/c5sc00685f [Paper]

http://www.rsc.org/suppdata/c5/sc/c5sc00685f/c5sc00685f1.pdf [Supplementary]


**Price2016**

'sucrose' [296K], 'water' [296K]

`Price2016(molefrac, compound, d_coeff)`

http://xlink.rsc.org/?DOI=C6CP03238A


**Zobrist2011**

sucrose [160-313K]

`Zobrist2011(kp, t, molefrac, compound, d_coeff)`

http://pubs.rsc.org/en/content/articlehtml/2011/cp/c0cp01273d


**Shiraiwa2013**

Uses perculation theory to estimate the diffusion of water through 'alpha-pinene', depends upon the self-diffussion coefficients and the variables f [0.65-1.] and Z [8-16]. The diffussion coefficients are valid for 'alpha-pinene' between the water mole fractions [0.3-1.]

`Shiraiwa2013(kp, molefrac, d_self, compound, d_coeff)`

http://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp51595h



## running the code
