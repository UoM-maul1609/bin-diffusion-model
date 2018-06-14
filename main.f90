program main

use diffusion_coefficients

implicit none

integer :: kp, n_comp
real :: t 
real, allocatable, dimension(:) :: molefrac
real, allocatable, dimension (:) :: d_self
real, allocatable, dimension (:) :: d_coeff
character (len=20) :: param, compound


! read in namelists
call read_in_dc_namelist('namelist.in', kp, n_comp, molefrac, t, d_self, param, compound, d_coeff)


! calculate and print diffusion coefficients
call diffusion_coefficient(kp, molefrac, t, d_self, param, compound, d_coeff)


end program main
