program main

use diffusion_coefficients

implicit none

integer :: kp, n_comp, i, param, compound
real :: t 
real, allocatable, dimension(:) :: molefrac
real, allocatable, dimension (:) :: d_self
real, allocatable, dimension (:) :: d_coeff


! read in namelists
call read_in_dc_namelist('namelist.in', kp, n_comp, molefrac, t, d_self, param, compound, d_coeff)


! calculate and print diffusion coefficients
call diffusion_coefficient(kp, molefrac, t, d_self, param, compound, d_coeff)

print*, 'molefrac', ':', 'd_coeff'
	do i=1,kp
		write(*,10) molefrac(i), ' : ', d_coeff(i)
	end do
10	format (f4.2,a,e9.3)


end program main
