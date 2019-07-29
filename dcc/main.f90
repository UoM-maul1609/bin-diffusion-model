	!> @mainpage
	!>@author
	!>Kathryn Fowler, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Test diffusion-coefficient calculator



	!>@author
	!>Kathryn Fowler, The University of Manchester
	!>@brief
	!>main programme reads in information, then calculates diffusion coefficients

    program main

        use diffusion_coefficients
        use numerics_type
        
        implicit none

        integer(i4b) :: kp, n_comp, i, param, compound
        real(wp) :: t 
        real(wp), allocatable, dimension(:) :: molefrac
        real(wp), allocatable, dimension (:) :: d_self
        real(wp), allocatable, dimension (:) :: d_coeff


        ! read in namelists
        call read_in_dc_namelist('namelist.in', kp, n_comp, &
            molefrac, t, d_self, param, compound, d_coeff)


        ! calculate and print diffusion coefficients
        call diffusion_coefficient(kp, 2,molefrac, t, d_self, param, compound, d_coeff)

        print*, 'molefrac', ':', 'd_coeff'
            do i=1,kp
                write(*,10) molefrac(i), ' : ', d_coeff(i)
            end do
        10	format (f4.2,a,e9.3)


    end program main
