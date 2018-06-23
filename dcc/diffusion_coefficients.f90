!>@author
!>Kathryn Fowler, The University of Manchester
!>@brief
!>module for SOA diffusion coefficient parameterisations


    module diffusion_coefficients
    use nrtype

    private
    public :: read_in_dc_namelist, diffusion_coefficient

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! read_in_namelist							       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] nmlfile: namelist file
    !>@param[out] kp: number of grid points / shells
    !>@param[out] n_comp: number of components
    !>@param[out] molefrac: water molefraction
    !>@param[out] t: temperature
    !>@param[out] d_self: self diffusion coefficients
    !>@param[out] param: parameterisation type
    !>@param[out] compound: organic component of aerosol
    !>@param[out] d_coeff: mutual diffusion coefficient
    subroutine read_in_dc_namelist(nmlfile, kp, n_comp, &
                                molefrac, t, d_self, param, compound, d_coeff)
        implicit none
        integer(i4b), intent(out) :: kp, n_comp, param, compound
        real(sp), intent(out) :: t
        real(sp), allocatable, dimension(:), intent(out) :: molefrac
        real(sp), allocatable, dimension (:), intent(out) :: d_coeff
        real(sp), allocatable, dimension (:), intent(out) :: d_self
        character (len=*), intent(in) :: nmlfile

        ! define namelists
        namelist /dc_setup/ kp, n_comp
        namelist /dc_vars/ molefrac, t, d_self, param, compound

        ! read in namelist
        open(10, file=nmlfile, status='old', recl=80, delim='apostrophe')
        read(10, nml=dc_setup)			
            allocate (molefrac(1:kp))
            allocate (d_self(1:n_comp))
            allocate (d_coeff(1:kp))
        read(10, nml=dc_vars)	

        close(10)

    end subroutine read_in_dc_namelist


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! diffusion_coefficient							       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] t: temperature
    !>@param[in] d_self: self diffusion coefficients
    !>@param[in] param: parameterisation type
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine diffusion_coefficient(kp, molefrac, t, &
                                d_self, param, compound, d_coeff)
        implicit none
        integer(i4b), intent(in) :: kp, param, compound
        real(sp), intent(in) :: t
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp), dimension (kp), intent(in) :: d_self


        select case (param)

            case (1) !constant
            ! constant diffusion coefficient
                d_coeff(:)=d_self(1)

            case (2) !darken
            ! linear relation between water molefraction and diffusion coefficient
                d_coeff=d_self(1)*(1-molefrac)+d_self(2)*(molefrac)

            case (3) !vignes
            ! logarithmic relation between water mole fraction and diffusion coefficient
                d_coeff=d_self(1)**(1-molefrac)*d_self(2)**(molefrac)

            case(4) !Lienhard2014
            ! http://pubs.rsc.org/en/content/articlehtml/2014/cp/c4cp01939c
                call Lienhard2014(kp, t, molefrac, compound, d_coeff)

            case(5) !Lienhard2015
            ! http://www.atmos-chem-phys-discuss.net/15/24473/2015/acpd-15-24473-2015.pdf [discussion]
            ! http://www.atmos-chem-phys.net/15/13599/2015/acp-15-13599-2015.pdf [final paper]
                call Lienhard2015(kp, t, molefrac, compound, d_coeff)

            case(6) !Price2014
            ! http://www.atmos-chem-phys.net/14/3817/2014/
                call Price2014(kp,molefrac, compound, d_coeff)

            case(7) !Price2015
            ! http://pubs.rsc.org/en/content/articlehtml/2015/sc/c5sc00685f [Paper]
            ! http://www.rsc.org/suppdata/c5/sc/c5sc00685f/c5sc00685f1.pdf [Supplementary]
                call Price2015(kp, t, molefrac, compound, d_coeff)

            case(8) !Price2016
            ! http://xlink.rsc.org/?DOI=C6CP03238A	
                call Price2016(kp,molefrac, compound, d_coeff)

            case(9) !Zobrist2011
            ! http://pubs.rsc.org/en/content/articlehtml/2011/cp/c0cp01273d
                call Zobrist2011(kp, t, molefrac, compound, d_coeff)

            case(10) !Shiraiwa2013
            ! http://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp51595h 
                call Shiraiwa2013(kp, molefrac, d_self, compound, d_coeff)

            case default
                print*, "selected param not found"

        end select


    end subroutine diffusion_coefficient


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lienhard2014								       !
    ! http://pubs.rsc.org/en/content/articlehtml/2014/cp/c4cp01939c 	       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] t: temperature
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine Lienhard2014(kp, t, molefrac, compound, d_coeff)
        implicit none
        integer(i4b), intent(in) :: kp, compound
        real(sp), intent(in) :: t
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: d_coeffit, d_molefrac, c, d
        real(sp), dimension(kp) :: alpha


        ! citric acid
        if (compound == 1) then
            d_coeffit = 10.**(-15.-(175./(t-208.)))
            d_molefrac = 10.**(-6.514-(387.4/(t-118.)))
            if (t>265) then
                c = -41.+0.143*265.
            else
                c = -41.+0.143*t
            end if
            if (t>255) then
                d = -69.+0.28*255.
                else
                    d = -69.+0.28*t
            end if
            alpha = exp((1.-molefrac)**2.*(c+3.*d-4.*d*(1.-molefrac)))
                d_coeff = d_molefrac**(molefrac*alpha)*d_coeffit**(1.-molefrac*alpha)
        else
            print*, 'selected compound not found'
        end if

    end subroutine Lienhard2014


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lienhard2015								       !	
    ! http://www.atmos-chem-phys-discuss.net/15/24473/2015/acpd-15-24473-2015.pdf  !
    ! [discussion]								       !
    ! http://www.atmos-chem-phys.net/15/13599/2015/acp-15-13599-2015.pdf 	       !
    ! [final paper]								       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] t: temperature
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine Lienhard2015(kp, t, molefrac, compound, d_coeff)
        implicit none
        integer(i4b), intent(in) :: kp, compound
        real(sp), intent(in) :: t
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: logdwtg0, eact, tg, a, a1, a2, b, b1, b2, &
            t0, t1, t2, ta, tb, s, zeta_a_0_aw0, zeta_a_aw0, &
            zeta_v_0, zeta_v, dwt0
        real(sp), dimension(kp) :: alpha, dwt1

    

        select case (compound)
    
        case (1) !levoglucosan
                logdwtg0=-27.06
                eact=157.9
                tg=249.0
                a1=-52.77
                a2=0.211
                ta=243.0
                b1=8.561
                b2=0.027
                tb=243.0        
            case (2) !levoglucosan/NH4HSO4
                logdwtg0=-45.58
                eact=142.2
                tg=206.5
                a1=-44.96
                a2=0.185
                ta=243.0
                b1=16.57
                b2=-0.063
                tb=243.0        
            case (3) !raffinose
                logdwtg0=-15.76
                eact=80.1
                tg=378.3
                a1=-139.9
                a2=0.347
                ta=273.5
                b1=17.00
                b2=0.00
                tb=273.5       
            case (4) !3-MBTCA
                logdwtg0=-24.86
                eact=64.5
                tg=305.0
                a1=-5.033
                a2=0.015
                ta=300.0
                b1=0.00
                b2=0.00
                tb=300.0        
            case (5) !alpha-pinene
                logdwtg0=-26.60
                eact=65.5
                tg=270.0
                a1=-18.31
                a2=0.063
                ta=273.0
                b1=-10.65
                b2=0.039
                tb=273.0        
            case (6) !sucrose
                logdwtg0=-18.22
                eact=190.3
                tg=335.7
                a1=-16.65
                a2=0.050
                ta=253.0
                b1=-14.65
                b2=0.050
                tb=253.0    
            case (7) !citric acid
                logdwtg0=-29.67
                eact=122.3
                tg=280.1
                a1=-41.00
                a2=0.143
                ta=265.0
                b1=-69.00
                b2=0.280
                tb=255.0
            case (8) !shikimic acid
                logdwtg0=-18.21
                eact=204.9
                tg=326.8
                a1=-16.30
                a2=0.062
                ta=263.0
                b1=16.3
                b2=-0.062
                tb=263.0
        case default
            print*, "selected compound not found"
        end select

        zeta_a_0_aw0=logdwtg0+1.e3*eact/8.314/tg
        zeta_a_aw0=zeta_a_0_aw0-1.e3*eact/8.314/t

        t1=t
        t2=t
        if (t1>ta) then 
            t1=ta 
        end if
        if (t2>tb) then 
            t2=tb 
        end if
        a=(a1+a2*t1)
        b=(b1+b2*t2)
    
        !alpha=exp( (1-molefrac).^2) 
        !alpha=1.
        alpha(:)=exp( (1.-molefrac)**2.*(a+3.*b-4.*b*(1-molefrac))  )

        t0=118
        s=892
        zeta_v_0=log(3.06e-3)
        zeta_v=zeta_v_0-s/(t-t0)

        dwt0=exp(zeta_a_aw0)
        dwt1=exp(zeta_a_aw0+(molefrac)*alpha*(zeta_v-zeta_a_aw0))
    
        ! cm**2 s**-1
        d_coeff=dwt0**(1-molefrac*alpha)*dwt1**(molefrac*alpha)
        ! m**2 s**-1
        d_coeff=d_coeff*1.e-4

    end subroutine Lienhard2015


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Price2014								       !	
    ! http://www.atmos-chem-phys.net/14/3817/2014/ 				       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine Price2014(kp, molefrac, compound, d_coeff)
        implicit none
        real(sp), dimension(kp), intent(in) :: molefrac
        integer(i4b), intent(in) :: kp, compound
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: a, b, c, d

        select case (compound)
            case (1) !sucrose
                a = -20.89
                b = 25.92
                c = -26.97
                d = 13.35
            case (2) !levoglucosan
                a = -18.41
                b = 31.10
                c = -44.43
                d = 23.12
            case (3) !MgSO4
                print*, 'data plotted, but no values present in paper'
            case (4) !raffinose
                a = -17.21
                b = 24.00
                c = -32.50
                d = 17.02
            case default
            print*, "selected compound not found"
        end select
    
        d_coeff = 10.**(a+b*molefrac+c*molefrac**2+d*molefrac**3)

    end subroutine Price2014


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Price2015								       !	
    ! http://pubs.rsc.org/en/content/articlehtml/2015/sc/c5sc00685f [Paper]	       !
    ! http://www.rsc.org/suppdata/c5/sc/c5sc00685f/c5sc00685f1.pdf [Supplementary] !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] t: temperature
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine Price2015(kp, t, molefrac, compound, d_coeff)
        implicit none
        integer(i4b), intent(in) :: kp, compound
        real(sp), intent(in) :: t
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: do_som, do_wat, c, d
        real(sp), dimension(kp) :: alpha, dwt1


        ! alpha-pinene
        if (compound==1) then
            do_som = 10.**( -(7.4+(650./(t-165))) )
            do_wat = 10.**( -(6.514+(387.4/(t-118))) )
                if (t>230) then
                c = -13.+0.043*230
                d = -10.5+0.035*230
                else
                c = -13.+0.043*t
                d = -10.5+0.035*t
                end if
        alpha = exp((1-molefrac)**2.*(c+3.*d-4.*d*(1-molefrac)))
            d_coeff = do_wat**(molefrac*alpha)*do_som**(1-molefrac*alpha)
        else
                print*, "selected compound not found"
        end if

    end subroutine Price2015

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Price2016								       !	
    ! http://xlink.rsc.org/?DOI=C6CP03238A 					       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine Price2016(kp, molefrac, compound, d_coeff)
        implicit none
        real(sp), dimension(kp), intent(in) :: molefrac
        integer(i4b), intent(in) :: kp, compound
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: a, b, c, d

        select case (compound)
            case (1) !sucrose
                a = -30.97
                b = 54.89
                c = -62.34
                d = 29.12
        case (2) !water
                a = -20.89
                b = 25.92
                c = -26.97
                d = 13.35
            case default
            print*, "selected compound not found"
        end select
    
        d_coeff = 10.**(a+b*molefrac+c*molefrac**2.+d*molefrac**3.)

    end subroutine Price2016


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Zobrist2011								       !	
    ! http://pubs.rsc.org/en/content/articlehtml/2011/cp/c0cp01273d 	       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] t: temperature
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine Zobrist2011(kp, t, molefrac, compound, d_coeff)
        implicit none
        integer(i4b), intent (in) :: kp, compound
        real(sp), intent(in) :: t
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: a, b, c, d, e, f, g, t_theta
        real(sp), dimension(kp) :: molefracd, aw, ad, bd, to


        ! sucrose
        if (compound==1) then  
                a = -1.
                b = -0.99721
                c = 0.13599
                d = 0.001688
                e = -0.005151
                f = 0.009607
                g = -0.006142
                t_theta = 298.15		! kelvin, 160<T<313
            molefracd=1-molefrac
                aw = (1.+a*molefracd)/(1.+b*molefracd+c*molefracd**2.) &
                +(t-t_theta)*(d*molefracd+e*molefracd**2.+ &
                f*molefracd**3+g*molefracd**4)
                ad = 7.+0.175*(1.-46.46*(1.-aw))
                bd = 262.867*(1.+10.53*(1.-aw)-0.3*(1.-aw)**2.)
                to = 127.9*(1.+0.4514*(1.-aw)-0.51*(1.-aw)**1.7) 
                d_coeff = 10.**( -(ad+(bd/(t-to))) )
        else
            print*, 'selected compound not found'
        end if

    end subroutine Zobrist2011


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Shiraiwa2013								       !	
    ! http://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp51595h 	       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] kp: number of grid points / shells
    !>@param[in] molefrac: water molefraction
    !>@param[in] t: temperature
    !>@param[in] d_self: self diffusion coefficients
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine Shiraiwa2013(kp, molefrac, d_self, compound, d_coeff)
        implicit none
        integer(i4b), intent(in) :: kp, compound
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(in) :: d_self
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: rho1, rho2, m1, m2, f, z, d12, d11 
        real(sp), dimension (kp) :: v1, v2, volf, d12d, d11d
    

        ! function variables	
        f = 1. 	! [0.65-1]
        z = 16 	! [8 - 16]

        ! alpha-pinene
        if (compound==1) then
    
                rho1 = 1000. 	! kg/m^3
                rho1 = 1. 	! g/cm^3
                rho2 = 858. 	! kg/m^3
                rho2 = rho2/1.e3 ! g/cm^3
                m1 = 18.02 	! g/mol
                m2 = 136.23 	! g/mol
    
                v1 = molefrac*m1/rho1
                v2 = (1.-molefrac)*m2/rho2
    
                volf = v1/(v1+v2)
                !volf = molefrac
    
                d12 = d_self(1)
                d11 = d_self(2)

                d12d = (z*(1.-volf)/2./f-1.)*d12
                d11d = (z*volf/2./f-1.)*d11

                d_coeff = (d12d+d11d+sqrt((d12d+d11d)**2.+2.*(z-2.)*d12*d11))/(z-2)

        else
            print*, 'selected compound not found'

        end if

    end subroutine Shiraiwa2013


    end module diffusion_coefficients

