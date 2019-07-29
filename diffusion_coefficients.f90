!>@author
!>Kathryn Fowler, The University of Manchester
!>@brief
!>module for SOA diffusion coefficient parameterisations


    module diffusion_coefficients
    use numerics_type

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
        real(wp), intent(out) :: t
        real(wp), allocatable, dimension(:), intent(out) :: molefrac
        real(wp), allocatable, dimension (:), intent(out) :: d_coeff
        real(wp), allocatable, dimension (:), intent(out) :: d_self
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
    !>@param[in] n_comps: number of compositions (including water)
    !>@param[in] molefrac: water molefraction
    !>@param[in] t: temperature
    !>@param[in] d_self: self diffusion coefficients
    !>@param[in] param: parameterisation type
    !>@param[in] compound: organic component of aerosol
    !>@param[inout] d_coeff: mutual diffusion coefficient
    subroutine diffusion_coefficient(kp, n_comps, molefrac, t, &
                                d_self, param, compound, d_coeff)
        implicit none
        integer(i4b), intent(in) :: kp, n_comps, param, compound
        real(wp), intent(in) :: t
        real(wp), dimension(kp), intent(in) :: molefrac
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp), dimension (n_comps), intent(in) :: d_self


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
                stop

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
        real(wp), intent(in) :: t
        real(wp), dimension(kp), intent(in) :: molefrac
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp) :: d_coeffit, d_molefrac, c, d
        real(wp), dimension(kp) :: alpha


        ! citric acid
        if (compound == 1) then
            d_coeffit = 10._wp**(-15._wp-(175._wp/(t-208._wp)))
            d_molefrac = 10._wp**(-6.514_wp-(387.4_wp/(t-118._wp)))
            if (t>265._wp) then
                c = -41._wp+0.143_wp*265._wp
            else
                c = -41._wp+0.143_wp*t
            end if
            if (t>255._wp) then
                d = -69._wp+0.28_wp*255._wp
                else
                    d = -69._wp+0.28_wp*t
            end if
            alpha = exp((1._wp-molefrac)**2._wp*(c+3._wp*d-4._wp*d*(1._wp-molefrac)))
                d_coeff = d_molefrac**(molefrac*alpha)*d_coeffit**(1._wp-molefrac*alpha)
        else
            print*, 'selected compound not found'
            stop
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
        real(wp), intent(in) :: t
        real(wp), dimension(kp), intent(in) :: molefrac
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp) :: logdwtg0, eact, tg, a, a1, a2, b, b1, b2, &
            t0, t1, t2, ta, tb, s, zeta_a_0_aw0, zeta_a_aw0, &
            zeta_v_0, zeta_v, dwt0
        real(wp), dimension(kp) :: alpha, dwt1

    

        select case (compound)
    
        case (1) !levoglucosan
                logdwtg0=-27.06_wp
                eact=157.9_wp
                tg=249.0_wp
                a1=-52.77_wp
                a2=0.211_wp
                ta=243.0_wp
                b1=8.561_wp
                b2=0.027_wp
                tb=243.0_wp     
            case (2) !levoglucosan/NH4HSO4
                logdwtg0=-45.58_wp
                eact=142.2_wp
                tg=206.5_wp
                a1=-44.96_wp
                a2=0.185_wp
                ta=243.0_wp
                b1=16.57_wp
                b2=-0.063_wp
                tb=243.0_wp    
            case (3) !raffinose
                logdwtg0=-15.76_wp
                eact=80.1_wp
                tg=378.3_wp
                a1=-139.9_wp
                a2=0.347_wp
                ta=273.5_wp
                b1=17.00_wp
                b2=0.00_wp
                tb=273.5_wp     
            case (4) !3-MBTCA
                logdwtg0=-24.86_wp
                eact=64.5_wp
                tg=305.0_wp
                a1=-5.033_wp
                a2=0.015_wp
                ta=300.0_wp
                b1=0.00_wp
                b2=0.00_wp
                tb=300.0_wp      
            case (5) !alpha-pinene
                logdwtg0=-26.60_wp
                eact=65.5_wp
                tg=270.0_wp
                a1=-18.31_wp
                a2=0.063_wp
                ta=273.0_wp
                b1=-10.65_wp
                b2=0.039_wp
                tb=273.0_wp     
            case (6) !sucrose
                logdwtg0=-18.22_wp
                eact=190.3_wp
                tg=335.7_wp
                a1=-16.65_wp
                a2=0.050_wp
                ta=253.0_wp
                b1=-14.65_wp
                b2=0.050_wp
                tb=253.0_wp 
            case (7) !citric acid
                logdwtg0=-29.67_wp
                eact=122.3_wp
                tg=280.1_wp
                a1=-41.00_wp
                a2=0.143_wp
                ta=265.0_wp
                b1=-69.00_wp
                b2=0.280_wp
                tb=255.0_wp
            case (8) !shikimic acid
                logdwtg0=-18.21_wp
                eact=204.9_wp
                tg=326.8_wp
                a1=-16.30_wp
                a2=0.062_wp
                ta=263.0_wp
                b1=16.3_wp
                b2=-0.062_wp
                tb=263.0_wp
        case default
            print*, "selected compound not found"
            stop
        end select

        zeta_a_0_aw0=logdwtg0+1.e3_wp*eact/8.314_wp/tg
        zeta_a_aw0=zeta_a_0_aw0-1.e3_wp*eact/8.314_wp/t

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
        alpha(:)=exp( (1._wp-molefrac)**2._wp*(a+3._wp*b-4._wp*b*(1._wp-molefrac))  )

        t0=118._wp
        s=892._wp
        zeta_v_0=log(3.06e-3_wp)
        zeta_v=zeta_v_0-s/(t-t0)

        dwt0=exp(zeta_a_aw0)
        dwt1=exp(zeta_a_aw0+(molefrac)*alpha*(zeta_v-zeta_a_aw0))
    
        ! cm**2 s**-1
        d_coeff=dwt0**(1._wp-molefrac*alpha)*dwt1**(molefrac*alpha)
        ! m**2 s**-1
        d_coeff=d_coeff*1.e-4_wp

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
        real(wp), dimension(kp), intent(in) :: molefrac
        integer(i4b), intent(in) :: kp, compound
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp) :: a, b, c, d

        select case (compound)
            case (1) !sucrose
                a = -20.89_wp
                b = 25.92_wp
                c = -26.97_wp
                d = 13.35_wp
            case (2) !levoglucosan
                a = -18.41_wp
                b = 31.10_wp
                c = -44.43_wp
                d = 23.12_wp
            case (3) !MgSO4
                print*, 'data plotted, but no values present in paper'
                stop
            case (4) !raffinose
                a = -17.21_wp
                b = 24.00_wp
                c = -32.50_wp
                d = 17.02_wp
            case default
                print*, "selected compound not found"
                stop
        end select
    
        d_coeff = 10._wp**(a+b*molefrac+c*molefrac**2+d*molefrac**3)

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
        real(wp), intent(in) :: t
        real(wp), dimension(kp), intent(in) :: molefrac
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp) :: do_som, do_wat, c, d
        real(wp), dimension(kp) :: alpha, dwt1


        ! alpha-pinene
        if (compound==1) then
            do_som = 10._wp**( -(7.4_wp+(650._wp/(t-165._wp))) )
            do_wat = 10._wp**( -(6.514_wp+(387.4_wp/(t-118._wp))) )
                if (t>230._wp) then
                c = -13._wp+0.043_wp*230._wp
                d = -10.5_wp+0.035_wp*230._wp
                else
                c = -13._wp+0.043_wp*t
                d = -10.5_wp+0.035_wp*t
                end if
        alpha = exp((1._wp-molefrac)**2.*(c+3._wp*d-4._wp*d*(1._wp-molefrac)))
            d_coeff = do_wat**(molefrac*alpha)*do_som**(1._wp-molefrac*alpha)
        else
                print*, "selected compound not found"
                stop
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
        real(wp), dimension(kp), intent(in) :: molefrac
        integer(i4b), intent(in) :: kp, compound
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp) :: a, b, c, d

        select case (compound)
            case (1) !sucrose
                a = -30.97_wp
                b = 54.89_wp
                c = -62.34_wp
                d = 29.12_wp
        case (2) !water
                a = -20.89_wp
                b = 25.92_wp
                c = -26.97_wp
                d = 13.35_wp
        case default
            print*, "selected compound not found"
            stop
        end select
    
        d_coeff = 10._wp**(a+b*molefrac+c*molefrac**2.+d*molefrac**3.)

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
        real(wp), intent(in) :: t
        real(wp), dimension(kp), intent(in) :: molefrac
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp) :: a, b, c, d, e, f, g, t_theta
        real(wp), dimension(kp) :: molefracd, aw, ad, bd, to


        ! sucrose
        if (compound==1) then  
                a = -1._wp
                b = -0.99721_wp
                c = 0.13599_wp
                d = 0.001688_wp
                e = -0.005151_wp
                f = 0.009607_wp
                g = -0.006142_wp
                t_theta = 298.15_wp		! kelvin, 160<T<313
            molefracd=1-molefrac
                aw = (1._wp+a*molefracd)/(1._wp+b*molefracd+c*molefracd**2._wp) &
                +(t-t_theta)*(d*molefracd+e*molefracd**2._wp+ &
                f*molefracd**3._wp+g*molefracd**4._wp)
                ad = 7._wp+0.175_wp*(1._wp-46.46_wp*(1.-aw))
                bd = 262.867_wp*(1.+10.53_wp*(1._wp-aw)-0.3_wp*(1._wp-aw)**2._wp)
                to = 127.9_wp*(1._wp+0.4514_wp*(1._wp-aw)-0.51_wp*(1._wp-aw)**1.7_wp) 
                d_coeff = 10._wp**( -(ad+(bd/(t-to))) )
        else
            print*, 'selected compound not found'
            stop
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
        real(wp), dimension(kp), intent(in) :: molefrac
        real(wp), dimension (kp), intent(in) :: d_self
        real(wp), dimension (kp), intent(inout) :: d_coeff
        real(wp) :: rho1, rho2, m1, m2, f, z, d12, d11 
        real(wp), dimension (kp) :: v1, v2, volf, d12d, d11d
    

        ! function variables	
        f = 1._wp 	! [0.65-1]
        z = 16._wp 	! [8 - 16]

        ! alpha-pinene
        if (compound==1) then
    
                rho1 = 1000._wp 	! kg/m^3
                rho1 = 1._wp 	! g/cm^3
                rho2 = 858._wp 	! kg/m^3
                rho2 = rho2/1.e3_wp ! g/cm^3
                m1 = 18.02_wp 	! g/mol
                m2 = 136.23_wp 	! g/mol
    
                v1 = molefrac*m1/rho1
                v2 = (1._wp-molefrac)*m2/rho2
    
                volf = v1/(v1+v2)
                !volf = molefrac
    
                d12 = d_self(1)
                d11 = d_self(2)

                d12d = (z*(1._wp-volf)/2._wp/f-1._wp)*d12
                d11d = (z*volf/2._wp/f-1._wp)*d11

                d_coeff = (d12d+d11d+sqrt((d12d+d11d)**2._wp+ &
                    2._wp*(z-2._wp)*d12*d11))/(z-2._wp)

        else
            print*, 'selected compound not found'
            stop
        end if

    end subroutine Shiraiwa2013


    end module diffusion_coefficients


