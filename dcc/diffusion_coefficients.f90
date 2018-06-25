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
        real(sp), intent(in) :: t
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: d_coeffit, d_molefrac, c, d
        real(sp), dimension(kp) :: alpha


        ! citric acid
        if (compound == 1) then
            d_coeffit = 10._sp**(-15._sp-(175._sp/(t-208._sp)))
            d_molefrac = 10._sp**(-6.514_sp-(387.4_sp/(t-118._sp)))
            if (t>265._sp) then
                c = -41._sp+0.143_sp*265._sp
            else
                c = -41._sp+0.143_sp*t
            end if
            if (t>255._sp) then
                d = -69._sp+0.28_sp*255._sp
                else
                    d = -69._sp+0.28_sp*t
            end if
            alpha = exp((1._sp-molefrac)**2._sp*(c+3.*d-4._sp*d*(1._sp-molefrac)))
                d_coeff = d_molefrac**(molefrac*alpha)*d_coeffit**(1._sp-molefrac*alpha)
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
        real(sp), intent(in) :: t
        real(sp), dimension(kp), intent(in) :: molefrac
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: logdwtg0, eact, tg, a, a1, a2, b, b1, b2, &
            t0, t1, t2, ta, tb, s, zeta_a_0_aw0, zeta_a_aw0, &
            zeta_v_0, zeta_v, dwt0
        real(sp), dimension(kp) :: alpha, dwt1

    

        select case (compound)
    
        case (1) !levoglucosan
                logdwtg0=-27.06_sp
                eact=157.9_sp
                tg=249.0_sp
                a1=-52.77_sp
                a2=0.211_sp
                ta=243.0_sp
                b1=8.561_sp
                b2=0.027_sp
                tb=243.0_sp     
            case (2) !levoglucosan/NH4HSO4
                logdwtg0=-45.58_sp
                eact=142.2_sp
                tg=206.5_sp
                a1=-44.96_sp
                a2=0.185_sp
                ta=243.0_sp
                b1=16.57_sp
                b2=-0.063_sp
                tb=243.0_sp    
            case (3) !raffinose
                logdwtg0=-15.76_sp
                eact=80.1_sp
                tg=378.3_sp
                a1=-139.9_sp
                a2=0.347_sp
                ta=273.5_sp
                b1=17.00_sp
                b2=0.00_sp
                tb=273.5_sp     
            case (4) !3-MBTCA
                logdwtg0=-24.86_sp
                eact=64.5_sp
                tg=305.0_sp
                a1=-5.033_sp
                a2=0.015_sp
                ta=300.0_sp
                b1=0.00_sp
                b2=0.00_sp
                tb=300.0_sp      
            case (5) !alpha-pinene
                logdwtg0=-26.60_sp
                eact=65.5_sp
                tg=270.0_sp
                a1=-18.31_sp
                a2=0.063_sp
                ta=273.0_sp
                b1=-10.65_sp
                b2=0.039_sp
                tb=273.0_sp     
            case (6) !sucrose
                logdwtg0=-18.22_sp
                eact=190.3_sp
                tg=335.7_sp
                a1=-16.65_sp
                a2=0.050_sp
                ta=253.0_sp
                b1=-14.65_sp
                b2=0.050_sp
                tb=253.0_sp 
            case (7) !citric acid
                logdwtg0=-29.67_sp
                eact=122.3_sp
                tg=280.1_sp
                a1=-41.00_sp
                a2=0.143_sp
                ta=265.0_sp
                b1=-69.00_sp
                b2=0.280_sp
                tb=255.0_sp
            case (8) !shikimic acid
                logdwtg0=-18.21_sp
                eact=204.9_sp
                tg=326.8_sp
                a1=-16.30_sp
                a2=0.062_sp
                ta=263.0_sp
                b1=16.3_sp
                b2=-0.062_sp
                tb=263.0_sp
        case default
            print*, "selected compound not found"
            stop
        end select

        zeta_a_0_aw0=logdwtg0+1.e3_sp*eact/8.314_sp/tg
        zeta_a_aw0=zeta_a_0_aw0-1.e3_sp*eact/8.314_sp/t

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
        alpha(:)=exp( (1._sp-molefrac)**2._sp*(a+3._sp*b-4._sp*b*(1._sp-molefrac))  )

        t0=118._sp
        s=892._sp
        zeta_v_0=log(3.06e-3_sp)
        zeta_v=zeta_v_0-s/(t-t0)

        dwt0=exp(zeta_a_aw0)
        dwt1=exp(zeta_a_aw0+(molefrac)*alpha*(zeta_v-zeta_a_aw0))
    
        ! cm**2 s**-1
        d_coeff=dwt0**(1._sp-molefrac*alpha)*dwt1**(molefrac*alpha)
        ! m**2 s**-1
        d_coeff=d_coeff*1.e-4_sp

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
                a = -20.89_sp
                b = 25.92_sp
                c = -26.97_sp
                d = 13.35_sp
            case (2) !levoglucosan
                a = -18.41_sp
                b = 31.10_sp
                c = -44.43_sp
                d = 23.12_sp
            case (3) !MgSO4
                print*, 'data plotted, but no values present in paper'
                stop
            case (4) !raffinose
                a = -17.21_sp
                b = 24.00_sp
                c = -32.50_sp
                d = 17.02_sp
            case default
                print*, "selected compound not found"
                stop
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
            do_som = 10._sp**( -(7.4_sp+(650._sp/(t-165._sp))) )
            do_wat = 10._sp**( -(6.514_sp+(387.4_sp/(t-118._sp))) )
                if (t>230._sp) then
                c = -13._sp+0.043_sp*230._sp
                d = -10.5_sp+0.035_sp*230._sp
                else
                c = -13._sp+0.043_sp*t
                d = -10.5_sp+0.035_sp*t
                end if
        alpha = exp((1._sp-molefrac)**2.*(c+3._sp*d-4._sp*d*(1._sp-molefrac)))
            d_coeff = do_wat**(molefrac*alpha)*do_som**(1._sp-molefrac*alpha)
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
        real(sp), dimension(kp), intent(in) :: molefrac
        integer(i4b), intent(in) :: kp, compound
        real(sp), dimension (kp), intent(inout) :: d_coeff
        real(sp) :: a, b, c, d

        select case (compound)
            case (1) !sucrose
                a = -30.97_sp
                b = 54.89_sp
                c = -62.34_sp
                d = 29.12_sp
        case (2) !water
                a = -20.89_sp
                b = 25.92_sp
                c = -26.97_sp
                d = 13.35_sp
        case default
            print*, "selected compound not found"
            stop
        end select
    
        d_coeff = 10._sp**(a+b*molefrac+c*molefrac**2.+d*molefrac**3.)

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
                a = -1._sp
                b = -0.99721_sp
                c = 0.13599_sp
                d = 0.001688_sp
                e = -0.005151_sp
                f = 0.009607_sp
                g = -0.006142_sp
                t_theta = 298.15_sp		! kelvin, 160<T<313
            molefracd=1-molefrac
                aw = (1._sp+a*molefracd)/(1._sp+b*molefracd+c*molefracd**2._sp) &
                +(t-t_theta)*(d*molefracd+e*molefracd**2._sp+ &
                f*molefracd**3._sp+g*molefracd**4._sp)
                ad = 7._sp+0.175_sp*(1._sp-46.46_sp*(1.-aw))
                bd = 262.867_sp*(1.+10.53_sp*(1._sp-aw)-0.3_sp*(1._sp-aw)**2._sp)
                to = 127.9_sp*(1._sp+0.4514_sp*(1._sp-aw)-0.51_sp*(1._sp-aw)**1.7_sp) 
                d_coeff = 10._sp**( -(ad+(bd/(t-to))) )
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
    
                rho1 = 1000._sp 	! kg/m^3
                rho1 = 1._sp 	! g/cm^3
                rho2 = 858._sp 	! kg/m^3
                rho2 = rho2/1.e3_sp ! g/cm^3
                m1 = 18.02_sp 	! g/mol
                m2 = 136.23_sp 	! g/mol
    
                v1 = molefrac*m1/rho1
                v2 = (1._sp-molefrac)*m2/rho2
    
                volf = v1/(v1+v2)
                !volf = molefrac
    
                d12 = d_self(1)
                d11 = d_self(2)

                d12d = (z*(1._sp-volf)/2._sp/f-1._sp)*d12
                d11d = (z*volf/2._sp/f-1._sp)*d11

                d_coeff = (d12d+d11d+sqrt((d12d+d11d)**2._sp+ &
                    2._sp*(z-2._sp)*d12*d11))/(z-2._sp)

        else
            print*, 'selected compound not found'
            stop
        end if

    end subroutine Shiraiwa2013


    end module diffusion_coefficients


