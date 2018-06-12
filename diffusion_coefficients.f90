module diffusion_coefficients

contains

! read_in_namelist
subroutine read_in_dc_namelist(nmlfile, kp, n_comp, molefrac, t, d_s, param, compound, d_c)
	implicit none
	integer, intent(out) :: kp, n_comp
	real, intent(out) :: t
	real, allocatable, dimension(:), intent(out) :: molefrac
	real, allocatable, dimension (:), intent(out) :: d_c
	real, allocatable, dimension (:), intent(out) :: d_s
	character (len=50), intent(out) :: param, compound
	character (len=11), intent(in) :: nmlfile

	! define namelists
	namelist /dc_setup/ kp, n_comp
	namelist /dc_vars/ molefrac, t, d_s, param, compound

	! read in namelist
	open(10, file=nmlfile, status='old', recl=80, delim='apostrophe')
	read(10, nml=dc_setup)			
		allocate (molefrac(1:kp))
		allocate (d_s(1:n_comp))
		allocate (d_c(1:kp))
	read(10, nml=dc_vars)	

	close(10)

end subroutine read_in_dc_namelist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! d_coeff
subroutine d_coeff(kp, molefrac, t, d_s, param, compound, d_c)
	implicit none
	integer, intent(in) :: kp
	real, intent(in) :: t
	real, allocatable, dimension(:), intent(in) :: molefrac
	real, allocatable, dimension (:), intent(inout) :: d_c
	real, allocatable, dimension (:), intent(in) :: d_s
	character (len=50), intent(in) :: param, compound

	select case (param)
	case ('constant')
	! constant diffusion coefficient
		d_c(:)=d_s(1)
	case ('darken')
	! linear relation between water molefraction and diffusion coefficient
		d_c=d_s(1)*(1-molefrac)+d_s(2)*(molefrac)
	case ('vignes')
	! logarithmic relation between water mole fraction and diffusion coefficient
		d_c=d_s(1)**(1-molefrac)*d_s(2)**(molefrac)
	case('Lienhard2014')
	! http://pubs.rsc.org/en/content/articlehtml/2014/cp/c4cp01939c
		call Lienhard2014(t, molefrac, compound, d_c)
	case('Lienhard2015')
	! http://www.atmos-chem-phys-discuss.net/15/24473/2015/acpd-15-24473-2015.pdf [discussion]
	! http://www.atmos-chem-phys.net/15/13599/2015/acp-15-13599-2015.pdf [final paper]
		call Lienhard2015(t, molefrac, compound, d_c)
	case default
		print*, "selected param not found"
	end select

	print*, "diffusion coefficient = ", d_c

end subroutine d_coeff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lienhard2014
! http://pubs.rsc.org/en/content/articlehtml/2014/cp/c4cp01939c
subroutine Lienhard2014(t, molefrac, compound, d_c)
	implicit none
	real, intent(in) :: t
	real, allocatable, dimension(:), intent(in) :: molefrac
	character (len=50), intent(in) :: compound
	real, allocatable, dimension (:), intent(inout) :: d_c
	real :: d_cit, d_w, c, d
	real, allocatable, dimension(:) :: alpha

	if (compound == 'citric acid') then
		d_cit = 10.**(-15.-(175./(t-208.)))
		d_w = 10.**(-6.514-(387.4/(t-118.)))
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
    		d_c = d_w**(molefrac*alpha)*d_cit**(1.-molefrac*alpha)
	else
		print*, 'selected compound not found'
	end if

end subroutine Lienhard2014


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lienhard2015
! http://www.atmos-chem-phys-discuss.net/15/24473/2015/acpd-15-24473-2015.pdf [discussion]
! http://www.atmos-chem-phys.net/15/13599/2015/acp-15-13599-2015.pdf [final paper]
subroutine Lienhard2015(t, molefrac, compound, d_c)
	implicit none
	real, intent(in) :: t
	real, allocatable, dimension(:), intent(in) :: molefrac
	character (len=50), intent(in) :: compound
	real, allocatable, dimension (:), intent(inout) :: d_c
	real :: logdwtg0, eact, tg, a, a1, a2, b, b1, b2, &
		t0, t1, t2, ta, tb, s, zeta_a_0_aw0, zeta_a_aw0, &
		zeta_v_0, zeta_v, dwt0
	real, allocatable, dimension(:) :: alpha, dwt1

	select case (compound)
	case ('levoglucosan')
        	logdwtg0=-27.06
        	eact=157.9
        	tg=249.0
        	a1=-52.77
        	a2=0.211
        	ta=243.0
        	b1=8.561
        	b2=0.027
        	tb=243.0        
    	case ('levoglucosan/NH4HSO4')
        	logdwtg0=-45.58
        	eact=142.2
        	tg=206.5
        	a1=-44.96
        	a2=0.185
        	ta=243.0
        	b1=16.57
        	b2=-0.063
        	tb=243.0        
    	case ('raffinose')
        	logdwtg0=-15.76
        	eact=80.1
        	tg=378.3
        	a1=-139.9
        	a2=0.347
        	ta=273.5
        	b1=17.00
        	b2=0.00
        	tb=273.5       
    	case ('3-MBTCA')
        	logdwtg0=-24.86
        	eact=64.5
        	tg=305.0
        	a1=-5.033
        	a2=0.015
        	ta=300.0
        	b1=0.00
        	b2=0.00
        	tb=300.0        
    	case ('alpha-pinene')
        	logdwtg0=-26.60
        	eact=65.5
        	tg=270.0
        	a1=-18.31
        	a2=0.063
        	ta=273.0
        	b1=-10.65
        	b2=0.039
        	tb=273.0        
    	case ('sucrose')
        	logdwtg0=-18.22
        	eact=190.3
        	tg=335.7
        	a1=-16.65
        	a2=0.050
        	ta=253.0
        	b1=-14.65
        	b2=0.050
        	tb=253.0    
    	case ('citric acid')
        	logdwtg0=-29.67
        	eact=122.3
        	tg=280.1
        	a1=-41.00
        	a2=0.143
        	ta=265.0
        	b1=-69.00
        	b2=0.280
        	tb=255.0
    	case ('shikimic acid')
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
	alpha=exp( (1.-molefrac)**2.*(a+3.*b-4.*b*(1-molefrac))  )

	t0=118
	s=892
	zeta_v_0=log(3.06e-3)
	zeta_v=zeta_v_0-s/(t-t0)

	dwt0=exp(zeta_a_aw0)
	dwt1=exp(zeta_a_aw0+(molefrac)*alpha*(zeta_v-zeta_a_aw0))
	
	! cm**2 s**-1
	d_c=dwt0**(1-molefrac*alpha)*dwt1**(molefrac*alpha)
	! m**2 s**-1
	d_c=d_c*1.e-4

end subroutine Lienhard2015

end module diffusion_coefficients


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program main

use diffusion_coefficients
implicit none

integer :: kp, n_comp
real :: t 
real, allocatable, dimension(:) :: molefrac
real, allocatable, dimension (:) :: d_s
real, allocatable, dimension (:) :: d_c
character (len=50) :: param, compound


call read_in_dc_namelist('namelist.in', kp, n_comp, molefrac, t, d_s, param, compound, d_c)


call d_coeff(kp, molefrac, t, d_s, param, compound, d_c)


end program main


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

