	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>module for bin microphysics with diffusion model
	module bdm
    use numerics_type
    use numerics, only : find_pos, poly_int
    use diffusion, only : grid, backward_euler, move_boundary
    use bmm
    use diffusion_coefficients 
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the bin diffusion model

    implicit none


        character (len=200) :: outputfile='output', bin_model_file='', &
                            diffusion_model_file='', diffusion_coeff_file=''
                            
        type(grid), allocatable, dimension(:) :: grida
        
        real(wp), allocatable, dimension(:) :: nwo, jw, nw,ns, aw
        real(wp), allocatable, dimension(:,:) :: nso
        integer(i4b) :: koehler_shell_flag
        integer(i4b) :: diffusion_type
        
        ! for calculating diffusion coefficients
        real(wp), allocatable, dimension (:) :: d_self
        integer(i4b) :: n_compsw, & ! number of compositions (including water)
                        param, & ! type of parameterisation for diffusion coefficients
                        compound ! compound for diffusion coefficients


	private 
	public :: read_in_bdm_namelist, initialise_bdm_arrays, bdm_driver, grida

	contains	
	
		

	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! read in the namelist                                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the bin diffusion model
	!>@param[in] nmlfile
	subroutine read_in_bdm_namelist(nmlfile)
	    use bmm
	    use diffusion, only : nmd
	    use diffusion_coefficients
		implicit none
        character (len=200), intent(in) :: nmlfile
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ outputfile, bin_model_file, diffusion_model_file, &
                            diffusion_coeff_file, diffusion_type, &
                            koehler_shell_flag
        namelist /run_vars/ nmd
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelist for model - pointing to bin and diffusion models    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! bin-microphysics namelist                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call read_in_bmm_namelist(bin_model_file)
        call validate_bdm_bmm_configuration()
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! diffusion namelist                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(8,file=diffusion_model_file,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! diffusion coefficients namelist                                      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call read_in_dc_namelist2(diffusion_coeff_file, n_compsw, &
            d_self, param, compound)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end subroutine read_in_bdm_namelist

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! validate the deliberately restricted BDM/BMM coupling                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine validate_bdm_bmm_configuration()
        use bmm
        implicit none

        if (bin_scheme_flag /= BIN_FULL_MOVING) &
            error stop 'BDM requires BMM bin_scheme_flag=0 (full_moving)'
        if (sce_flag /= 0_i4b) &
            error stop 'BDM requires sce_flag=0 (no aggregation/collisions)'
        if (.not.adiabatic_prof) &
            error stop 'BDM does not support BMM entrainment; set adiabatic_prof=.true.'
        if (n_comps /= 1_i4b) &
            error stop 'BDM currently supports exactly one soluble aerosol component'
        if (any(nu_core1(1:n_comps) <= 0._wp)) &
            error stop 'BDM requires a soluble aerosol component with nu_core1 > 0'
        if (sv_flag /= 0_i4b) &
            error stop 'BDM does not support semivolatile aerosol; set sv_flag=0'
        if (kappa_flag /= 0_i4b) &
            error stop 'BDM radial activity requires molecular Koehler; set kappa_flag=0'
        if (n_inp_classes /= 0_i4b) &
            error stop 'BDM supports homogeneous nucleation only; set n_inp_classes=0'

        if (ice_flag == 1_i4b) then
            if (.not.ice_nucleation_mech(INUC_KOOP) .or. &
                any(ice_nucleation_mech(2:N_INUC_MECH))) &
                error stop 'BDM requires Koop-only homogeneous nucleation'
            if (mode1_flag .or. mode2_flag .or. hm_flag .or. break_flag /= 0_i4b) &
                error stop 'BDM does not support secondary-ice/breakup mechanisms'
        endif

        if (entrain_period /= 0_i4b .or. entrain_aerosol .or. release_aerosol) &
            error stop 'BDM does not support entrainment/aerosol exchange'
        if (abs(ent_rate) > tiny(1._wp)) &
            error stop 'BDM does not support entrainment; set ent_rate=0'
        if (chamber_forcing_active .or. chamber_bl_mix /= 0_i4b .or. &
            chamber_fan_loss /= 0_i4b .or. chamber_wall_loss /= 0_i4b) &
            error stop 'BDM does not support BMM chamber forcing/loss options'

    end subroutine validate_bdm_bmm_configuration

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! read in dcc namelist							                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>@author
    !>Kathryn Fowler and Paul J. Connolly, The University of Manchester
    !>@brief
    !>read in the data from the namelists for diffusion coefficients
    !>@param[in] nmlfile: namelist file
    !>@param[out] n_comp: number of components
    !>@param[out] d_self: self diffusion coefficients
    !>@param[out] param: parameterisation type
    !>@param[out] compound: organic component of aerosol
    subroutine read_in_dc_namelist2(nmlfile, n_comp, &
                d_self, param, compound)
        implicit none
        integer(i4b), intent(out) :: n_comp, param, compound
        real(wp), allocatable, dimension (:), intent(out) :: d_self
        character (len=*), intent(in) :: nmlfile

        ! define namelists
        namelist /dc_setup/ n_comp
        namelist /dc_vars/ d_self, param, compound

        ! read in namelist
        open(10, file=nmlfile, status='old', recl=80, delim='apostrophe')
        read(10, nml=dc_setup)			
        allocate (d_self(1:n_comp))
        read(10, nml=dc_vars)	

        close(10)

    end subroutine read_in_dc_namelist2
    
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialise arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>sets up the arrays
    subroutine initialise_bdm_arrays()
        use numerics_type
        use numerics, only : find_pos, poly_int, vode_integrate, zeroin, fmin
        use bmm
        use diffusion, only : allocate_and_set_diff, nmd, gridd

        implicit none
        integer(i4b) :: AllocateStatus,i
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! bin-microphysics arrays                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call initialise_bmm_arrays(psurf, tsurf, q_read, theta_read, rh_read, z_read, &
                    time_chamber, press_chamber, temp_chamber, qtot_chamber, &
                    runtime, dt, zinit, tpert, use_prof_for_tprh, &
                    winit, tinit, pinit, &
                    rhinit, radinit, bubble_flag, &
                    microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, adiabatic_prof, vert_ent, z_ctop, &
                    ent_rate, n_levels_s, n_levels_c, &
                    alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_comps, &
                    n_aer1,d_aer1,sig_aer1,dmina,dmaxa,mass_frac_aer1,molw_core1, &
                    density_core1, nu_core1, kappa_core1, org_content1, molw_org1, &
                    kappa_org1, density_org1, delta_h_vap1,nu_org1, log_c_star1,sce_flag)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! work out the diameter, including the water...                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call wetdiam(parcel1%mbin(:,n_comps+1),parcel1%mbin,parcel1%rhobin, &
                    parcel1%n_bin_modew,parcel1%dw) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays for diffusion model         		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( grida(1:parcel1%n_bin_modew), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nwo(1:parcel1%n_bin_modew), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nso(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

		allocate( jw(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nw(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( ns(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( aw(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

        
        do i=1,parcel1%n_bin_modew
            call allocate_and_set_diff(nmd%kp,nmd%dt,nmd%runtime, &
!                parcel1%dw(i)/2._wp, nmd%rad_min, nmd%rad_max, nmd%t,nmd%p, parcel1%rh, &
                parcel1%dw(i)/2._wp, nmd%rad_min, nmd%rad_max, nmd%t,nmd%p, &
                parcel1%mbin(i,n_comps+1)/molw_water/ &
                    (parcel1%mbin(i,n_comps+1)/molw_water + &
                     parcel1%mbin(i,1)/parcel1%molwbin(i,1)*parcel1%nubin(i,1)), &
                parcel1%molwbin(i,1), &
                parcel1%rhobin(i,1), nu_core1(1), nmd%d_coeff, &
                grida(i)%kp,grida(i)%kp_cur, &
                grida(i)%ntim,grida(i)%dt, grida(i)%rad, &
                grida(i)%rad_min,grida(i)%rad_max, &
                grida(i)%t,grida(i)%p,&
                grida(i)%rh, &
                grida(i)%mwsol,grida(i)%rhosol,grida(i)%d_coeff,grida(i)%r, &
                grida(i)%r_old,grida(i)%r05,grida(i)%r05_old, &
                grida(i)%u,grida(i)%d,grida(i)%d05,grida(i)%dr,&
                grida(i)%dr_old,grida(i)%dr05,grida(i)%dr05_old,grida(i)%vol, &
                grida(i)%vol_old,grida(i)%c,grida(i)%cold)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine initialise_bdm_arrays
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! driver for bdm                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>driver for the bin diffusion model
    subroutine bdm_driver()
    use numerics_type
    use bmm, only : parcel, parcel1
    implicit none
    integer(i4b) :: i, j,nt
    logical :: new_file=.true.
    real(wp) :: flux,deltaV, radius, radiusold, t_tstep
    
    
    ! Keep the BMM/DVODE tolerances selected by initialise_bmm_arrays.
    ! BDM used to tighten these to rtol=1e-6 and water atol=1e-30 kg,
    ! which is unnecessarily severe for the parcel-water ODE and can force
    ! DVODE to take O(1e4) internal steps near evaporation/reactivation.
    
    
    
    nt=ceiling(runtime / real(dt,kind=wp))
    do i=1,nt
        t_tstep=parcel1%y(parcel1%ite) ! use a constant temperature for the diffusion 
                                       ! over the time-step
    
        ! output to file
        call output(io1%new_file,outputfile)
        call outputdiff(new_file,outputfile)
        
        
        
        ! store old aerosol state
        do j=1,parcel1%n_bin_modew
            grida(j)%cold=grida(j)%c
            grida(j)%kp_cur_old=grida(j)%kp_cur
            grida(j)%rad_old=grida(j)%rad
            grida(j)%r_old=grida(j)%r
            grida(j)%r05_old=grida(j)%r05
            grida(j)%dr_old=grida(j)%dr
            grida(j)%dr05_old=grida(j)%dr05
            grida(j)%vol_old=grida(j)%vol
            grida(j)%t=t_tstep ! will use this for the temperature

        enddo
        
        ! one time-step of model
        call bin_microphysics(fparcelwarmdiff,fparcelcold,icenucleation, &
                              noncollisional_iceformation_diff)
        
        
        ! diffusion out-side of solver
!         do j=1,n_bins
!             grida(j)%c=grida(j)%cold
!             grida(j)%kp_cur=grida(j)%kp_cur_old
!             grida(j)%rad=grida(j)%rad_old
!             grida(j)%r=grida(j)%r_old
!             grida(j)%r05=grida(j)%r05_old
!             grida(j)%dr=grida(j)%dr_old
!             grida(j)%dr05=grida(j)%dr05_old
!             grida(j)%vol=grida(j)%vol_old
!             
!             
! !             deltaV=max(parcel1%y(j)-parcel1%yold(j),-parcel1%y(j))/rhow
!             deltaV=(parcel1%y(j)-&
!                 sum(grida(j)%c(1:grida(j)%kp_cur,1)* &
!                 grida(j)%vol(1:grida(j)%kp_cur))*molw_water) /rhow
!         
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! shift radii and calculate the velocity of boundaries                       !
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  			call move_boundary(grida(j)%kp,grida(j)%kp_cur,parcel1%dt, &
!  			    radiusold,radius,grida(j)%r,grida(j)%r05,grida(j)%dr,grida(j)%dr05, &
!  			    grida(j)%vol,grida(j)%u,grida(j)%c,flux, &
!  			    grida(j)%rad_min,grida(j)%rad_max, grida(j)%mwsol, grida(j)%rhosol, &
!  			    deltaV)
!             radiusold=radius
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
! 
! 
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! Set diffusion coefficients - inc. zero at boundary                         !
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			select case (diffusion_type)
! 			    case(0)
! 			        grida(j)%d05(:)=grida(j)%d_coeff
! 			    case(1)
! 			        call diffusion_coefficient(grida(j)%kp_cur, &
! 			                grida(j)%c(1:grida(j)%kp_cur,1) / &
! 			                    sum(grida(j)%c(1:grida(j)%kp_cur,:),2), &
! 			                grida(j)%t, d_self, param, &
! 			                compound, grida(j)%d05(1:grida(j)%kp_cur))
! 			    case default
! 			        print *,'error diffusion type'
! 			        stop
! 			end select
! 			grida(j)%d05(grida(j)%kp_cur:grida(j)%kp) = 0._wp
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
! 
! 
! 
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			! solve diffusion equation                                                   !
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			call backward_euler(grida(j)%kp,grida(j)%kp_cur,parcel1%dt, &
! 			    grida(j)%r,grida(j)%r05,grida(j)%u,grida(j)%d,grida(j)%d05,&
! 			    grida(j)%dr,grida(j)%dr05,grida(j)%c,grida(j)%cold,flux)
! 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         
! !             parcel1%y(j)=sum(grida(j)%vol(1:grida(j)%kp_cur) * &
! !                     grida(j)%c(1:grida(j)%kp_cur,1) )*molw_water
!         enddo
        
        ! check there are no negative values
!         where(parcel1%y(1:parcel1%n_bin_modew).le.0.e1_wp)
!             parcel1%y(1:parcel1%n_bin_modew)=1.e-22_wp
!         end where

        ! The ODE solution is authoritative.  Never resurrect water from grida:
        ! grida is mutated during provisional DVODE RHS evaluations and is not
        ! itself the accepted ODE state.  Keep the duplicate BMM warm state
        ! explicitly synchronised; this is redundant with current full_moving
        ! BMM, but also keeps the older BMM subtree consistent.
        where (parcel1%y(1:parcel1%n_bin_modew) < 0._wp)
            parcel1%y(1:parcel1%n_bin_modew)=0._wp
        end where
        parcel1%mbin(1:parcel1%n_bin_modew,n_comps+1)= &
            parcel1%y(1:parcel1%n_bin_modew)
        if (ice_flag.eq.1) then
            parcel1%moments(1:parcel1%n_bin_modew,n_comps+4)= &
                parcel1%npart(1:parcel1%n_bin_modew)* &
                parcel1%y(1:parcel1%n_bin_modew)
            parcel1%moments(1:parcel1%n_bin_modew,n_comps+5)= &
                parcel1%npart(1:parcel1%n_bin_modew)* &
                parcel1%y(1:parcel1%n_bin_modew)
        endif

        ! Rebuild the radial BDM state at the accepted end-of-step water masses.
        call finalise_bdm_radial_state(parcel1%dt)
        
        ! break-out if flag has been set 
        if(parcel1%break_flag) exit
    enddo
    ! output to file
    call output(io1%new_file,outputfile)
    call outputdiff(new_file,outputfile)
    
    
    
    
    end subroutine bdm_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Reconstruct the radial BDM state at the accepted parcel water masses.         !
    ! DVODE may evaluate the RHS at provisional states and may interpolate to tout, !
    ! so grida must not be assumed to contain the accepted end-of-step state.        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine finalise_bdm_radial_state(dt_step)
        use numerics_type
        implicit none
        real(wp), intent(in) :: dt_step
        integer(i4b) :: i
        real(wp) :: deltaV, flux, radius, radiusold, water_old

        do i=1,parcel1%n_bin_modew
            grida(i)%c=grida(i)%cold
            grida(i)%kp_cur=grida(i)%kp_cur_old
            grida(i)%rad=grida(i)%rad_old
            grida(i)%r=grida(i)%r_old
            grida(i)%r05=grida(i)%r05_old
            grida(i)%dr=grida(i)%dr_old
            grida(i)%dr05=grida(i)%dr05_old
            grida(i)%vol=grida(i)%vol_old

            water_old=sum(grida(i)%c(1:grida(i)%kp_cur,1)* &
                          grida(i)%vol(1:grida(i)%kp_cur))*molw_water
            deltaV=(max(parcel1%y(i),0._wp)-water_old)/rhow

            flux=0._wp
            radius=grida(i)%rad
            radiusold=radius
            call move_boundary(grida(i)%kp,grida(i)%kp_cur,dt_step, &
                radiusold,radius,grida(i)%r,grida(i)%r05,grida(i)%dr,grida(i)%dr05, &
                grida(i)%vol,grida(i)%u,grida(i)%c,flux, &
                grida(i)%rad_min,grida(i)%rad_max,grida(i)%mwsol,grida(i)%rhosol, &
                deltaV)
            grida(i)%rad=radius

            select case (diffusion_type)
            case(0)
                grida(i)%d05(:)=grida(i)%d_coeff
                grida(i)%d05(grida(i)%kp_cur:grida(i)%kp+1)=0._wp
            case(1)
                call diffusion_coefficient(grida(i)%kp_cur,n_comps+1, &
                    grida(i)%c(1:grida(i)%kp_cur,1) / &
                    max(sum(grida(i)%c(1:grida(i)%kp_cur,:),2),tiny(1._wp)), &
                    grida(i)%t,d_self,param,compound, &
                    grida(i)%d05(1:grida(i)%kp_cur))
                grida(i)%d05(grida(i)%kp_cur:grida(i)%kp+1)=0._wp
                grida(i)%d05(0)=grida(i)%d05(1)
            case default
                error stop 'BDM: invalid diffusion_type'
            end select

            call backward_euler(grida(i)%kp,grida(i)%kp_cur,dt_step, &
                grida(i)%r,grida(i)%r05,grida(i)%u,grida(i)%d,grida(i)%d05, &
                grida(i)%dr,grida(i)%dr05,grida(i)%c,grida(i)%cold,flux)
        enddo
    end subroutine finalise_bdm_radial_state
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the equilibrium humidity over a particle        				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the equilibrium humidity over a particle
	!>@param[in] t: temperature
	!>@param[in] mwat: mass of water
	!>@param[in] mbin: mass of aerosol components in each bin
	!>@param[in] nwo: number of moles of water in outer shell in each bin
	!>@param[in] nso: number of moles of solute in outer shell in each bin
	!>@param[in] rhobin: density of each component
	!>@param[in] nubin: van hoff factor in each bin
	!>@param[in] molwbin: molecular weight in each bin
	!>@param[in] sz: length of array
	!>@param[inout] rh_eq: equilibrium humidity
	!>@param[inout] rhoat: density of particle
	!>@param[inout] dw: wet diameter
    subroutine koehler01_diff(t,mwat,mbin,nwo,nso,rhobin,nubin,molwbin,sz,rh_eq,rhoat,dw) 
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat, nwo
      real(wp), dimension(:,:), intent(in) :: mbin,rhobin,nubin,molwbin,nso
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: nw, fac
      real(wp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(wp), intent(in) :: t
      real(wp) :: sigma

      ! calculate the diameter and radius
      nw(:)=mwat(:)/molw_water
      fac(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/max(fac(:),tiny(1._wp))
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))*6._wp / &
             (pi*max(rhoat(:),tiny(1._wp))))**(1._wp/3._wp)
  
      ! calculate surface tension
      sigma=surface_tension(t)

      ! Outer-shell Raoult activity.  Empty numerical shells should not
      ! produce 0/0 or NaNs; assigning zero activity there drives condensation
      ! if the bin is populated and is ignored when npart is negligible.
      fac(:)=nwo(:)+sum(nso(:,:)*nubin(:,:),2)
      rh_eq(:)=0._wp
      where (fac(:) > tiny(1._wp) .and. dw(:) > tiny(1._wp))
          ! Use safe denominators even inside WHERE.  Fortran processors may
          ! evaluate a vector RHS for masked-out elements, so bare /dw and
          ! /fac can still raise IEEE divide-by-zero/invalid flags.
          rh_eq(:)=exp(4._wp*molw_water*sigma/r_gas/t / &
                       max(rhoat(:),tiny(1._wp))/max(dw(:),tiny(1._wp))) * &
                       nwo(:)/max(fac(:),tiny(1._wp))
      end where
       
       
!       rh_eq(:)=exp(4._wp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))/ fac(:)
       

    end subroutine koehler01_diff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! derivatives for a warm parcel model with diffusion                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates rates of change for a warm parcel model
	!>@param[in] neq: length of solution vector
	!>@param[in] tt: time
	!>@param[in] y: solution vector
	!>@param[inout] ydot: derivates calculated
	!>@param[in] rpar: real data coming in
	!>@param[in] ipar: integer data coming in
    subroutine fparcelwarmdiff(neq, tt, y, ydot, rpar, ipar)
        use numerics_type
        use numerics, only : dfsid1,find_pos

        implicit none
        real(wp), intent(inout) :: tt
        real(wp), intent(inout), dimension(neq) :: y, ydot
        integer(i4b), intent(inout) :: neq
        real(wp), intent(inout) :: rpar
        integer(i4b), intent(inout) :: ipar

        ! local variables
        real(wp) :: wv=0._wp, wl=0._wp, wi=0._wp, rm, cpm, &
                  drv=0._wp, dri=0._wp,dri2=0._wp, &
                  rh,t,p,err,sl, w, &
                  te, qve, pe, var, dummy, rhoe, rhop, b
        ! diffusion:
        real(wp) :: tstart, deltaV,flux,radius,radiusold

        integer(i4b) :: i, j,iloc, ipart, ipr, ite, irh, iz,iw
        real(wp), dimension(neq) :: y_eval

        ipart=parcel1%n_bin_modew
        ipr=parcel1%ipr
        ite=parcel1%ite
        irh=parcel1%irh
        iz =parcel1%iz
        iw =parcel1%iw

        rh=y(irh)
        t=y(ite)
        p=y(ipr)
        w=y(iw)

        ! Never modify DVODE's trial state inside the RHS.  Use a tiny positive
        ! evaluation mass for thermodynamics when a finite-difference/Jacobian
        ! trial crosses y=0.
        y_eval=y
        where (y(1:ipart) > 0._wp)
            y_eval(1:ipart)=y(1:ipart)
        elsewhere
            y_eval(1:ipart)=1.e-30_wp
        end where
        
        tstart=parcel1%tout-parcel1%dt ! starting time
    

        ! check there are no negative values
!         where(y(1:ipart).le.0._wp)
!             y(1:ipart)=abs(y(1:ipart))
!         end where


        ! calculate mixing ratios from rh, etc
        sl=svp_liq(t)*rh/(p-svp_liq(t)) ! saturation ratio
        sl=(sl*p/(1._wp+sl))/svp_liq(t)
        wv=eps1*rh*svp_liq(t) / (p-svp_liq(t)) ! vapour mixing ratio
        wl=sum(parcel1%npart*y_eval(1:ipart))        ! liquid mixing ratio
        if (ice_flag == 1) then
            wi=sum(parcel1%npartice*parcel1%yice(1:ipart))
        else
            wi=0._wp
        endif

        ! calculate the moist gas constants and specific heats
        rm=ra+wv*rv
        cpm=cp+wv*cpv+wl*cpw+wi*cpi

        ! now calculate derivatives
        ! adiabatic parcel model
        ydot(iz )=w                         ! vertical wind
        ydot(ipr)=-p/rm/t*grav*ydot(iz)      ! hydrostatic equation


        ! diffusion stuff - part 1
        do i=1,ipart
            
            grida(i)%c=grida(i)%cold
            grida(i)%kp_cur=grida(i)%kp_cur_old
            grida(i)%rad=grida(i)%rad_old
            grida(i)%r=grida(i)%r_old
            grida(i)%r05=grida(i)%r05_old
            grida(i)%dr=grida(i)%dr_old
            grida(i)%dr05=grida(i)%dr05_old
            grida(i)%vol=grida(i)%vol_old
            
            nwo(i)=grida(i)%c(grida(i)%kp_cur,1)
            nso(i,1)=grida(i)%c(grida(i)%kp_cur,2)

            !if(parcel1%npart(i).le. 1.e-9_wp) cycle
            
!             deltaV=max(y(i)- &
!                 sum(grida(i)%c(1:grida(i)%kp_cur,1)*grida(i)%vol(1:grida(i)%kp_cur))*molw_water &
!                 ,-y(i))/rhow
            deltaV=(max(y(i),0._wp)- &
                sum(grida(i)%c(1:grida(i)%kp_cur,1)*grida(i)%vol(1:grida(i)%kp_cur))*molw_water &
                )/rhow
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! shift radii and calculate the velocity of boundaries                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			flux=0._wp
            radius=grida(i)%rad
            radiusold=radius
 			call move_boundary(grida(i)%kp,grida(i)%kp_cur,tt-tstart, &
 			    radiusold,radius,grida(i)%r,grida(i)%r05,grida(i)%dr,grida(i)%dr05, &
 			    grida(i)%vol,grida(i)%u,grida(i)%c,flux, &
 			    grida(i)%rad_min,grida(i)%rad_max, grida(i)%mwsol, grida(i)%rhosol, &
 			    deltaV)
            radiusold=radius
            grida(i)%rad=radius
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			do j=1,1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Set diffusion coefficients - inc. zero at boundary                     !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                select case (diffusion_type)
                    case(0)
                        grida(i)%d05(:)=grida(i)%d_coeff
                        grida(i)%d05(grida(i)%kp_cur:grida(i)%kp+1) = 0._wp
                    case(1)
                        call diffusion_coefficient(grida(i)%kp_cur,n_comps+1, &
                                grida(i)%c(1:grida(i)%kp_cur,1) / &
                                    max(sum(grida(i)%c(1:grida(i)%kp_cur,:),2), &
                                        tiny(1._wp)), &
                                grida(i)%t, d_self, param, &
                                compound, grida(i)%d05(1:grida(i)%kp_cur))
                                
!                         grida(i)%d05(0) = grida(i)%d05(1)
!                         grida(i)%d05(1:grida(i)%kp_cur-1)= &
!                             0.5_wp*(grida(i)%d05(1:grida(i)%kp_cur-1)+ &
!                                     grida(i)%d05(2:grida(i)%kp_cur))
                        grida(i)%d05(grida(i)%kp_cur:grida(i)%kp+1) = 0._wp
                        grida(i)%d05(0) = grida(i)%d05(1)
                    case default
                        print *,'error diffusion type'
                        stop
                end select
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! solve diffusion equation                                               !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call backward_euler(grida(i)%kp,grida(i)%kp_cur,(tt-tstart)/1._wp, &
                    grida(i)%r,grida(i)%r05,grida(i)%u,grida(i)%d,grida(i)%d05,&
                    grida(i)%dr,grida(i)%dr05,grida(i)%c,grida(i)%cold,flux)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo        


 !            nwo(i)=grida(i)%c(grida(i)%kp_cur,1)*grida(i)%vol(grida(i)%kp_cur)/ &
!                 sum(grida(i)%c(1:grida(i)%kp_cur,1)*grida(i)%vol(grida(i)%kp_cur)) *y(i)/molw_water
!             nso(i,1)=grida(i)%c(grida(i)%kp_cur,2)*grida(i)%vol(grida(i)%kp_cur)/ &
!                 sum(grida(i)%c(1:grida(i)%kp_cur,2)*grida(i)%vol(grida(i)%kp_cur)) *parcel1%mbin(i,1)/parcel1%molwbin(i,1)
            nwo(i)=grida(i)%c(grida(i)%kp_cur,1)*grida(i)%vol(grida(i)%kp_cur)
            nso(i,1)=grida(i)%c(grida(i)%kp_cur,2)*grida(i)%vol(grida(i)%kp_cur)
        enddo

        ! calculate equilibrium rhs
        select case (kappa_flag)
            case (0)
                select case (koehler_shell_flag) 
                    case(0) ! standard koehler eq
                        call koehler01(t,y_eval(1:ipart),parcel1%mbin(:,1:n_comps), &
                           parcel1%rhobin(:,1:n_comps), parcel1%nubin(:,1:n_comps), &
                           parcel1%molwbin(:,1:n_comps),ipart, &
                           parcel1%rh_eq,parcel1%rhoat, parcel1%dw) 
                           
                    case(1) ! just use water in outer shell
                        call koehler01_diff(t,y_eval(1:ipart),parcel1%mbin(1:ipart,1:n_comps), &
                           nwo,nso, &
                           parcel1%rhobin(1:ipart,1:n_comps), parcel1%nubin(1:ipart,1:n_comps), &
                           parcel1%molwbin(1:ipart,1:n_comps),ipart, &
                           parcel1%rh_eq,parcel1%rhoat, parcel1%dw) 
                           
                    case default
                        print *,'error koehler_shell_flag'
                        stop
                end select
            case (1)
              call kkoehler01(t,y_eval(1:ipart),parcel1%mbin(:,1:n_comps), &
                   parcel1%rhobin(:,1:n_comps), parcel1%kappabin(:,1:n_comps), &
                   parcel1%molwbin(:,1:n_comps),ipart, &
                   parcel1%rh_eq,parcel1%rhoat, parcel1%dw)
        case default
            print *,'error kappa_flag'
        end select
        
        
        

        ! particle growth rate - radius growth rate
        parcel1%da_dt=dropgrowthrate01(t,p,sl,parcel1%rh_eq, &
            parcel1%rhoat,parcel1%dw,ipart)
        ! do not bother if number concentration too small
        do i=1,ipart
            if(isnan(parcel1%da_dt(i)).or.(parcel1%npart(i).le.1.e-9_wp)) then
              parcel1%da_dt(i)=0._wp
            endif
!             if(isnan(parcel1%da_dt(i)).or.(parcel1%npart(i).le. 1.e-9_wp)) then
!               parcel1%da_dt(i)=0._wp
!             endif
        enddo

!         if((tt > 200._wp) ) then
!             print *,nwo(40)/(nwo(40)+nso(40,1)), y(40),parcel1%da_dt(40),tt,y(irh), parcel1%rh_eq(40)
!         endif

        
        ! mass growth rate
        ydot(1:ipart)=pi*parcel1%rhoat*parcel1%dw**2 * parcel1%da_dt

        ! Keep y=0 as an absorbing evaporation boundary during one DVODE call.
        ! The fixed mbin water mass identifies bins that START this external
        ! step dry, which avoids a discontinuous finite-difference Jacobian at
        ! y=0 while still allowing positive condensation/reactivation.
        where ((parcel1%mbin(1:ipart,n_comps+1) <= 0._wp .or. &
                y(1:ipart) <= 0._wp) .and. ydot(1:ipart) < 0._wp)
            ydot(1:ipart)=0._wp
        end where

        ! change in vapour content
        drv = -sum(ydot(1:ipart)*parcel1%npart)
        if((.not. adiabatic_prof) .and. (.not. vert_ent)) then ! entraining?
            !calculate the environmental p, qv, te, density
            ! parcel p, density
            ! buoyancy...
            ! locate position
            iloc=find_pos(parcel1%z_sound(1:n_levels_s),y(iz))
            iloc=min(n_levels_s-1,iloc)
            iloc=max(1,iloc)
            ! linear interp p
            call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            pe=var
            ! linear interp qv
            call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%q_sound(1,iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            qve=var
            ! linear interp te
            call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            te=var
            ! env density:
            rhoe=pe/(rm*te)
            ! parcel density:
            rhop=p/(rm*t)
            !buoyancy
            if((parcel1%z_sound(n_levels_s) .lt. y(iz)) .or. &
                (parcel1%z_sound(1) .gt. y(iz))) then
                b=0._wp
            else
                b=grav*(rhoe-rhop)/rhoe
            endif
            ! forcing
            drv=drv+w*ent_rate*(qve-wv)

        endif
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in temperature of parcel                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(ite)=rm/p*ydot(ipr)*t/cpm  ! temperature change: expansion
        ydot(ite)=ydot(ite)-lv/cpm*drv ! temp change: condensation
        if((.not. adiabatic_prof) .and. (.not. vert_ent)) then ! entraining?
            ydot(ite)=ydot(ite)+w*ent_rate*(te-y(ite) + lv/cpm*(qve-wv))
            !ydot(iw) = b -w*ent_rate*y(iw)
            ydot(iw) = 0._wp
        else
            ydot(iw) = 0._wp
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in rh of parcel                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(irh)=(p-svp_liq(t))*svp_liq(t)*drv
        ydot(irh)=ydot(irh)+svp_liq(t)*wv*ydot(ipr)
        ydot(irh)=ydot(irh)-wv*p*dfsid1(svp_liq,t,1.e0_wp,1.e-8_wp,err)*ydot(ite)
        ydot(irh)=ydot(irh) / (eps1*svp_liq(t)**2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      
        
    end subroutine fparcelwarmdiff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ice nucleation                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine noncollisional_iceformation_diff(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,moments,medges, &
                         t,p,nbins1,ncomps,nbinw,nmoms,nmodes,yice,rh,dt,sce_flag_in, &
                         mode1_flag_in, ice_nucleation_mech_in)
      ! BDM version of the current BMM func4 callback.  BDM deliberately supports
      ! only full-moving, non-collisional, soluble-aerosol homogeneous freezing.
      ! The freezing probability is evaluated from the radially resolved water
      ! activity in grida; the transfer to the BMM ice population follows the
      ! current full-moving noncollisional_iceformation bookkeeping.
      use numerics_type
      use sce, only : sce_receiving_bin
      implicit none

      real(wp), intent(inout) :: t
      real(wp), intent(in) :: p,dt
      integer(i4b), intent(in) :: nbins1,ncomps,nbinw,nmoms,nmodes,sce_flag_in
      logical, intent(in) :: mode1_flag_in
      logical, dimension(N_INUC_MECH), intent(in) :: ice_nucleation_mech_in

      real(wp), dimension(nbinw), intent(inout) :: npart,npartice
      real(wp), dimension(nbinw), intent(in) :: mwat
      real(wp), dimension(nbinw,ncomps), intent(in) :: mbin2, &
                                              rhobin,nubin,kappabin,molwbin
      real(wp), dimension(2*nbinw,nmoms), intent(inout) :: moments
      real(wp), dimension(nbins1+1,nmodes), intent(in) :: medges
      real(wp), dimension(nbinw,ncomps+1), intent(inout) :: mbin2_ice
      real(wp), dimension(nbinw), intent(inout) :: yice
      real(wp), intent(inout) :: rh

      real(wp), dimension(nbinw) :: dn01,m01
      real(wp), dimension(nmoms) :: momtemp
      real(wp) :: fracinliq, exponent_arg, mleft
      integer(i4b) :: i,j,k,jl,jh,inew

      ! These checks duplicate the startup validation intentionally: they protect
      ! the callback contract if BMM controls are ever changed during a run.
      if (sce_flag_in /= 0_i4b) &
          error stop 'BDM: SCE/aggregation is not supported'
      if (mode1_flag_in) &
          error stop 'BDM: mode-1 secondary ice is not supported'
      if (parcel1%bin_scheme_flag /= BIN_FULL_MOVING) &
          error stop 'BDM: bin_scheme_flag must be full_moving (0)'
      if (ncomps /= 1_i4b) &
          error stop 'BDM: radial diffusion currently supports one soluble aerosol component'
      if (n_inp_classes /= 0_i4b) &
          error stop 'BDM: heterogeneous/INP classes are not supported'
      if (.not.ice_nucleation_mech_in(INUC_KOOP) .or. &
          any(ice_nucleation_mech_in(2:N_INUC_MECH))) &
          error stop 'BDM: only Koop homogeneous ice nucleation is supported'
      if (nbinw /= size(grida)) &
          error stop 'BDM: BMM warm-bin count does not match radial-grid count'

      ! DVODE may finish by interpolation and its last RHS evaluation need
      ! not correspond to the accepted state at tout.  Reconstruct grida from
      ! the accepted warm-water masses before using the radial profile for Koop.
      call finalise_bdm_radial_state(dt)

      dn01=0._wp
      if (t <= ttr) then
          do k=1,nbinw
              if (npart(k) <= qsmall2 .or. mwat(k) <= tiny(1._wp)) cycle
              if (grida(k)%kp_cur < 1) cycle

              ! Radially resolved water activity.  BDM has one soluble material
              ! component in addition to water, hence c(:,1)=water and c(:,2)=solute.
              aw(1:grida(k)%kp_cur)=grida(k)%c(1:grida(k)%kp_cur,1) / &
                  max(grida(k)%c(1:grida(k)%kp_cur,1) + &
                      grida(k)%c(1:grida(k)%kp_cur,2)*nubin(k,1), tiny(1._wp))

              ! Water moles in each radial shell and Koop nucleation rate.
              nw(1:grida(k)%kp_cur)=grida(k)%c(1:grida(k)%kp_cur,1) * &
                                     grida(k)%vol(1:grida(k)%kp_cur)
              jw(1:grida(k)%kp_cur)=koopnucrate(aw(1:grida(k)%kp_cur), &
                                                 t,p,grida(k)%kp_cur)

              ! Integral J dV dt, using water volume to retain the original BDM
              ! interpretation of Koop's volumetric nucleation rate.
              exponent_arg=sum(jw(1:grida(k)%kp_cur)*nw(1:grida(k)%kp_cur)) * &
                           molw_water/rhow*dt
              dn01(k)=npart(k)*(1._wp-exp(-max(exponent_arg,0._wp)))
              dn01(k)=min(max(dn01(k),0._wp),npart(k))
          enddo
      endif

      if (maxval(dn01) <= qsmall2) return

      ! Existing extensive ice water in each moving ice category.
      m01=yice*npartice

      ! Same primary-freezing transfer used by the current BMM full-moving
      ! noncollisional_iceformation routine, with fragmentation disabled.
      do i=1,nmodes
          jl=(i-1)*nbins1+1
          jh=i*nbins1
          do j=1,nbins1
              k=j+(i-1)*nbins1
              if (dn01(k) <= qsmall2) cycle

              npart(k)=max(npart(k)-dn01(k),0._wp)
              fracinliq=npart(k)/max(npart(k)+dn01(k),1.e-30_wp)

              momtemp=0._wp
              momtemp(1:ncomps)=(1._wp-fracinliq)*moments(k,1:ncomps)
              moments(k,1:ncomps)=moments(k,1:ncomps)*fracinliq
              moments(k,ncomps+1)=moments(k,ncomps+1)*fracinliq
              moments(k,ncomps+2)=moments(k,ncomps+2)*fracinliq
              moments(k,ncomps+4)=npart(k)*mwat(k)
              moments(k,ncomps+5)=npart(k)*mwat(k)

              inew=sce_receiving_bin(mwat(k),jl,jh,yice,.true.)

              moments(inew+nbinw,1:ncomps)= &
                  moments(inew+nbinw,1:ncomps)+momtemp(1:ncomps)

              mleft=mwat(k)*dn01(k)
              m01(inew)=m01(inew)+mleft
              npartice(inew)=npartice(inew)+dn01(k)
              if (npartice(inew) > qsmall2) yice(inew)=m01(inew)/npartice(inew)

              ! Ice phi, monomer-number and volume moments.
              moments(inew+nbinw,ncomps+1)= &
                  moments(inew+nbinw,ncomps+1)+dn01(k)
              moments(inew+nbinw,ncomps+2)= &
                  moments(inew+nbinw,ncomps+2)+dn01(k)
              moments(inew+nbinw,ncomps+3)= &
                  moments(inew+nbinw,ncomps+3)+mleft/rhoice
          enddo
      enddo

      where ((m01 > 0._wp) .and. (npartice > 0._wp))
          ! Keep the divisor safe even for masked-out elements; see the note in
          ! koehler01_diff regarding vector evaluation of WHERE expressions.
          yice=m01/max(npartice,tiny(1._wp))
      end where

      mbin2_ice(:,1:ncomps)=moments(nbinw+1:2*nbinw,1:ncomps) / &
          (1.e-50_wp+spread(npartice,2,ncomps))
      mbin2_ice(:,ncomps+1)=yice

      ! Match the current BMM thermodynamic bookkeeping: freezing changes T,
      ! then RH is recomputed from unchanged vapour mass at the updated T.
      call adjust_t_rh(sum(mwat*dn01),t,rh,p)

    end subroutine noncollisional_iceformation_diff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! output to netcdf                                                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>output 1 time-step of model
	!>@param[inout] new_file
    subroutine outputdiff(new_file,outputfile)

    use numerics_type
    use netcdf

    implicit none
    logical, intent(inout) :: new_file
    character(len=*), intent(in) :: outputfile
    
    integer(i4b) :: i
    
    ! output to netcdf file
    if(new_file) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! open / create the netcdf file                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
        call check( nf90_redef(io1%ncid) )
        ! define dimensions (netcdf hands back a handle)

            ! diffusion stuff
            call check( nf90_def_dim(io1%ncid, "kp", grida(1)%kp, io1%y_dimid) )
            call check( nf90_def_dim(io1%ncid, "ncomp", 2, io1%z_dimid) )   
        ! close the file, freeing up any internal netCDF resources
        ! associated with the file, and flush any buffers
        call check( nf90_close(io1%ncid) )


        ! now define some variables, units, etc
        call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
        ! define mode
        call check( nf90_redef(io1%ncid) )

            ! define variable: c
            call check( nf90_def_var(io1%ncid, "c", NF90_DOUBLE, &
                    (/io1%y_dimid, io1%z_dimid, io1%bin_dimid,io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "c", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "moles per m3") )

            ! define variable: r
            call check( nf90_def_var(io1%ncid, "r", NF90_DOUBLE, &
                        (/io1%y_dimid, io1%bin_dimid, io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "r", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "m") )

            ! define variable: vol
            call check( nf90_def_var(io1%ncid, "vol", NF90_DOUBLE, &
                        (/io1%y_dimid, io1%bin_dimid, io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "vol", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "m3") )


                   
        call check( nf90_enddef(io1%ncid) )
        call check( nf90_close(io1%ncid) )

        new_file=.false.
    endif
    io1%icur=io1%icur-1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write data to file                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
    do i=1,n_bins
        ! write variable: c
        call check( nf90_inq_varid(io1%ncid, "c", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, grida(i)%c(1:grida(i)%kp,1:2), &
                    start = (/1,1,i,io1%icur/)))

        ! write variable: r
        call check( nf90_inq_varid(io1%ncid, "r", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, grida(i)%r(1:grida(i)%kp), &
                    start = (/1,i,io1%icur/)))

        ! write variable: vol
        call check( nf90_inq_varid(io1%ncid, "vol", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, grida(i)%vol(1:grida(i)%kp), &
                    start = (/1,i,io1%icur/)))
    enddo


    call check( nf90_close(io1%ncid) )


    io1%icur=io1%icur+1
    end subroutine outputdiff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	end module bdm

