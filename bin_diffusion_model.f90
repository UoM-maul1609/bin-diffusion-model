	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>module for bin microphysics with diffusion model
	module bdm
    use nrtype
    use nr, only : locate, polint
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
        
        real(sp), allocatable, dimension(:) :: nwo, jw, nw,ns, aw
        real(sp), allocatable, dimension(:,:) :: nso
        integer(i4b) :: koehler_shell_flag
        integer(i4b) :: diffusion_type
        
        ! for calculating diffusion coefficients
        real(sp), allocatable, dimension (:) :: d_self
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
        real(sp), allocatable, dimension (:), intent(out) :: d_self
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
        use nrtype
        use nr, only : locate, polint, rkqs, odeint, zbrent, brent
        use bmm, only : initialise_bmm_arrays, nu_core1
        use diffusion, only : allocate_and_set_diff, nmd, gridd

        implicit none
        integer(i4b) :: AllocateStatus,i
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! bin-microphysics arrays                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call initialise_bmm_arrays()
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! work out the diameter, including the water...                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call wetdiam(parcel1%mbin(:,n_comps+1),parcel1%mbin,parcel1%rhobin, &
                    parcel1%n_bin_mode,parcel1%dw) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays for diffusion model         		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		allocate( grida(1:n_bins), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nwo(1:n_bins), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nso(1:n_bins,1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

		allocate( jw(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nw(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( ns(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( aw(1:nmd%kp), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

        
        do i=1,n_bins
            call allocate_and_set_diff(nmd%kp,nmd%dt,nmd%runtime, &
!                parcel1%dw(i)/2._sp, nmd%rad_min, nmd%rad_max, nmd%t,nmd%p, parcel1%rh, &
                parcel1%dw(i)/2._sp, nmd%rad_min, nmd%rad_max, nmd%t,nmd%p, &
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
    use nrtype
    use bmm, only : parcel, parcel1
    implicit none
    integer(i4b) :: i, j,nt
    logical :: new_file=.true.
    real(sp) :: flux,deltaV, radius, radiusold, t_tstep
    
    
    ! reduce tolerances for this model:
    parcel1%rtol=1.e-6_sp
    parcel1%atol(1:parcel1%n_bin_mode)=1.e-30_sp
    
    
    
    nt=ceiling(runtime / real(dt,kind=sp))
    do i=1,nt
        t_tstep=parcel1%y(parcel1%ite) ! use a constant temperature for the diffusion 
                                       ! over the time-step
    
        ! output to file
        call output(io1%new_file,outputfile)
        call outputdiff(new_file,outputfile)
        
        
        
        ! store old aerosol state
        do j=1,n_bins
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
        call bin_microphysics(fparcelwarmdiff,fparcelcold,icenucleation_diff)
        
        
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
! 			grida(j)%d05(grida(j)%kp_cur:grida(j)%kp) = 0._sp
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
!         where(parcel1%y(1:parcel1%n_bin_mode).le.0.e1_sp)
!             parcel1%y(1:parcel1%n_bin_mode)=1.e-22_sp
!         end where

        do j=1,parcel1%n_bin_mode
            if(parcel1%y(j).le.0.e1_sp) then
                parcel1%y(j)=sum(grida(j)%c(1:grida(j)%kp_cur,1)* &
                                            grida(j)%vol(1:grida(j)%kp_cur))*molw_water
            endif
        enddo
        
        ! break-out if flag has been set 
        if(parcel1%break_flag) exit
    enddo
    ! output to file
    call output(io1%new_file,outputfile)
    call outputdiff(new_file,outputfile)
    
    
    
    
    end subroutine bdm_driver
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
      use nrtype
      implicit none
      real(sp), dimension(:), intent(in) :: mwat, nwo
      real(sp), dimension(:,:), intent(in) :: mbin,rhobin,nubin,molwbin,nso
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz) :: nw, fac
      real(sp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(sp), intent(in) :: t
      real(sp) :: sigma

      ! calculate the diameter and radius
      nw(:)=mwat(:)/molw_water
      rhoat(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat(:);
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))*6._sp/(pi*rhoat(:)))**(1._sp/3._sp)
  
      ! calculate surface tension
      sigma=surface_tension(t)
        
      !fac(:)=1._sp+max(sum(nso(:,:)*nubin(:,:),2)/nwo(:),0.1_sp)
      
      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      rh_eq(:)=exp(4._sp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))* &
           (nwo(:))/(nwo(:)+sum(nso(:,:)*nubin(:,:),2) ) 
       
       
!       rh_eq(:)=exp(4._sp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))/ fac(:)
       

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
        use nrtype
        use nr, only : dfridr,locate

        implicit none
        real(sp), intent(inout) :: tt
        real(sp), intent(inout), dimension(neq) :: y, ydot
        integer(i4b), intent(inout) :: neq
        real(sp), intent(inout) :: rpar
        integer(i4b), intent(inout) :: ipar

        ! local variables
        real(sp) :: wv=0._sp, wl=0._sp, wi=0._sp, rm, cpm, &
                  drv=0._sp, dri=0._sp,dri2=0._sp, &
                  rh,t,p,err,sl, w, &
                  te, qve, pe, var, dummy, rhoe, rhop, b
        ! diffusion:
        real(sp) :: tstart, deltaV,flux,radius,radiusold

        integer(i4b) :: i, j,iloc, ipart, ipr, ite, irh, iz,iw

        ipart=parcel1%n_bin_mode
        ipr=parcel1%ipr
        ite=parcel1%ite
        irh=parcel1%irh
        iz =parcel1%iz
        iw =parcel1%iw

        rh=y(irh)
        t=y(ite)
        p=y(ipr)
        w=y(iw)
        
        tstart=parcel1%tout-parcel1%dt ! starting time
    

        ! check there are no negative values
!         where(y(1:ipart).le.0._sp)
!             y(1:ipart)=abs(y(1:ipart))
!         end where


        ! calculate mixing ratios from rh, etc
        sl=svp_liq(t)*rh/(p-svp_liq(t)) ! saturation ratio
        sl=(sl*p/(1._sp+sl))/svp_liq(t)
        wv=eps1*rh*svp_liq(t) / (p-svp_liq(t)) ! vapour mixing ratio
        wl=sum(parcel1%npart*y(1:ipart))             ! liquid mixing ratio

        ! calculate the moist gas constants and specific heats
        rm=ra+wv*rv
        cpm=cp+wv*cpv+wl*cpw+wi*cpi

        ! now calculate derivatives
        ! adiabatic parcel model
        ydot(iz )=w                         ! vertical wind
        ydot(ipr)=-p/rm/t*grav*ydot(iz)      ! hydrostatic equation


        ! diffusion stuff - part 1
        do i=1,n_bins
            
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

            !if(parcel1%npart(i).le. 1.e-9_sp) cycle
            
!             deltaV=max(y(i)- &
!                 sum(grida(i)%c(1:grida(i)%kp_cur,1)*grida(i)%vol(1:grida(i)%kp_cur))*molw_water &
!                 ,-y(i))/rhow
            deltaV=(y(i)- &
                sum(grida(i)%c(1:grida(i)%kp_cur,1)*grida(i)%vol(1:grida(i)%kp_cur))*molw_water &
                )/rhow
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! shift radii and calculate the velocity of boundaries                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			flux=0._sp
 			call move_boundary(grida(i)%kp,grida(i)%kp_cur,tt-tstart, &
 			    radiusold,radius,grida(i)%r,grida(i)%r05,grida(i)%dr,grida(i)%dr05, &
 			    grida(i)%vol,grida(i)%u,grida(i)%c,flux, &
 			    grida(i)%rad_min,grida(i)%rad_max, grida(i)%mwsol, grida(i)%rhosol, &
 			    deltaV)
            radiusold=radius
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			do j=1,1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Set diffusion coefficients - inc. zero at boundary                     !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                select case (diffusion_type)
                    case(0)
                        grida(i)%d05(:)=grida(i)%d_coeff
                        grida(i)%d05(grida(i)%kp_cur:grida(i)%kp+1) = 0._sp
                    case(1)
                        call diffusion_coefficient(grida(i)%kp_cur,n_comps+1, &
                                grida(i)%c(1:grida(i)%kp_cur,1) / &
                                    sum(grida(i)%c(1:grida(i)%kp_cur,:),2), &
                                grida(i)%t, d_self, param, &
                                compound, grida(i)%d05(1:grida(i)%kp_cur))
                                
!                         grida(i)%d05(0) = grida(i)%d05(1)
!                         grida(i)%d05(1:grida(i)%kp_cur-1)= &
!                             0.5_sp*(grida(i)%d05(1:grida(i)%kp_cur-1)+ &
!                                     grida(i)%d05(2:grida(i)%kp_cur))
                        grida(i)%d05(grida(i)%kp_cur:grida(i)%kp+1) = 0._sp
                        grida(i)%d05(0) = grida(i)%d05(1)
                    case default
                        print *,'error diffusion type'
                        stop
                end select
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! solve diffusion equation                                               !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call backward_euler(grida(i)%kp,grida(i)%kp_cur,(tt-tstart)/1._sp, &
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
                        call koehler01(t,y(1:ipart),parcel1%mbin(:,1:n_comps), &
                           parcel1%rhobin(:,1:n_comps), parcel1%nubin(:,1:n_comps), &
                           parcel1%molwbin(:,1:n_comps),ipart, &
                           parcel1%rh_eq,parcel1%rhoat, parcel1%dw) 
                           
                    case(1) ! just use water in outer shell
                        call koehler01_diff(t,abs(y(1:ipart)),parcel1%mbin(1:ipart,1:n_comps), &
                           nwo,nso, &
                           parcel1%rhobin(1:ipart,1:n_comps), parcel1%nubin(1:ipart,1:n_comps), &
                           parcel1%molwbin(1:ipart,1:n_comps),ipart, &
                           parcel1%rh_eq,parcel1%rhoat, parcel1%dw) 
                           
                    case default
                        print *,'error koehler_shell_flag'
                        stop
                end select
            case (1)
              call kkoehler01(t,y(1:ipart),parcel1%mbin(:,1:n_comps), &
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
            if(isnan(parcel1%da_dt(i))) then
              parcel1%da_dt(i)=0._sp
            endif
!             if(isnan(parcel1%da_dt(i)).or.(parcel1%npart(i).le. 1.e-9_sp)) then
!               parcel1%da_dt(i)=0._sp
!             endif
        enddo

!         if((tt > 200._sp) ) then
!             print *,nwo(40)/(nwo(40)+nso(40,1)), y(40),parcel1%da_dt(40),tt,y(irh), parcel1%rh_eq(40)
!         endif

        
        ! mass growth rate
        ydot(1:ipart)=pi*parcel1%rhoat*parcel1%dw**2 * parcel1%da_dt
        ! change in vapour content
        drv = -sum(ydot(1:ipart)*parcel1%npart)
        if((.not. adiabatic_prof) .and. (.not. vert_ent)) then ! entraining?
            !calculate the environmental p, qv, te, density
            ! parcel p, density
            ! buoyancy...
            ! locate position
            iloc=locate(parcel1%z_sound(1:n_levels_s),y(iz))
            iloc=min(n_levels_s-1,iloc)
            iloc=max(1,iloc)
            ! linear interp p
            call polint(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            pe=var
            ! linear interp qv
            call polint(parcel1%z_sound(iloc:iloc+1), parcel1%q_sound(1,iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            qve=var
            ! linear interp te
            call polint(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            te=var
            ! env density:
            rhoe=pe/(rm*te)
            ! parcel density:
            rhop=p/(rm*t)
            !buoyancy
            if((parcel1%z_sound(n_levels_s) .lt. y(iz)) .or. &
                (parcel1%z_sound(1) .gt. y(iz))) then
                b=0._sp
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
            ydot(iw) = 0._sp
        else
            ydot(iw) = 0._sp
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in rh of parcel                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(irh)=(p-svp_liq(t))*svp_liq(t)*drv
        ydot(irh)=ydot(irh)+svp_liq(t)*wv*ydot(ipr)
        ydot(irh)=ydot(irh)-wv*p*dfridr(svp_liq,t,1.e0_sp,err)*ydot(ite)
        ydot(irh)=ydot(irh) / (eps1*svp_liq(t)**2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      
        
    end subroutine fparcelwarmdiff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ice nucleation                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine icenucleation_diff(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,t,p,sz,sz2,yice,rh,dt) 
      use nrtype
      implicit none
      real(sp), intent(inout) :: t
      real(sp), intent(in) :: p,rh,dt
      real(sp), dimension(sz2), intent(inout) :: npart,npartice
      real(sp), dimension(:), intent(in) :: mwat
      real(sp), dimension(:,:), intent(in) :: mbin2, &
                                              rhobin,nubin,kappabin,molwbin
      integer(i4b), intent(in) :: sz,sz2
      real(sp), dimension(sz2) :: dn01,m01,dw,dd,kappa,rhoat
      real(sp), dimension(sz2,sz) :: dmaer01
      real(sp), dimension(sz2,sz), intent(inout) :: mbin2_ice
      
      real(sp), intent(inout), dimension(sz2) :: yice
      integer(i4b) :: i
      
      
      ! loop over each particle:
      do i=1,sz2
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! first calculate the ice formation over dt using koop et al. 2000       !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! number of moles of water
          nw(1:grida(i)%kp_cur)=grida(i)%c(1:grida(i)%kp_cur,1)* &
                                grida(i)%vol(1:grida(i)%kp_cur)
          ! number of moles of solute
          ns(1:grida(i)%kp_cur)=grida(i)%c(1:grida(i)%kp_cur,2)* &
                                grida(i)%vol(1:grida(i)%kp_cur)
          ! activity of water
          select case(kappa_flag)
            case(0)
              ! beware of underflow here:
              aw(1:grida(i)%kp_cur)=(grida(i)%c(1:grida(i)%kp_cur,1))/ &
                (grida(i)%c(1:grida(i)%kp_cur,1)+ &
                (grida(i)%c(1:grida(i)%kp_cur,2))*nubin(i,1) )
              !aw(:)=(nw(:))/(nw(:)+(ns(:))*nubin(i,1) )
              !aw(:)=(mwat(i)/molw_water)/(mwat(i)/molw_water+mbin2(i,1)/molwbin(i,1)*nubin(i,1) )
            case(1)
              rhoat(i)=mwat(i)/rhow+sum(mbin2(i,:)/rhobin(i,:))
              rhoat(i)=(mwat(i)+sum(mbin2(i,:)))/rhoat(i);
  
              dw(i)=((mwat(i)+sum(mbin2(i,:)))*6._sp/(pi*rhoat(i)))**(1._sp/3._sp)
  
              dd(i)=((sum(mbin2(i,:)/rhobin(i,:)))* &
                     6._sp/(pi))**(1._sp/3._sp) ! dry diameter
                                  ! needed for eqn 6, petters and kreidenweis (2007)
              kappa(i)=sum((mbin2(i,:)+1.e-60_sp)/rhobin(i,:)*kappabin(i,:)) &
                     / sum((mbin2(i,:)+1.e-60_sp)/rhobin(i,:))
                     ! equation 7, petters and kreidenweis (2007)
              aw(i)=(dw(i)**3-dd(i)**3)/(dw(i)**3-dd(i)**3*(1._sp-kappa(i))) ! from eq 6,p+k(acp,2007)
            case default
              print *,'error kappa_flag'
              stop
          end select
          ! koop et al. (2000) nucleation rate - due to homogeneous nucleation.
          jw(1:grida(i)%kp_cur)=koopnucrate(aw(1:grida(i)%kp_cur),t,p,grida(i)%kp_cur)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
          ! the number of ice crystals nucleated:
          dn01(i)=max( npart(i)* &
                (1._sp-exp(-sum(jw(1:grida(i)%kp_cur)* &
                nw(1:grida(i)%kp_cur))*molw_water/rhow*dt)) ,0._sp)
      enddo
      
      



      if(t.gt.ttr) then
          dn01=0._sp
      endif
      !!!!
      ! total aerosol mass in each bin added together:
      dmaer01(:,:)=(mbin2_ice(:,:)*(spread(npartice(:),2,sz)+1.e-50_sp)+ &
                      mbin2(:,:)*spread(dn01(:),2,sz) ) 
      ! total water mass that will be in the ice bins:
      m01=(yice*npartice+mwat(:)*dn01(:)) 

      ! number conc. of liquid bins:
      npart(:)=npart(:)-dn01(:)
      ! number conc. of ice bins:
      npartice(:)=npartice(:)+dn01(:)
      ! new ice mass in bin:
      m01=m01/(npartice) 
      
      
      where(m01.gt.0._sp.and.npartice.gt.0._sp)
        yice=m01
      elsewhere
        yice=yice
      end where
      
      ! aerosol mass in ice bins
      mbin2_ice(:,:)=dmaer01(:,:)/(1.e-50_sp+spread(npartice,2,sz))

      ! latent heat of fusion:
      t=t+lf/cp*sum(mwat(:)*dn01(:))
    end subroutine icenucleation_diff
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

    use nrtype
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

