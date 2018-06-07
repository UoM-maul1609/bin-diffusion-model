	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>module for bin microphysics with diffusion model
	module bdm
    use nrtype
    use nr, only : locate, polint
    use diffusion, only : grid, backward_euler, move_boundary
    use bmm
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the bin diffusion model

    implicit none


        character (len=200) :: outputfile='output', bin_model_file='', &
                            diffusion_model_file=''
                            
        type(grid), allocatable, dimension(:) :: grida
        
        real(sp), allocatable, dimension(:) :: nwo
        real(sp), allocatable, dimension(:,:) :: nso
        integer(i4b) :: koehler_shell_flag


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
		implicit none
        character (len=200), intent(in) :: nmlfile
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ outputfile, bin_model_file, diffusion_model_file, &
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
         
	end subroutine read_in_bdm_namelist







	
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

        
        do i=1,n_bins
            call allocate_and_set_diff(nmd%kp,nmd%dt,nmd%runtime, &
                parcel1%dw(i)/2._sp, nmd%rad_min, nmd%rad_max, nmd%t,nmd%p, parcel1%rh, &
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
    implicit none
    integer(i4b) :: i, j,nt
    logical :: new_file=.true.
    
    nt=ceiling(runtime / real(dt,kind=sp))
    do i=1,nt
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
        enddo
        
        ! one time-step of model
        call bin_microphysics(fparcelwarmdiff,fparcelcold)
        
        ! check there are no negative values
        where(parcel1%y(1:parcel1%n_bin_mode).le.0.e1_sp)
            parcel1%y(1:parcel1%n_bin_mode)=1.e-22_sp
        end where

             
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
      real(sp), dimension(sz) :: nw
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

      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      rh_eq(:)=exp(4._sp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))* &
           (nwo(:))/(nwo(:)+sum(nso(:,:)*nubin(:,:),2) ) 

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
        where(y(1:ipart).le.0.e1_sp)
            y(1:ipart)=1.e-22_sp
        end where


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
            
            deltaV=(y(i)-parcel1%yold(i))/rhow
        
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! shift radii and calculate the velocity of boundaries                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 			call move_boundary(grida(i)%kp,grida(i)%kp_cur,tt-tstart, &
 			    radiusold,radius,grida(i)%r,grida(i)%r05,grida(i)%dr,grida(i)%dr05, &
 			    grida(i)%vol,grida(i)%u,grida(i)%c,flux, &
 			    grida(i)%rad_min,grida(i)%rad_max, grida(i)%mwsol, grida(i)%rhosol, &
 			    deltaV)
            radiusold=radius
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Set diffusion coefficient to zero at boundary                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			grida(i)%d05(:)=grida(i)%d_coeff
			grida(i)%d05(grida(i)%kp_cur:grida(i)%kp) = 0._sp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! solve diffusion equation                                                   !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call backward_euler(grida(i)%kp,grida(i)%kp_cur,tt-tstart, &
			    grida(i)%r,grida(i)%r05,grida(i)%u,grida(i)%d,grida(i)%d05,&
			    grida(i)%dr,grida(i)%dr05,grida(i)%c,grida(i)%cold,flux)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            nwo(i)=grida(i)%c(grida(i)%kp_cur,1)
            nso(i,1)=grida(i)%c(grida(i)%kp_cur,2)
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
                        call koehler01_diff(t,y(1:ipart),parcel1%mbin(:,1:n_comps), &
                           nwo,nso, &
                           parcel1%rhobin(:,1:n_comps), parcel1%nubin(:,1:n_comps), &
                           parcel1%molwbin(:,1:n_comps),ipart, &
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
            if(isnan(parcel1%da_dt(i)).or.parcel1%npart(i).le. 1.e-9_sp) then
              parcel1%da_dt(i)=0._sp
            endif
        enddo


        
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
    enddo


    call check( nf90_close(io1%ncid) )


    io1%icur=io1%icur+1
    end subroutine outputdiff
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	end module bdm

