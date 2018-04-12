  !> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

  double precision :: M_sun = 1.989d33
  double precision :: L_sun = 3.827d33
  double precision :: R_sun = 6.96d10

  double precision :: M_star
  double precision :: L_star
  double precision :: T_star
  double precision :: R_star

  double precision :: Flux0, c_sound0, T_star0, kappa0
  double precision :: c_light0, g0, geff0, heff0, Gamma
  double precision :: L_star0, R_star0, M_star0
  double precision :: tau_bound,  P_bound, rho_bound

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_constants

    call set_coordinate_system("Cartesian_2D")

    M_star = 150*M_sun
    L_star = (M_star/M_sun)**3.d0*L_sun !*unit_time/unit_pressure*unit_length**3.d0
    R_star = 30*R_sun
    T_star = (L_star/(4d0*dpi*R_star**2*5.67051d-5))**0.25d0
    tau_bound = 50.d0

    call initglobaldata_usr()

    !Fix dimensionless stuff here
    unit_length        = R_star
    unit_numberdensity = 4.6d-8/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature   = T_star

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Special Boundary conditions
    usr_special_bc => special_bound

    ! Keep the internal energy constant with internal bound
    usr_internal_bc => constant_e

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call hd_activate()

    print*, 'unit_time', unit_time
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_numberdensity', unit_numberdensity
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_pressure', unit_pressure
    print*, '================================================================'

  end subroutine usr_init

!========================== ================================================================

subroutine initglobaldata_usr
  use mod_global_parameters

  L_star0 = L_star*unit_time/(unit_pressure*unit_length**3.d0)
  R_star0 = R_star/unit_length
  M_star0 = M_star/(unit_density*unit_length**3.d0)

  Flux0 = L_star0/(4*dpi*R_star0**2)
  T_star0 = T_star/unit_temperature
  c_sound0 = dsqrt((1.38d-16*T_star/(0.6*mp_cgs)))/unit_velocity

  kappa0 = 0.34*(unit_density*unit_length**3.d0)/unit_length**2.d0
  c_light0 = const_c/unit_velocity
  g0 = 6.67e-8*M_star/R_star**2*(unit_density*unit_length/unit_pressure)
  geff0 = g0*(one - (kappa0*Flux0)/(c_light0*g0))
  heff0 = c_sound0**2/geff0
  Gamma = (kappa0*Flux0)/(c_light0*g0)

  P_bound = geff0*tau_bound/kappa0
  rho_bound = P_bound/c_sound0**two

end subroutine initglobaldata_usr

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_variables
    use mod_hd_phys, only: hd_get_pthermal

    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
       ixmin2,ixmax1,ixmax2
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)
    double precision :: density(ixGmin1:ixGmax1,ixGmin2:ixGmax2),&
        pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2), pert(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2), amplitude
    integer :: i


    amplitude = 5.d-2
    pressure(:,ixGmin2) = p_bound
    density(:,ixGmin2) = rho_bound

    do i=ixGmin2,ixGmax2
      pressure(:,i) = p_bound*exp(-x(:,i,2)/heff0)
      density(:,i) = rho_bound*exp(-x(:,i,2)/heff0)
    enddo

    ! Set initial values for w
    call RANDOM_NUMBER(pert)
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = density(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2) !*(one + amplitude*pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2))
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero

    call RANDOM_NUMBER(pert)
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = pressure(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2)/(hd_gamma - one) !*(one + amplitude*pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2))
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = &
       3.d0*Gamma/(one-Gamma)*pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2)

    print*, "R_star", R_star0, L_star0
    print*, "R_star", R_star, L_star
    print*, "Flux", Flux0

    print*, "g0", g0 *unit_length/unit_time**2
    print*, "geff0", geff0 *unit_length/unit_time**2
    print*, "c_sound0", c_sound0 *unit_length/unit_time
    print*, "Gamma", Gamma
    print*, "heff0", heff0 *unit_length/ R_star, "Stellar Radii"
    print*, "heff0", heff0 *unit_length, "cm"
    print*, "Tstar0", T_star0
    print*, "Tstar", T_star

    ! print*, "density", w(5,3:10,rho_) *unit_density
    ! print*, "energy", w(5,3:10,e_) *unit_pressure
    ! print*, "rad_energy", w(5,3:10,r_e) *unit_pressure

    print*, rho_bound*unit_density, p_bound*unit_pressure
    print*, "factor", 3.d0*Gamma/(one-Gamma)

    do i=ixGmin2,ixGmax2
      print*, x(5,i,2),w(5,i,rho_)*unit_density
    enddo

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

  subroutine special_bound(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,ixBmin2,&
     ixBmax1,ixBmax2,iB,w,x)

    use mod_global_parameters
    use mod_variables
    use mod_physics, only: phys_get_pthermal

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,ixBmin2,&
       ixBmax1,ixBmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision :: velocity(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndir),&
        pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2)

    select case (iB)

    case(3)

      w(:,ixBmax2, rho_) = rho_bound
      w(:,ixBmax2, mom(1)) = zero

      velocity(:,ixBmax2,2) = 2*(w(:,ixBmax2+1,mom(2))/w(:,ixBmax2+1,&
         rho_) - w(:,ixBmax2+2,mom(2))/w(:,ixBmax2+2,rho_))
      velocity(:,ixBmin2,2) = 2*(w(:,ixBmax2+1,mom(2))/w(:,ixBmax2+1,&
         rho_) - 2*w(:,ixBmax2+2,mom(2))/w(:,ixBmax2+2,rho_))

      ! w(:,ixBmax2, mom(2)) =  w(:,ixBmax2+1, mom(2))

      w(:,ixBmax2, e_) = p_bound/(hd_gamma-one)
      w(:,ixBmax2, r_e) = 3.d0*Gamma/(one-Gamma)*p_bound

      w(:,ixBmin2,:) =  w(:,ixBmax2,:)

      w(:,ixBmax2, mom(2)) = velocity(:,ixBmax2,2)*rho_bound
      w(:,ixBmin2, mom(2)) = velocity(:,ixBmin2,2)*rho_bound

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound


  !==========================================================================================

    !> internal boundary, user defined
    !
    !> This subroutine can be used to artificially overwrite ALL conservative
    !> variables in a user-selected region of the mesh, and thereby act as
    !> an internal boundary region. It is called just before external (ghost cell)
    !> boundary regions will be set by the BC selection. Here, you could e.g.
    !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
    !> which can be used to identify the internal boundary region location.
    !> Its effect should always be local as it acts on the mesh.

    subroutine constant_e(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,w,x)

      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nw)
      double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndim)
      double precision :: pressure(ixImin1:ixImax1,ixImin2:ixImax2)

      pressure(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,&
         ixImin2:ixImax2,rho_)*c_sound0**2
      w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = pressure(ixImin1:ixImax1,&
         ixImin2:ixImax2)/(hd_gamma - one)

    end subroutine constant_e

  !==========================================================================================

  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ndim)

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       2) = -6.67e-8*M_star/R_star**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field


  !==========================================================================================

  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: rad_flux(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim), rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2),&
        fld_lambda(ixImin1:ixImax1,ixImin2:ixImax2), fld_R(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, rad_flux)
    call fld_get_radpress(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, rad_pressure)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, fld_lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)/(unit_pressure*unit_velocity)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)/(unit_pressure*unit_velocity)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=rad_pressure(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4)=fld_lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5)=fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'F1 F2 RP lam fld_R'

  end subroutine specialvarnames_output

  !==========================================================================================

  end module mod_usr

  !==========================================================================================
