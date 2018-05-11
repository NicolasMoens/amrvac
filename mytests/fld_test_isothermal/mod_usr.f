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
  tau_bound = 100.d0

  call initglobaldata_usr()

  !Fix dimensionless stuff here
  unit_length        = R_star
  unit_numberdensity = 8.955d-8/((1.d0+4.d0*He_abundance)*mp_cgs)
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

!==========================================================================================

subroutine initglobaldata_usr
use mod_global_parameters

L_star0 = L_star*unit_time/(unit_pressure*unit_length**3.d0)
R_star0 = R_star/unit_length
M_star0 = M_star/(unit_density*unit_length**3.d0)

Flux0 = L_star0/(4*dpi*R_star0**2)
T_star0 = T_star/unit_temperature
c_sound0 = dsqrt((1.38d-16*T_star/(0.6*mp_cgs)))/unit_velocity

kappa0 = fld_kappa0 !0.34*(unit_density*unit_length**3.d0)/unit_length**2.d0
c_light0 = const_c/unit_velocity
g0 = 6.67e-8*M_star/R_star**2*(unit_density*unit_length/unit_pressure)
geff0 = g0*(one - (kappa0*Flux0)/(c_light0*g0))
heff0 = c_sound0**2/geff0
Gamma = (kappa0*Flux0)/(c_light0*g0)

P_bound = geff0*tau_bound/kappa0
rho_bound = P_bound/c_sound0**two

end subroutine initglobaldata_usr

!==========================================================================================

<<<<<<< HEAD
!> A routine for specifying initial conditions
subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
   ixmax1,ixmax2, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables
  use mod_hd_phys, only: hd_get_pthermal

  integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
     ixmin2,ixmax1,ixmax2
  double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)
  double precision :: density(ixGmin1:ixGmax1,ixGmin2:ixGmax2),&
      pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2), pert(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2), amplitude
  double precision :: rad_Flux(ixmin1:ixmax1,ixmin2:ixmax2, 1:ndim)
  double precision :: opacity(ixmin1:ixmax1,ixmin2:ixmax2),&
      Gamma_dep(ixmin1:ixmax1,ixmin2:ixmax2)
  integer :: i

  amplitude = 5.d-2 !1.d-5 !3.d-2

  pressure(:,ixGmin2) = p_bound
  density(:,ixGmin2) = rho_bound

  do i=ixGmin2,ixGmax2
    pressure(:,i) = p_bound*dexp(-x(:,i,2)/heff0)
    density(:,i) = pressure(:,i)/c_sound0**2 !rho_bound*dexp(-x(:,i,2)/heff0) !
  enddo

  ! Set initial values for w
  !call RANDOM_NUMBER(pert)
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = density(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2) !*(one + amplitude*pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2))
  !> Updqte pressure due to density fluctuations
  !pressure(ixG^S) = c_sound0**two*w(ixG^S, rho_) !> FIX THIS
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)/(hd_gamma - one)
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = &
     3.d0*Gamma/(one-Gamma)*pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2) !> CHANGEd

  !---------------------------------------------------------------------------
  ! Call fld_kappa to calculate correct, Opacity dependent Gamma for initial conditions
  call fld_get_radflux(w,x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,rad_Flux) !> CHANGEd
  call fld_get_opacity(w,x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,opacity)

  Gamma_dep(ixmin1:ixmax1,ixmin2:ixmax2) = opacity(ixmin1:ixmax1,&
     ixmin2:ixmax2)*rad_Flux(ixmin1:ixmax1,ixmin2:ixmax2,2)/(c_light0*g0) !> CHANGEd

  w(ixmin1:ixmax1,ixmin2:ixmax2,r_e) = 3.d0*Gamma_dep(ixmin1:ixmax1,&
     ixmin2:ixmax2)/(one-Gamma_dep(ixmin1:ixmax1,&
     ixmin2:ixmax2))*pressure(ixmin1:ixmax1,ixmin2:ixmax2)
  !---------------------------------------------------------------------------


  !> perturb rho
  call RANDOM_NUMBER(pert)

  ! do i=ixGmin2,ixGmax2
  !   if (i .gt. 0.5d0*ixmax2) then
  !     pert(:,i) = zero
  !   endif
  !   if (i .lt. 0.25d0*ixmax2) then
  !     pert(:,i) = zero
  !   endif
  !   if (i .gt. 0.75d0*ixmax1) then
  !     pert(i,:) = zero
  !   endif
  !   if (i .lt. 0.25d0*ixmax1) then
  !     pert(i,:) = zero
  !   endif
  ! enddo

  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = density(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)*(one + amplitude*pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2))

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

end subroutine initial_conditions
=======
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
    double precision :: Gamma_dep(ixmin1:ixmax1,ixmin2:ixmax2),&
        fld_kappa(ixmin1:ixmax1,ixmin2:ixmax2)
    integer :: i

    double precision                   :: fld_lambda(ixmin1:ixmax1,&
       ixmin2:ixmax2), fld_R(ixmin1:ixmax1,ixmin2:ixmax2)

    amplitude = zero !1.d0-5 !3.d-2

    pressure(:,ixGmin2) = p_bound
    density(:,ixGmin2) = rho_bound

    do i=ixGmin2,ixGmax2
      pressure(:,i) = p_bound*dexp(-x(:,i,2)/heff0)
      density(:,i) = pressure(:,i)/c_sound0**2 !rho_bound*dexp(-x(:,i,2)/heff0) !
    enddo

    ! Set initial values for w
    call RANDOM_NUMBER(pert)
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = density(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2)*(one + amplitude*pert(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2))
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = pressure(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2)/(hd_gamma - one)

    call fld_get_opacity(w,x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
       ixmax1,ixmax2,fld_kappa)

    Gamma_dep = (fld_kappa*Flux0)/(c_light0*g0)

    w(ixmin1:ixmax1,ixmin2:ixmax2,r_e) = 3.d0*Gamma_dep(ixmin1:ixmax1,&
       ixmin2:ixmax2)/(one-Gamma_dep(ixmin1:ixmax1,&
       ixmin2:ixmax2))*pressure(ixmin1:ixmax1,ixmin2:ixmax2)

    !w(ixG^S,r_e) = w(ixG^S,r_e)*(one + dsin(x(ixG^S,1)))

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

  end subroutine initial_conditions
>>>>>>> 924324e1c3d570855f9b58bd582c5da8c596bec5

!==========================================================================================

! Extra routines can be placed here
! ...

!==========================================================================================

<<<<<<< HEAD
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
  double precision :: fld_R(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: fld_kappa(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: Gamma_dep(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: fld_lambda(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: a(1:nw), b(1:nw), c(1:nw)
  integer :: i,j

  select case (iB)

  case(3)

    do i = ixBmin2,ixBmax2
      w(:,i, rho_) = p_bound*dexp(-x(:,i,2)/heff0)/c_sound0**2
      w(:,i, mom(1)) = zero
      velocity(:,i,2) = two*w(:,i+1,mom(2))/w(:,i+1,rho_) - w(:,i+2,&
         mom(2))/w(:,i+2,rho_)
      w(:,i, mom(2)) = velocity(:,i,2)*w(:,i, rho_)
      w(:,i, e_) = p_bound*dexp(-x(:,i,2)/heff0)/(hd_gamma-one)
      !w(:,i, r_e) = 3.d0*Gamma/(one-Gamma)*p_bound*dexp(-x(:,i,2)/heff0)
    enddo

    !> Fixing the R_E boundary for correct gamma due to opacity fluctuation
    !-------------------------------------------------------------------------
    call fld_get_fluxlimiter(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixGmin1+2,&
       ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_lambda, fld_R)

    Gamma_dep(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2) = one/(one + &
       (3.d0*p_bound*dexp(-x(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2,&
       2)/heff0)/w(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2, r_e)))

    do i = ixBmax2,ixBmin2,-1
      w(:,i, r_e) = (x(:,i+2,2)-x(:,i,2))*g0*w(:,i+1,rho_)*Gamma_dep(:,&
         i+1)/fld_lambda(:,i+1) + w(:,i+2, r_e)
    enddo
    !-------------------------------------------------------------------------
=======
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
    double precision :: fld_R(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
    double precision :: fld_kappa(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
    double precision :: fld_lambda(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
    double precision :: a(1:nw), b(1:nw), c(1:nw)
    integer :: i,j

    select case (iB)

    case(3)

      call fld_get_fluxlimiter(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixGmin1+2,ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_lambda, fld_R)
      call fld_get_opacity(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixGmin1+2,&
         ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_kappa)


      do i = ixBmin2,ixBmax2
        w(:,i, rho_) = p_bound*dexp(-x(:,i,2)/heff0)/c_sound0**2
        w(:,i, mom(1)) = zero
        velocity(:,i,2) = two*w(:,i+1,mom(2))/w(:,i+1,rho_) - w(:,i+2,&
           mom(2))/w(:,i+2,rho_)
        w(:,i, mom(2)) = zero !velocity(:,i,2)*w(:,i, rho_)
        w(:,i, e_) = p_bound*dexp(-x(:,i,2)/heff0)/(hd_gamma-one)
        w(:,i, r_e) = 3.d0*Gamma/(one-Gamma)*p_bound*exp(-x(:,i,2)/heff0)
      enddo
>>>>>>> 924324e1c3d570855f9b58bd582c5da8c596bec5


<<<<<<< HEAD
  case(4)

    ! !> Do quadratic interpolation
    ! do i = ixGmin1,ixGmax1
    !   a(:) = 0.5d0*w(i,ixGmax2-2, :) - w(i,ixGmax2-3, :) - 0.5d0*w(i,ixGmax2-4, :)
    !   b(:)  = -3.5d0*w(i,ixGmax2-2, :) + 6.d0*w(i,ixGmax2-3, :) - 2.5d0*w(i,ixGmax2-4, :)
    !   c(:)   = 6.d0*w(i,ixGmax2-2, :) - 8.d0*w(i,ixGmax2-3, :) + 3d0*w(i,ixGmax2-4, :)
    !
    !   w(i,ixGmax2-1, :) = a(:) + b(:) + c(:)
    !   w(i,ixGmax2, :) = c(:)
    ! enddo
    ! w(:,ixBmin2, :) = 2*w(:,ixBmin2-1, :) - w(:,ixBmin2-2, :)
    ! w(:,ixBmax2, :) = 2*w(:,ixBmax2-1, :) - w(:,ixBmax2-2, :)


    ! !> loglineair interpolation
    ! do i = ixGmin1,ixGmax1
    !   b(:) = dlog(w(i,ixGmax2-3,:)/w(i,ixGmax2-2,:))/dlog(x(i,ixGmax2-3,2)/x(i,ixGmax2-2,2))
    !   a(:)  = w(i,ixGmax2-2,:)*x(i,ixGmax2-2,2)**(-b(:))
    !
    !   !> Density
    !   w(i,ixGmax2-1, rho_) = a(rho_)*x(i,ixGmax2-1,2)**b(rho_)
    !   w(i,ixGmax2, rho_) = a(rho_)*x(i,ixGmax2,2)**b(rho_)
    !
    !   !> Gas Energy
    !   w(i,ixGmax2-1, e_) = a(e_)*x(i,ixGmax2-1,2)**b(e_)
    !   w(i,ixGmax2, e_) = a(e_)*x(i,ixGmax2,2)**b(e_)
    !
    !   !> Radiation Energy
    !   w(i,ixGmax2-1, r_e) = max(a(r_e)*x(i,ixGmax2-1,2)**b(r_e),zero)
    !   w(i,ixGmax2, r_e) = max(a(r_e)*x(i,ixGmax2,2)**b(r_e),zero)
    ! enddo

    ! !> exponential interpolation
    ! do i = ixGmin1,ixGmax1
    !   b(:) = dlog(w(i,ixGmax2-2,:)/w(i,ixGmax2-3,:))*(x(i,ixGmax2-2,2)/x(i,ixGmax2-3,2))
    !   a(:)  = w(i,ixGmax2-2,:)*dexp(x(i,ixGmax2-2,2)*(-b(:)))
    !
    !   !> Density
    !   w(i,ixGmax2-1, rho_) = a(rho_)*dexp(-x(i,ixGmax2-1,2)*b(rho_))
    !   w(i,ixGmax2, rho_) = a(rho_)*dexp(-x(i,ixGmax2,2)*b(rho_))
    !
    !   !> Gas Energy
    !   w(i,ixGmax2-1, e_) = a(e_)*dexp(-x(i,ixGmax2-1,2)*b(e_))
    !   w(i,ixGmax2, e_) = a(e_)*dexp(-x(i,ixGmax2,2)*b(e_))
    !
    !   !> Radiation Energy
    !   w(i,ixGmax2-1, r_e) = max(a(r_e)*dexp(-x(i,ixGmax2-1,2)*b(r_e)),zero)
    !   w(i,ixGmax2, r_e) = max(a(r_e)*dexp(-x(i,ixGmax2,2)*b(r_e)),zero)
    ! enddo

    !> Conservation law
    call fld_get_fluxlimiter(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixGmin1+2,&
       ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_lambda, fld_R)
    call fld_get_opacity(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixGmin1+2,&
       ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_kappa)

    do i = ixBmin2-1, ixBmax2-1
      w(:,i+1,r_e) = ((w(:,i-1,mom(2))*w(:,i-1,r_e)/w(:,i+1,&
         rho_)- fld_lambda(:,ixBmin2)*c_light0/(w(:,i-1,rho_)*fld_kappa(:,&
         ixBmin2)) * (w(:,i-2,r_e) - w(:,i,r_e))/abs(x(:,i+1,2)- x(:,i-1,&
         2))- w(:,i,mom(2))*w(:,i,r_e)/w(:,i,rho_))* w(:,i,rho_)*fld_kappa(:,&
         ixBmin2)/(fld_lambda(:,ixBmin2)*c_light0)*abs(x(:,i+1,2)- x(:,i-1,&
         2)))+ w(:,i-1,r_e)
       do j = ixGmin2,ixGmax2
         w(j,i+1,r_e) = min(w(j,i+1,r_e), w(j,i,r_e))
       enddo
    enddo

  case default
    call mpistop("BC not specified")
  end select
end subroutine special_bound

=======
      ! !> Do quadratic interpolation
      ! do i = ixGmin1,ixGmax1
      !   a(:) = 0.5d0*w(i,ixGmax2-2, :) - w(i,ixGmax2-3, :) - 0.5d0*w(i,ixGmax2-4, :)
      !   b(:)  = -3.5d0*w(i,ixGmax2-2, :) + 6.d0*w(i,ixGmax2-3, :) - 2.5d0*w(i,ixGmax2-4, :)
      !   c(:)   = 6.d0*w(i,ixGmax2-2, :) - 8.d0*w(i,ixGmax2-3, :) + 3d0*w(i,ixGmax2-4, :)
      !
      !   w(i,ixGmax2-1, :) = a(:) + b(:) + c(:)
      !   w(i,ixGmax2, :) = c(:)
      ! enddo
      ! w(:,ixBmin2, :) = 2*w(:,ixBmin2-1, :) - w(:,ixBmin2-2, :)
      ! w(:,ixBmax2, :) = 2*w(:,ixBmax2-1, :) - w(:,ixBmax2-2, :)


      ! !> loglineair interpolation
      ! do i = ixGmin1,ixGmax1
      !   b(:) = dlog(w(i,ixGmax2-3,:)/w(i,ixGmax2-2,:))/dlog(x(i,ixGmax2-3,2)/x(i,ixGmax2-2,2))
      !   a(:)  = w(i,ixGmax2-2,:)*x(i,ixGmax2-2,2)**(-b(:))
      !
      !   !> Density
      !   w(i,ixGmax2-1, rho_) = a(rho_)*x(i,ixGmax2-1,2)**b(rho_)
      !   w(i,ixGmax2, rho_) = a(rho_)*x(i,ixGmax2,2)**b(rho_)
      !
      !   !> Gas Energy
      !   w(i,ixGmax2-1, e_) = a(e_)*x(i,ixGmax2-1,2)**b(e_)
      !   w(i,ixGmax2, e_) = a(e_)*x(i,ixGmax2,2)**b(e_)
      !
      !   !> Radiation Energy
      !   w(i,ixGmax2-1, r_e) = max(a(r_e)*x(i,ixGmax2-1,2)**b(r_e),zero)
      !   w(i,ixGmax2, r_e) = max(a(r_e)*x(i,ixGmax2,2)**b(r_e),zero)
      ! enddo

      ! !> exponential interpolation
      ! do i = ixGmin1,ixGmax1
      !   b(:) = dlog(w(i,ixGmax2-2,:)/w(i,ixGmax2-3,:))*(x(i,ixGmax2-2,2)/x(i,ixGmax2-3,2))
      !   a(:)  = w(i,ixGmax2-2,:)*dexp(x(i,ixGmax2-2,2)*(-b(:)))
      !
      !   !> Density
      !   w(i,ixGmax2-1, rho_) = a(rho_)*dexp(-x(i,ixGmax2-1,2)*b(rho_))
      !   w(i,ixGmax2, rho_) = a(rho_)*dexp(-x(i,ixGmax2,2)*b(rho_))
      !
      !   !> Gas Energy
      !   w(i,ixGmax2-1, e_) = a(e_)*dexp(-x(i,ixGmax2-1,2)*b(e_))
      !   w(i,ixGmax2, e_) = a(e_)*dexp(-x(i,ixGmax2,2)*b(e_))
      !
      !   !> Radiation Energy
      !   w(i,ixGmax2-1, r_e) = max(a(r_e)*dexp(-x(i,ixGmax2-1,2)*b(r_e)),zero)
      !   w(i,ixGmax2, r_e) = max(a(r_e)*dexp(-x(i,ixGmax2,2)*b(r_e)),zero)
      ! enddo

      !> Conservation law
      call fld_get_fluxlimiter(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixGmin1+2,ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_lambda, fld_R)
      call fld_get_opacity(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixGmin1+2,&
         ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_kappa)

      do i = ixBmin2-1, ixBmax2-1
        ! w(:,i+1,r_e) = ((w(:,i-1,mom(2))*w(:,i-1,r_e)/w(:,i+1,rho_)&
        !  - fld_lambda(:,i-1)*c_light0/(w(:,i-1,rho_)*fld_kappa(:,i-1)) &
        !  * (w(:,i-2,r_e) - w(:,i,r_e))/abs(x(:,i+1,2)- x(:,i-1,2))&
        !  - w(:,i,mom(2))*w(:,i,r_e)/w(:,i,rho_))&
        !  * w(:,i,rho_)*fld_kappa(:,i-1)/(fld_lambda(:,i)*c_light0)*abs(x(:,i+1,2)- x(:,i-1,2)))&
        !  + w(:,i-1,r_e)
        w(:,i+1,r_e) = ((w(:,i-1,mom(2))*w(:,i-1,r_e)/w(:,i+1,&
           rho_)- fld_lambda(:,ixBmin2)*c_light0/(w(:,i-1,rho_)*fld_kappa(:,&
           ixBmin2)) * (w(:,i-2,r_e) - w(:,i,r_e))/abs(x(:,i+1,2)- x(:,i-1,&
           2))- w(:,i,mom(2))*w(:,i,r_e)/w(:,i,rho_))* w(:,i,rho_)*fld_kappa(:,&
           ixBmin2)/(fld_lambda(:,ixBmin2)*c_light0)*abs(x(:,i+1,2)- x(:,i-1,&
           2)))+ w(:,i-1,r_e)
         do j = ixGmin2,ixGmax2
           w(j,i+1,r_e) = min(w(j,i+1,r_e), w(j,i,r_e))
         enddo
      enddo

      print*," asdfafsfgvs"

      do i = ixGmin2, ixGmax2
        print*,  fld_kappa(5,i)/(unit_time*unit_velocity*unit_density)
      enddo
>>>>>>> 924324e1c3d570855f9b58bd582c5da8c596bec5

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

<<<<<<< HEAD
=======
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
         ixImin2:ixImax2)/(hd_gamma - one) + half*(w(ixImin1:ixImax1,&
         ixImin2:ixImax2,mom(1))**two + w(ixImin1:ixImax1,ixImin2:ixImax2,&
         mom(2))**two)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

    end subroutine constant_e

  !==========================================================================================

  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
>>>>>>> 924324e1c3d570855f9b58bd582c5da8c596bec5
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision :: pressure(ixImin1:ixImax1,ixImin2:ixImax2)

<<<<<<< HEAD
    pressure(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)*c_sound0**2
    w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = pressure(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(hd_gamma - one) + half*(w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(1))**two + w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(2))**two)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

  end subroutine constant_e
=======
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       2) = -6.67e-8*M_star/R_star**2*(unit_time**2/unit_length)
>>>>>>> 924324e1c3d570855f9b58bd582c5da8c596bec5

!==========================================================================================

!> Calculate gravitational acceleration in each dimension
subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
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
  double precision                   :: rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1:ndim), rad_pressure(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
      fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2), fld_R(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2), fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
  double precision                   :: g_rad(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2), D(ixImin1:ixImax1,&
     ixImin2:ixImax2,1:ndim)
  integer                            :: idim

  call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, rad_flux)
  call fld_get_radpress(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, rad_pressure)
  call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, fld_lambda, fld_R)
  call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, fld_kappa)
  call fld_get_diffcoef(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, D)

  do idim = 1,ndim
    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim) = fld_kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idim)/c_light0
  enddo
  big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)/(6.67e-8*M_star/R_star**2*(unit_time**2/unit_length))

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)*(unit_pressure*unit_velocity)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)*(unit_pressure*unit_velocity)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=rad_pressure(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*unit_pressure
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4)=fld_lambda(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5)=fld_R(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6)=g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)*unit_length/(unit_time**2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+7)=g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)*unit_length/(unit_time**2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+8)=big_gamma(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+9)=D(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+10)=D(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     2)

end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'F1 F2 RP lam fld_R ar1 ar2 Gam D1 D2'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
