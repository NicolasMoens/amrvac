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

double precision :: Gamma, c_sound, Flux, g_eff, g_grav, H_eff, kappa, c_light

double precision :: Flux0, c_sound0, T_star0, kappa0, T_bound0
double precision :: c_light0, g0, geff0, heff0, Gamma_edd
double precision :: L_star0, R_star0, M_star0
double precision :: tau_bound,  P_bound, rho_bound, T_bound

double precision :: lower_bc_rho(2), lower_bc_e(2), lower_bc_re(2)

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
  tau_bound = 250.d0

  call initglobaldata_usr

  ! Initialize units
  usr_set_parameters => initglobaldata_usr

  ! A routine for initial conditions is always required
  usr_init_one_grid => initial_conditions

  ! Special Boundary conditions
  usr_special_bc => special_bound

  ! Graviatational field
  usr_gravity => set_gravitation_field

  ! Output routines
  usr_aux_output    => specialvar_output
  usr_add_aux_names => specialvarnames_output

  ! Active the physics module
  call hd_activate()

end subroutine usr_init

!===============================================================================

subroutine initglobaldata_usr
use mod_global_parameters

T_bound = (one +(3.d0/4.d0*tau_bound)**(1.d0/4.d0))*T_star
c_sound =  dsqrt((1.38d-16*T_bound/(0.6*mp_cgs)))

g_grav = 6.67e-8*M_star/R_star**two

Flux =  L_star/(4*dpi*R_star**2)
kappa = 0.34d0
c_light = const_c
Gamma_edd = (kappa*Flux)/(c_light*g_grav)
g_eff = g_grav*(one - Gamma_edd)
H_eff = c_sound**2/g_eff

p_bound = g_eff*tau_bound/kappa
rho_bound = p_bound/c_sound**two

unit_length        = H_eff
unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*mp_cgs)
unit_temperature   = T_bound

L_star0 = L_star*unit_time/(unit_pressure*unit_length**3.d0)
R_star0 = R_star/unit_length
M_star0 = M_star/(unit_density*unit_length**3.d0)

Flux0 = L_star0/(4*dpi*R_star0**2)
T_star0 = T_star/unit_temperature
T_bound0 = T_bound/unit_temperature
c_sound0 = dsqrt((1.38d-16*T_star/(0.6*mp_cgs)))/unit_velocity

kappa0 = fld_kappa0
c_light0 = const_c/unit_velocity
g0 = 6.67e-8*M_star/R_star**2*(unit_density*unit_length/unit_pressure)
geff0 = g0*(one - (kappa0*Flux0)/(c_light0*g0))
heff0 = c_sound0**2/geff0
Gamma = (kappa0*Flux0)/(c_light0*g0)

P_bound = one
rho_bound = one



end subroutine initglobaldata_usr

!==========================================================================================

!> A routine for specifying initial conditions
subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
   ixmax1,ixmax2, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables

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

  double precision :: a,b,c

  double precision :: opt_depth(ixGmin2:ixGmax2)
  double precision :: temp_init(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
  double precision :: temperature(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
  double precision :: k1(ixGmin1:ixGmax1,ixGmin2:ixGmax2), k2(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)
  double precision :: k3(ixGmin1:ixGmax1,ixGmin2:ixGmax2), k4(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)

  integer :: i

<<<<<<< HEAD
  amplitude = 0.1d0 !5.d-1  !1.d-5 !3.d-2
=======
  amplitude = zero !5.d-2 !1.d-5 !3.d-2
>>>>>>> af21b50c2b25ea6bf0bf0bd00635ccf77318ac13

  pressure(:,ixGmin2) = p_bound

  ! Set pressure profile using RK
  a = Flux0*kappa0/(4.d0/3.d0*fld_sigma_0*geff0)
  b = T_bound0**4.d0 - a*p_bound
  c = -geff0*mp_cgs *0.6d0/kB_cgs * unit_pressure/(&
     unit_temperature*unit_density)

  pressure(:,ixGmin2) = p_bound

  do i=ixGmin2,ixGmax2-1
    k1(:,i) = (x(:,1+1,2) - x(:,1,2))*c*pressure(:,i)/((a*pressure(:,&
       i)+b)**(1.d0/4.d0))
    k2(:,i) = (x(:,1+1,2) - x(:,1,2))*c*(pressure(:,i)+half*k1(:,&
       i))/((a*(pressure(:,i)+half*k1(:,i))+b)**(1.d0/4.d0))
    k3(:,i) = (x(:,1+1,2) - x(:,1,2))*c*(pressure(:,i)+half*k2(:,&
       i))/((a*(pressure(:,i)+half*k2(:,i))+b)**(1.d0/4.d0))
    k4(:,i) = (x(:,1+1,2) - x(:,1,2))*c*(pressure(:,i)+k3(:,&
       i))/((a*(pressure(:,i)+k3(:,i))+b)**(1.d0/4.d0))

    pressure(:,i+1) = pressure(:,i) + one/6.d0 * (k1(:,i) + two*k2(:,&
       i) + two*k3(:,i) + k4(:,i))
  enddo

  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)/(hd_gamma - one)
  density(ixGmin1:ixGmax1,ixGmin2:ixGmax2) = pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)/((a*pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)+b)**(1.d0/4.d0))*mp_cgs*0.6d0/kB_cgs * &
     unit_pressure/(unit_temperature*unit_density)
  temp_init(ixGmin1:ixGmax1,ixGmin2:ixGmax2) = (a*pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2) + b)**(1.d0/4.d0)

  do i = ixGmin2,ixGmax2
    opt_depth(i) =  tau_bound - fld_kappa0*sum(density(5,:i))*(x(5,2,2)-x(5,1,&
       2))
  enddo

  ! Set initial values for w
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = density(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)/(hd_gamma - one)
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = &
     4.d0*fld_sigma_0/fld_speedofligt_0*(a*pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2) + b)


  lower_bc_rho(:) = w(5,1:2, rho_)
  lower_bc_e(:) = w(5,1:2, e_)
  lower_bc_re(:) = w(5,1:2, r_e)


  !> perturb rho
  call RANDOM_NUMBER(pert)
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
      rho_)*(one + amplitude*pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2))

  ! !> Write energy to file
  ! open(1,file='initial_cond')
  ! do i=ixGmin2,ixGmax2
  ! write(1,222) x(5,i,2), w(5,i,rho_), w(5,i,e_), w(5,i,r_e), temp_init(5,i), opt_depth(i)
  ! enddo
  ! 222 format(6e15.5E3)
  ! stop

end subroutine initial_conditions

!==========================================================================================

! Extra routines can be placed here
! ...

!==========================================================================================

subroutine special_bound(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,ixBmin2,&
   ixBmax1,ixBmax2,iB,w,x)

  use mod_global_parameters
  use mod_variables
  use mod_constants
  use mod_fld
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
  double precision :: kb0, mp0

  double precision :: num_flux, dens, op, dy
  double precision :: lb, ub, cb
  double precision :: e_ip, e_i, e_im

  integer :: i,j

  select case (iB)

  case(1)
    do i = ixGmin2, ixGmax2
      w(ixBmin1:ixBmax1,i,:) = w(ixGmax1-3:ixGmax1-2,i,:)
    enddo

  case(2)
    do i = ixGmin2, ixGmax2
      w(ixBmin1:ixBmax1,i,:) = w(ixGmin1+2:ixGmin1+3,i,:)
    enddo

  case(3)

    w(:,2, rho_) = lower_bc_rho(2)
    w(:,2,mom(2)) = w(:,3,mom(2))
    w(:,2, r_e) = w(:,3, r_e) - 3.0*w(:,2,&
        rho_)*fld_kappa0*Flux0/fld_speedofligt_0*(x(:,2,2)-x(:,3,2))
    pressure(:,2) = 1.38d-16/(0.6*mp_cgs)*w(:,2,&
        rho_)/unit_pressure*(unit_temperature*unit_density)*(&
       fld_speedofligt_0/(4.d0*fld_sigma_0)*w(:,2, r_e))**0.25d0
    w(:,2, e_) = pressure(:,2)/(hd_gamma - one) + (w(:,2, mom(2))*w(:,2,&
        mom(2))/(2*w(:,2, rho_)))



    w(:,1, rho_) = lower_bc_rho(1)
    w(:,1,mom(2)) = w(:,3,mom(2))
    w(:,1, r_e) = w(:,2, r_e) - 3.0*w(:,1,&
        rho_)*fld_kappa0*Flux0/fld_speedofligt_0*(x(:,1,2)-x(:,2,2))
    pressure(:,1) = 1.38d-16/(0.6*mp_cgs)*w(:,1,&
        rho_)/unit_pressure*(unit_temperature*unit_density)*(&
       fld_speedofligt_0/(4.d0*fld_sigma_0)*w(:,1, r_e))**0.25d0
    w(:,1, e_) = pressure(:,1)/(hd_gamma - one) + (w(:,1, mom(2))*w(:,1,&
        mom(2))/(2*w(:,1, rho_)))

  case(4)

    do i = ixGmin1,ixGmax1
      w(i, ixBmin2, rho_) = min(2*w(i, ixBmin2-1, rho_) - w(i, ixBmin2-2,&
          rho_),w(i, ixBmin2-1, rho_))
      w(i, ixBmax2, rho_) = min(2*w(i, ixBmax2-1, rho_) - w(i, ixBmax2-2,&
          rho_),w(i, ixBmax2-1, rho_))

      w(i, ixBmin2, mom(:)) = w(i, ixBmin2-1, mom(:))
      w(i, ixBmax2, mom(:)) = w(i, ixBmax2-1, mom(:))
    enddo

    !> Linear interpolation
    do i = ixGmin1,ixGmax1
      ! w(i, ixBmin2, r_e) = -(w(i, ixBmin2-2, r_e) - w(i, ixBmin2-1, r_e))/(x(i,ixBmin2-2,2) - x(i,ixBmin2-1,2)) &
      ! * (x(i,ixBmin2-1,2) - x(i,ixBmin2,2)) + w(i, ixBmin2-1, r_e)
      ! w(i, ixBmax2, r_e) = -(w(i, ixBmax2-2, r_e) - w(i, ixBmax2-1, r_e))/(x(i,ixBmax2-2,2) - x(i,ixBmax2-1,2)) &
      ! * (x(i,ixBmax2-1,2) - x(i,ixBmax2,2)) + w(i, ixBmax2-1, r_e)

      ! w(i, ixBmin2, r_e) = w(i, ixBmin2 - 2, rho_)/w(i, ixBmin2 - 1, rho_)&
      ! *(w(i, ixBmin2-1, r_e) - w(i, ixBmin2-3, r_e)) + w(i, ixBmin2 - 2, r_e)
      ! w(i, ixBmax2, r_e) = w(i, ixBmax2 - 2, rho_)/w(i, ixBmax2 - 1, rho_)&
      ! *(w(i, ixBmax2-1, r_e) - w(i, ixBmax2-3, r_e)) + w(i, ixBmax2 - 2, r_e)

      ! w(i, ixBmin2, r_e) = (w(i, ixBmin2-1, r_e) - w(i, ixBmin2-3, r_e))&
      ! *(w(i, ixBmin2-2, rho_)*w(i, ixBmin2-2, r_e))/(w(i, ixBmin2-1, rho_)*w(i, ixBmin2-1, r_e))&
      ! +w(i, ixBmin2 - 2, r_e)
      ! w(i, ixBmax2, r_e) = (w(i, ixBmax2-1, r_e) - w(i, ixBmax2-3, r_e))&
      ! *(w(i, ixBmax2-2, rho_)*w(i, ixBmax2-2, r_e))/(w(i, ixBmax2-1, rho_)*w(i, ixBmax2-1, r_e))&
      ! +w(i, ixBmax2 - 2, r_e)

      ! w(i, ixBmin2, r_e) = w(i, ixBmin2-1, r_e)
      ! w(i, ixBmax2, r_e) = w(i, ixBmax2-1, r_e)


      do j = ixBmin2, ixBmax2
        ! w(i,j,r_e) = min(w(i,j,r_e),w(i,j-1,r_e))
        ! w(i,j,r_e) = max(w(i,j,r_e), zero)

        lb = zero
        ub = two*w(i,j-1,r_e)
        cb = (lb + ub)/two

        dens = w(i,j-1,rho_)
        op = kappa0
        dy = x(i,j,2) - x(i,j-2,2)

        e_i = w(i,j-1,r_e)
        e_im = w(i,j-2,r_e)

        do while (abs(num_flux - Flux0)/Flux0 .gt. 1.d-5)
          if (Numerical_flux(lb, e_i, e_im, op, dens, dy)*Numerical_flux(ub,&
              e_i, e_im, op, dens, dy) .gt. zero) then
            e_ip = ub
            num_flux = Flux0
          else

            if (Numerical_flux(lb, e_i, e_im, op, dens, dy)*Numerical_flux(cb,&
                e_i, e_im, op, dens, dy) .lt. zero) then
              ub = cb
            else
              lb = cb
            endif
          cb = (lb + ub)/two
          e_ip = cb
          num_flux = Numerical_flux(e_ip, e_i, e_im, op, dens, dy)
          endif
        enddo
        w(i,j,r_e) = e_ip
      enddo
    enddo

    ! do i = ixGmin1,ixGmax1
    !   print*, i, w(i, ixBmax2-5:,r_e)
    ! enddo


    !> Corners?
    w(ixGmin1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-3,ixBmin2:ixBmax2,r_e)
    w(ixGmin1+1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-2,ixBmin2:ixBmax2,r_e)
    w(ixGmax1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+3,ixBmin2:ixBmax2,r_e)
    w(ixGmax1-1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+2,ixBmin2:ixBmax2,r_e)

  case default
    call mpistop("BC not specified")
  end select
end subroutine special_bound




function Numerical_flux(e_ip, e_i, e_im, op, dens, dy) result(f)
  use mod_global_parameters

  double precision, intent(in) :: e_ip, e_i, e_im
  double precision, intent(in) :: op,dens, dy
  double precision :: f

  f = -fld_speedofligt_0/(op*dens)*(2.d0 + (e_ip - &
     e_im)/(2.d0*dy*op*dens*e_i))/((2.d0*6.d0*dy)/(e_ip - e_im) + &
     (3.d0)/(op*dens*e_i) + (e_ip - e_im)/(op*dens*e_i)**2.d0)

end function Numerical_flux


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
     ixImin2:ixImax2,1:ndim), Temp(ixImin1:ixImax1,ixImin2:ixImax2)
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
  call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,Temp)

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
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+11)= fld_kappa(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+12)=Temp(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'F1 F2 P_rad lamnda fld_R ar1 ar2 Gamma D1 D2 kappa Temp'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
