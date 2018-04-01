  !> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

  double precision :: M_sun = 1.99d33
  double precision :: L_sun = 3.99d33
  double precision :: R_sun = 6.96d10

  double precision :: M_star
  double precision :: L_star
  double precision :: T_star
  double precision :: R_star
  double precision :: rho_star

  double precision :: Flux0, c_sound0, T_star0, kappa0, c_light0, g0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_constants

    call set_coordinate_system("Cartesian_2D")

    M_star = 1*M_sun
    L_star = 1*L_sun
    T_star = 6000
    R_star = 1*R_sun
    rho_star = 5.d2

    !Fix dimensionless stuff here
    unit_length        = R_sun
    unit_numberdensity = rho_star/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature   = T_star

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Special Boundary conditions
    usr_special_bc => special_bound

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Output routines
    ! usr_aux_output    => specialvar_output
    ! usr_add_aux_names => specialvarnames_output

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

  Flux0 = L_star/(4*dpi*R_star**2)&
  *(unit_density*unit_length**3/(unit_velocity*unit_pressure))
  T_star0 = T_star/unit_temperature
  c_sound0 = dsqrt((1.38d-16*T_star0/(0.6*mp_cgs))/unit_velocity)
  kappa0 = 0.34/unit_length**2
  c_light0 = const_c/unit_velocity
  g0 = 6.67e-8*M_star/R_star**2&
  *(unit_density*unit_length/unit_pressure)


end subroutine initglobaldata_usr

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixG^L, ix^L, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_hd_phys, only: hd_get_pthermal

    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, ndim)
    double precision, intent(inout) :: w(ixG^S, nw)

    ! Set initial values for w
    w(ixG^S, rho_) = one
    w(ixG^S, mom(:)) = zero
    w(ixG^S, e_) = one/(one-3.d0/5.d0)*(c_sound0**2&
        +(kappa0*Flux0/c_light0-g0)*(x(ixG^S,2)-x(1,1,2)))
    w(ixG^S,r_e) =  c_sound0*T_star0**4&
        -3*kappa0*Flux0/c_light0*(x(ixG^S,2)-x(1,1,2))

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

  subroutine special_bound(qt,ixG^L,ixB^L,iB,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: e_inflo

    select case (iB)

    case(3)
      w(:,ixBmax2, rho_) = one*perturbation(ixGmin2,ixGmax2,2)
      w(:,ixBmax2, mom(1)) = zero
      w(:,ixBmax2, mom(1)) = w(:,ixBmax2+1, mom(1))*perturbation(ixGmin2,ixGmax2,6)
      w(:,ixBmax2, e_) = one/(one-3.d0/5.d0)*c_sound0**2*perturbation(ixGmin2,ixGmax2,4)
      w(:,ixBmax2, r_e) = c_sound0*T_star0**4*perturbation(ixGmin2,ixGmax2,3)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  function perturbation(ixOmin,ixOmax,n) result(pert_ampl)
    use mod_global_parameters

    integer, intent(in) :: ixOmin,ixOmax, n
    double precision :: pert_ampl(ixOmin:ixOmax)

    integer :: i

    do i =  ixOmin,ixOmax
      pert_ampl(i) = sin(two*dpi*n*(i-ixOmin)/ixOmax)
      pert_ampl(i) = pert_ampl(i)*sin((global_time/dt)*1.d-1)
    enddo

    pert_ampl = one + 1.d-2*pert_ampl

  end function perturbation

!==========================================================================================

!> Calculate gravitational acceleration in each dimension
subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixI^L, ixO^L
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: wCT(ixI^S,1:nw)
  double precision, intent(out)   :: gravity_field(ixI^S,ndim)

  ! phi = -GM/R

  gravity_field(ixI^S,1) = zero
  gravity_field(ixI^S,2) = 6.67e-8*M_star/(R_star+x(ixI^S,2))&
  *(unit_density/unit_pressure)

end subroutine set_gravitation_field


!==========================================================================================

subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
! corresponding normalization values (default value 1)
  use mod_global_parameters
  use mod_physics

  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision                   :: rad_flux(ixI^S,1:ndim), rad_pressure(ixI^S), fld_lambda(ixI^S), fld_R(ixI^S)

  call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, rad_pressure)
  call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)

  w(ixO^S,nw+1)=rad_flux(ixO^S,1)
  w(ixO^S,nw+2)=rad_flux(ixO^S,2)
  w(ixO^S,nw+3)=rad_pressure(ixO^S)
  w(ixO^S,nw+4)=fld_lambda(ixO^S)
  w(ixO^S,nw+5)=fld_R(ixO^S)

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
