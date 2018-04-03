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

    call initglobaldata_usr()

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

  Flux0 = L_star/(4*dpi*R_star**2)*(unit_density*unit_length**&
     3/(unit_velocity*unit_pressure))
  T_star0 = T_star/unit_temperature
  c_sound0 = dsqrt((1.38d-16*T_star0/(0.6*mp_cgs))/unit_velocity)
  kappa0 = 0.34/unit_length**2
  c_light0 = const_c/unit_velocity
  g0 = 6.67e-8*M_star/R_star**2*(unit_density*unit_length/unit_pressure)


end subroutine initglobaldata_usr

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_hd_phys, only: hd_get_pthermal

    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
       ixmin2,ixmax1,ixmax2
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)

    ! Set initial values for w
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = one
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = &
       one/(one-3.d0/5.d0)*(c_sound0**2+(kappa0*Flux0/c_light0-&
       g0)*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,2)-x(1,1,2)))
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = &
       c_sound0*T_star0**4-3*kappa0*Flux0/c_light0*(x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,2)-x(1,1,2))

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

  subroutine special_bound(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,ixBmin2,&
     ixBmax1,ixBmax2,iB,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,ixBmin2,&
       ixBmax1,ixBmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision :: e_inflo

    select case (iB)

    case(3)
      w(:,ixBmax2, rho_) = one*perturbation(ixGmin2,ixGmax2,2)
      w(:,ixBmax2, mom(1)) = zero
      w(:,ixBmax2, mom(1)) = w(:,ixBmax2+1, mom(1))*perturbation(ixGmin2,&
         ixGmax2,6)
      w(:,ixBmax2, e_) = one/(one-3.d0/5.d0)*c_sound0**2*perturbation(ixGmin2,&
         ixGmax2,4)
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
subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
     ixImin2:ixImax2,ndim)

  ! phi = -GM/R

  gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
  gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
     2) = 6.67e-8*M_star/(R_star+x(ixImin1:ixImax1,ixImin2:ixImax2,&
     2))*(unit_density/unit_pressure)

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

  call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, rad_flux, rad_pressure)
  call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, fld_lambda, fld_R)

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)
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
