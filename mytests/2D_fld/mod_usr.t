!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

  ! Custom variables can be defined here
  double precision :: fld_boltzman_cgs = 5.67036713d-8
  double precision :: d_inflo = 1d3
  double precision :: T_eff = 5d3

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    !Fix dimensionless stuff here
    unit_length        = 1.d0                                         ! cm
    unit_numberdensity = 1.d3                                         ! cm^-3
    unit_temperature   = 1.d0                                         ! K

    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    usr_special_bc => special_bound
    !usr_gravity => set_gravitation_field

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call hd_activate()

    print*, 'unit_time', unit_time
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, '================================================================'

  end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
  use mod_global_parameters






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

    double precision :: temperature(ixG^S)
    double precision :: int_energy(ixG^S)

    integer :: i

    ! Set initial values for w
    w(ixG^S, rho_) = d_inflo
    w(ixG^S, mom(1)) = 0.d0
    w(ixG^S, mom(2)) = 0.d0

    !> Try some ~1/r^2 init condition for energy
    int_energy = T_eff * one/(hd_gamma-one) * w(ixG^S, rho_)*kb_cgs/(mp_cgs*fld_mu)*unit_time**2/(unit_length**2)
    w(ixG^S, e_) = int_energy
    w(ix^S, e_) = w(ix^S, e_) + 1.d-2/x(ix^S,2)**2

    ! !> Radiative Equilibrium, heating = cooling
    ! call hd_get_pthermal(w,x,ixG^L,ixG^L,temperature)
    !
    ! temperature(ixG^S) = temperature(ixG^S)/w(ixG^S,iw_rho)*mp_cgs*fld_mu/kb_cgs&
    ! *unit_length**2/(unit_time**2 * unit_temperature)
    !
    ! w(ixG^S,r_e) = 4*unit_velocity/const_c* temperature(ixG^S)**4*&
    ! fld_boltzman_cgs*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density)

    !> Radiative Equilibrium, heating = cooling
    w(ixG^S,r_e) = 4*unit_velocity/const_c* (T_eff/unit_temperature)**4*&
    fld_boltzman_cgs*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density)

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
      w(:,ixBmax2, rho_) =  d_inflo !w(:,ixBmax2+1, rho_)
      w(:,ixBmax2, mom(1)) = w(:,ixBmax2+1, mom(1))

      where (w(:,ixBmax2+1 , mom(2)) > zero)
        w(:,ixBmax2, mom(2)) = w(:,ixBmax2+1 , mom(2))
      elsewhere
        w(:,ixBmax2, mom(2)) = zero
      end where

      e_inflo = T_eff * one/(hd_gamma-one) * d_inflo *kb_cgs/(mp_cgs*fld_mu)*unit_time**2/(unit_length**2)
      w(:,ixBmax2, e_) = e_inflo !1.d0 !w(:,ixBmax2+1, e_)
      w(:,ixBmax2, r_e) = 4*unit_velocity/const_c* (T_eff/unit_temperature)**4*&
      fld_boltzman_cgs*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density) !w(:,ixBmax2+1, r_e)
      w(:,ixBmin2,:) = w(:,ixBmax2,:)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==========================================================================================

!> Calculate gravitational acceleration in each dimension
subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixI^L, ixO^L
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: wCT(ixI^S,1:nw)
  double precision, intent(out)   :: gravity_field(ixI^S,ndim)

  double precision :: mass, distance, cavendish

  mass = 1d35*unit_density*unit_length**3
  distance = 1d3*unit_length
  cavendish = 6.67d-8**unit_density*unit_time

  gravity_field(ixI^S,1) = zero
  gravity_field(ixI^S,2) = -cavendish*mass/((distance + x(ixI^S,2))**2)

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

  varnames = 'F1 F2 RP L R'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
