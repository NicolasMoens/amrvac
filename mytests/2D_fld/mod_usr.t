!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    usr_special_bc => special_bound

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

  !Fix dimensionless stuff here
  unit_length        = 1.d0                                         ! cm
  unit_numberdensity = 1.d3                                         ! cm^-3
  unit_temperature   = 1.d0                                         ! K


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
    double precision :: fld_boltzman_cgs = 5.67036713d-8

    integer :: i

    ! Set initial values for w
    w(ix^S, rho_) = 1.d0
    w(ix^S, mom(1)) = 0.d0
    w(ix^S, mom(2)) = 0.d0

    !> Try some ~1/r^2 init condition for energy
    w(ix^S, e_) = 1.d0
    w(ix^S, e_) = w(ix^S, e_) + 1.d0/x(ix^S,2)**2

    ! !> Radiative Equilibrium, heating = cooling
    call hd_get_pthermal(w,x,ix^L,ix^L,temperature)

    temperature(ix^S) = temperature(ix^S)/w(ix^S,iw_rho)*mp_cgs*fld_mu/kb_cgs&
    *unit_length**2/(unit_time**2 * unit_temperature)

    w(ix^S,r_e) = 4*unit_velocity/const_c* maxval(temperature(ix^S))**4*&
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

    select case (iB)

    ! case(1)
    !   w(:,ixBmax2,:) = w(:,ixBmax2+1,:)
    !   w(:,ixBmin2,:) = w(:,ixBmax2,:)
    !
    ! case(2)
    !   w(:,ixBmin1,:) = w(:,ixBmin1-1,:)
    !   w(:,ixBmax1,:) = w(:,ixBmin1,:)

    ! case(3)
    !   w(:,ixBmin2,rho_) = 1.d0
    !   w(:,ixBmin2,mom(1)) = 0.d0
    !   w(:,ixBmin2,mom(2)) = 0.d0
    !   w(:,ixBmin2,e_) = 1.d0
    !   w(:,ixBmin2,r_e) = 1.d0
    !
    !   w(:,ixBmax2,:) = w(:,ixBmin2,:)

    case(3)
      w(:,ixBmax2, rho_) =  w(:,ixBmax2+1, rho_)
      w(:,ixBmax2, mom(1)) = w(:,ixBmax2+1, mom(1))

      where (w(:,ixBmax2+1 , mom(2)) > zero)
        w(:,ixBmax2, mom(2)) = w(:,ixBmax2+1 , mom(2))
      elsewhere
        w(:,ixBmax2, mom(2)) = zero
      end where

      w(:,ixBmax2, e_) = w(:,ixBmax2+1, e_)
      w(:,ixBmax2, r_e) = 1.d0/(10*x(1,1,2))**2 !w(:,ixBmax2+1, r_e)
      w(:,ixBmin2,:) = w(:,ixBmax2,:)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==========================================================================================

end module mod_usr

!==========================================================================================
