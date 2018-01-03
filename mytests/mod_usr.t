!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions
    usr_special_bc => special_bound

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! ...

    ! Active the physics module
    call hd_activate()
  end subroutine usr_init

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixG^L, ix^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, ndim)
    double precision, intent(inout) :: w(ixG^S, nw)

    ! Set initial values for w
    w(ix^S, rho_) = 1./2*( 1 + x(ix^S, 1)**2 + x(ixG^S, 2)**2 )**2
    w(ix^S, mom(1)) = 1.
    w(ix^S, mom(2)) = 0.
    w(ix^S, e_) = 1.

    !w(ix^S,rad_e) = 2.
    !w(ix^S,rad_flux(1)) = 1./( 1 + x(ix^S, 1)**2 + x(ixG^S, 2)**2 )
    !w(ix^S,rad_flux(1)) = 1.

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

  subroutine special_bound(qt,ixI^L,ixO^L,iB,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ix1


    !---------------------------------------------

    ! Set boundary values for w

    do ix1 = ixImin1,ixImax1
      w(ix1,ixImin2:ixImax2, rho_) =  1./2.
      w(ix1,ixImin2:ixImax2, mom(1)) = 0.
      w(ix1,ixImin2:ixImax2, mom(2)) = 1.
      w(ix1,ixImin2:ixImax2, e_) = 1.

      !w(ix1,ixImin2:ixImax2, rad_flux(1)) = 0.
      !w(ix1,ixImin2:ixImax2, rad_flux(2)) = 1.
      !w(ix1,ixImin2:ixImax2, rad_e) = 2.
    end do


  end subroutine special_bound

!==========================================================================================

end module mod_usr

!==========================================================================================
