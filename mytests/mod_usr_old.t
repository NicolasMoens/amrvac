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

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    usr_special_bc => special_bound
    print*, "pointing toward special bound"
    print*,"##################################################################################################################################################################################################################################"

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

    ! double precision :: p(ixG^L)
    ! integer :: i

    ! Set initial values for w
    w(ix^S, rho_) = 1.
    w(ix^S, mom(1)) = 0.5
    w(ix^S, mom(2)) = 1.5
    w(ix^S, e_) = 20.

    w(ix^S,r_e) = 0.3
    w(ix^S,r_f(1)) = 2.0
    w(ix^S,r_f(2)) = 3.0

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

    print*, "####################################################################"
    print*,"iB" , iB

    do ix1 = 1,2
      w(ix1,ixImin2:ixImax2, rho_) =  1.d0
      w(ix1,ixImin2:ixImax2, mom(1)) = 0.5d0
      w(ix1,ixImin2:ixImax2, mom(2)) = 1.5d0
      w(ix1,ixImin2:ixImax2, e_) = 20.d0

      w(ix1,ixImin2:ixImax2, r_e) = 0.3d0
      w(ix1,ixImin2:ixImax2, r_f(1)) = 2.d0
      w(ix1,ixImin2:ixImax2, r_f(2)) = 3.d0
    end do

  end subroutine special_bound

!==========================================================================================

end module mod_usr

!==========================================================================================
