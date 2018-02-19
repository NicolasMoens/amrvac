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

    call set_coordinate_system("Cartesian_1D")

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9 !cm-3,cm-3

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
  subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
       ixmin2,ixmax1,ixmax2
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)

    ! double precision :: p(ixG^L)
    ! integer :: i

    ! Set initial values for w
    w(ixmin1:ixmax1,ixmin2:ixmax2, rho_) = 1.
    w(ixmin1:ixmax1,ixmin2:ixmax2, mom(1)) = 0.5
    w(ixmin1:ixmax1,ixmin2:ixmax2, e_) = 20.

    w(ixmin1:ixmax1,ixmin2:ixmax2,r_e) = 0.3
    w(ixmin1:ixmax1,ixmin2:ixmax2,r_f(1)) = 2.0

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

  subroutine special_bound(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iB,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iB
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: ix1

    do ix1 = 1,2
      w(ix1, rho_) =  1.
      w(ix1, mom(1)) = 0.5
      w(ix1, e_) = 20.

      w(ix1, r_e) = 0.3
      w(ix1, r_f(1)) = 2.
    end do

  end subroutine special_bound

!==========================================================================================

end module mod_usr

!==========================================================================================
