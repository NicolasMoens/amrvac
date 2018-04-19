  !> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld

  implicit none

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    use mod_constants

     double precision :: rho_0 = 1.d0
    ! double precision :: t_0 = 1.d-2
    ! double precision :: e_0 = 1.d0

    call set_coordinate_system("Cartesian_2D")

    ! !Fix dimensionless stuff here
    ! unit_length        = dsqrt(e_0/rho_0)*t_0                                        ! cm
    ! unit_numberdensity = rho_0/((1.d0+4.d0*He_abundance)*mp_cgs) !rho_0/(fld_mu*mp_cgs)                                      ! cm^-3
    ! unit_temperature   = e_0/(unit_numberdensity*(2.d0+3.d0*He_abundance)*kB_cgs) !e_0/(unit_numberdensity*hd_gamma*kB_cgs)                   ! K

    unit_velocity = const_c
    unit_numberdensity = rho_0/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_length = one/const_c

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Keep the radiative energy constant with internal bound
    usr_internal_bc => constant_var

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
    w(ixG^S, rho_) = 1.d2
    w(ixG^S, mom(:)) = zero
    w(ixG^S, e_) = 1.d0
    w(ixG^S,r_e) =  spotpattern(x,ixG^L,0.d0)


  end subroutine initial_conditions

  function spotpattern(x,ixG^L,t1) result(e0)
    use mod_global_parameters

    integer, intent(in) :: ixG^L
    double precision, intent(in) :: x(ixG^S, ndim), t1
    double precision :: e0(ixG^S)
    integer i,j

    do i = ixGmin1,ixGmax1
      do j = ixGmin2,ixGmax2
      e0(i,j) =  two + dexp(-8.d0 *dpi**2*t1)*sin(2*dpi*x(i,j,1))*sin(2*dpi*x(i,j,2))
      enddo
    enddo

  end function spotpattern

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

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

  subroutine constant_var(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    w(ixI^S,rho_) = 1.d2
    w(ixI^S,mom(:)) = zero
    w(ixI^S,e_) = 1.d0

  end subroutine constant_var

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
  double precision                   :: residual(ixI^S)
  ! double precision                   :: rad_flux(ixI^S,1:ndim), rad_pressure(ixI^S), fld_lambda(ixI^S), fld_R(ixI^S)

  ! call fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure)
  ! call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
  ! call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)

  residual(ixI^S) = spotpattern(x,ixI^L,global_time)
  residual(ixI^S) = (residual(ixI^S) - w(ixI^S,r_e))/residual(ixI^S)

  w(ixO^S,nw+1) = residual(ixO^S)

  ! w(ixO^S,nw+1)=rad_flux(ixO^S,1)
  ! w(ixO^S,nw+2)=rad_flux(ixO^S,2)
  ! w(ixO^S,nw+3)=rad_pressure(ixO^S)
  ! w(ixO^S,nw+4)=fld_lambda(ixO^S)
  ! w(ixO^S,nw+5)=fld_R(ixO^S)


end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'residual'
  ! varnames = 'F1 F2 RP lam fld_R'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
