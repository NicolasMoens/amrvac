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

     double precision :: rho_0 = one
     double precision :: t_0 = -dlog(1.d-1)/(8.d0*dpi**two)
     double precision :: e_0 = one

    call set_coordinate_system("Cartesian_2D")

    unit_velocity = const_c
    unit_numberdensity = rho_0/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_length = t_0*const_c

    ! unit_velocity = one
    ! unit_numberdensity = rho_0/((1.d0+4.d0*He_abundance)*mp_cgs)
    ! unit_length = one

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
    print*, t_0
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
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = 1.d2
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = one
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) =  spotpattern(x,ixGmin1,ixGmin2,&
       ixGmax1,ixGmax2,0.d0)
  end subroutine initial_conditions

  function spotpattern(x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,t1) result(e0)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, ndim),&
        t1
    double precision :: e0(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
    integer i,j

    do i = ixGmin1,ixGmax1
      do j = ixGmin2,ixGmax2
      e0(i,j) =  two + dexp(-8.d0 *dpi**two*t1*unit_time)*sin(two*dpi*x(i,j,&
         1))*sin(two*dpi*x(i,j,2))
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

  subroutine constant_var(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = 1.d2
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = 1.d0

    ! print*, global_time, dexp(-8.d0 *dpi**two*global_time*unit_time), dexp(dlog(1.d-2)/(8.d0*dpi**two)*8.d0*dpi**two)

  end subroutine constant_var

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
  double precision                   :: theoretical(ixImin1:ixImax1,&
     ixImin2:ixImax2)
  double precision                   :: residual(ixImin1:ixImax1,&
     ixImin2:ixImax2)
  ! double precision                   :: rad_flux(ixI^S,1:ndim), rad_pressure(ixI^S), fld_lambda(ixI^S), fld_R(ixI^S)

  ! call fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure)
  ! call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
  ! call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)

  theoretical(ixImin1:ixImax1,ixImin2:ixImax2) = spotpattern(x,ixImin1,ixImin2,&
     ixImax1,ixImax2,global_time)
  residual(ixImin1:ixImax1,ixImin2:ixImax2) = abs(theoretical(ixImin1:ixImax1,&
     ixImin2:ixImax2) - w(ixImin1:ixImax1,ixImin2:ixImax2,&
     r_e))/theoretical(ixImin1:ixImax1,ixImin2:ixImax2)

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = theoretical(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = residual(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)

  ! w(ixO^S,nw+1)=rad_flux(ixO^S,1)
  ! w(ixO^S,nw+2)=rad_flux(ixO^S,2)
  ! w(ixO^S,nw+3)=rad_pressure(ixO^S)
  ! w(ixO^S,nw+4)=fld_lambda(ixO^S)
  ! w(ixO^S,nw+5)=fld_R(ixO^S)


  !> Write error to file
  if (it == 0) open(1,file='3x3error_out1d-10')
  write(1,222) it, sum(residual(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2))/((ixOmax1-ixOmin1)*(ixOmax2-ixOmin2))
  if (it == it_max) close(1)
  222 format(i8,3e15.5E3)
  print*, it, sum(residual(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'theoretical residual'
  ! varnames = 'F1 F2 RP lam fld_R'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
