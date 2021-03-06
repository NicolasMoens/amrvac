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

  !Fix dimensionless stuff here
  unit_length        = 1.d0                                         ! cm
  unit_numberdensity = 1.d3 !cm-3,cm-3
  unit_temperature   = 1.d0                                         ! K

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

    double precision :: temperature(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
    double precision :: int_energy(ixGmin1:ixGmax1,ixGmin2:ixGmax2)

    integer :: i

    ! Set initial values for w
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = d_inflo
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(1)) = 0.d0
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(2)) = 0.d0

    !> Try some ~1/r^2 init condition for energy
    int_energy = T_eff * one/(hd_gamma-one) * w(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2, rho_)*kb_cgs/(mp_cgs*fld_mu)*unit_time**2/(unit_length**&
       2)
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = int_energy
    w(ixmin1:ixmax1,ixmin2:ixmax2, e_) = w(ixmin1:ixmax1,ixmin2:ixmax2,&
        e_) + 1.d-2/x(ixmin1:ixmax1,ixmin2:ixmax2,2)**2

    ! !> Radiative Equilibrium, heating = cooling
    ! call hd_get_pthermal(w,x,ixG^L,ixG^L,temperature)
    !
    ! temperature(ixG^S) = temperature(ixG^S)/w(ixG^S,iw_rho)*mp_cgs*fld_mu/kb_cgs&
    ! *unit_length**2/(unit_time**2 * unit_temperature)
    !
    ! w(ixG^S,r_e) = 4*unit_velocity/const_c* temperature(ixG^S)**4*&
    ! fld_boltzman_cgs*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density)

    !> Radiative Equilibrium, heating = cooling
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = 4*unit_velocity/const_c* &
       (T_eff/unit_temperature)**4*fld_boltzman_cgs*unit_time**3 * &
       unit_temperature**4 /(unit_length**3 *unit_density)

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
      w(:,ixBmax2, rho_) =  d_inflo !w(:,ixBmax2+1, rho_)
      w(:,ixBmax2, mom(1)) = w(:,ixBmax2+1, mom(1))

      where (w(:,ixBmax2+1 , mom(2)) > zero)
        w(:,ixBmax2, mom(2)) = w(:,ixBmax2+1 , mom(2))
      elsewhere
        w(:,ixBmax2, mom(2)) = zero
      end where

      e_inflo = T_eff * one/(hd_gamma-one) * d_inflo &
         *kb_cgs/(mp_cgs*fld_mu)*unit_time**2/(unit_length**2)
      w(:,ixBmax2, e_) = e_inflo !1.d0 !w(:,ixBmax2+1, e_)
      w(:,ixBmax2, r_e) = 4*unit_velocity/const_c* &
         (T_eff/unit_temperature)**4*fld_boltzman_cgs*unit_time**3 * &
         unit_temperature**4 /(unit_length**3 *unit_density) !w(:,ixBmax2+1, r_e)
      w(:,ixBmin2,:) = w(:,ixBmax2,:)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

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

  double precision :: mass, distance, cavendish

  mass = 1d35*unit_density*unit_length**3
  distance = 1d3*unit_length
  cavendish = 6.67d-8**unit_density*unit_time

  gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
  gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
     2) = -cavendish*mass/((distance + x(ixImin1:ixImax1,ixImin2:ixImax2,&
     2))**2)

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

  varnames = 'F1 F2 RP L R'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
