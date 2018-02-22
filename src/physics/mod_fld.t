!> Module for including flux limited diffusion in hydrodynamics simulations
module mod_fld
    implicit none

    !> source split or not
    logical :: fld_split = .false.

    !> Opacity per unit of unit_density
    double precision, public :: fld_kappa = 1.d0

    !> mean particle mass
    double precision, public :: fld_mu = 1.d0

    !> NICOLAS MOENS
    !> Index of the radiation energy density
    integer, protected :: iw_r_e = -1

    !> Index of the radiation energy
    integer, public, protected :: r_e

    !> public methods
    public :: fld_add_source
    public :: fld_get_radflux
    !public :: fld_get_flux
    public :: fld_get_flux_cons
    public :: fld_get_dt

  contains

  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa, fld_mu, fld_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do
  end subroutine fld_params_read

  !> Set radiation energy variable
  function var_set_radiation_energy() result(iw)
    use mod_variables
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_r_e              = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'r_e'
    prim_wnames(nwflux) = 'r_e'
  end function var_set_radiation_energy

  subroutine fld_init()

    use mod_global_parameters
    use mod_variables

    !> read par files
    call fld_params_read(par_files)

    !> Make kappa dimensionless
    fld_kappa = fld_kappa*unit_length*unit_density

    !> set radiation energy mod_variables
    r_e = var_set_radiation_energy()

  end subroutine fld_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine fld_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)

    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_hd_phys, only: hd_get_pthermal  !needed to get temp
    !use mod_phys, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision :: temperature(ixI^S)
    double precision :: rad_flux(ixI^S,1:ndim)
    double precision :: rad_pressure(ixI^S)
    double precision :: div_v(ixI^S)

    double precision :: radiation_force(ixI^S,1:ndim)
    double precision :: radiation_cooling(ixI^S)
    double precision :: radiation_heating(ixI^S)
    double precision :: photon_tiring(ixI^S)

    double precision :: fld_boltzman_cgs = 5.67036713d-8 !*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density) !fix units !AREN'T THESE csts KNOWN ANYWHERE ELSE?
    integer :: idir

    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      !> Calculate the radiative flux using the FLD Approximation
      call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux, rad_pressure)

      do idir = 1,ndir
        !> Radiation force = kappa*rho/c *Flux
        radiation_force(ixI^S,idir) = fld_kappa*wCT(ixI^S,iw_rho)*rad_flux(ixI^S, idir) *unit_velocity/const_c

        !> Momentum equation source term
        w(ixI^S,iw_mom(idir)) = w(ixI^S,iw_mom(idir)) &
            + qdt * radiation_force(ixI^S,idir)
      end do

      if(energy .and. .not.block%e_is_internal) then

        !> calc Temperature as p/rho * (m_h*mu)/k
        call hd_get_pthermal(w,x,ixI^L,ixO^L,temperature)
        temperature(ixO^S)=(temperature(ixO^S)/wCT(ixO^S,iw_rho))*mp_cgs*fld_mu/kb_cgs&
        *unit_length**2/(unit_time**2 * unit_temperature)

        !> Cooling = -4pi kappa rho B = 4 kappa rho sigma T**4
        radiation_cooling(ixI^S) = 4*fld_kappa*wCT(ixI^S,iw_rho)*temperature(ixI^S)&
        *fld_boltzman_cgs*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density)
        !> Heating = c kappa rho E_rad
        radiation_heating(ixI^S) = fld_kappa*wCT(ixI^S,iw_rho)*wCT(ixI^S,r_e) *const_c/unit_velocity

        !> Energy equation source terms
        w(ixO^S,iw_e) = w(ixO^S,iw_e) &
           + qdt * radiation_heating(ixI^S) &
           - qdt * radiation_cooling(ixI^S)

      end if

      !> Photon tiring
      call divvector(wCT(ixI^S,iw_mom(:)),ixI^L,ixO^L,div_v)
      photon_tiring(ixI^S) = div_v(ixI^S)/rad_pressure(ixI^S)

      !> Radiation Energy source term
      w(ixO^S,iw_e) = w(ixO^S,iw_e) &
         - qdt * photon_tiring(ixI^S) &
         - qdt * radiation_heating(ixI^S) &
         + qdt * radiation_cooling(ixI^S)

    end if

  end subroutine fld_add_source

  !> Calculate Radiation Flux
  !> Returns Radiation flux and radiation pressure
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, rad_pressure)
    use mod_global_parameters
    !use geometry, only: gradient

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_flux(ixI^S, 1:ndim), rad_pressure(ixI^S)
    double precision :: fld_lambda(ixI^S), fld_R(ixI^S), normgrad2(ixI^S)
    double precision :: grad_r_e(ixI^S, 1:ndim)
    integer :: idir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixI^S) = zero
    do idir = 1,ndir
      !call gradient(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir)) !!! IS IT WRONG TO USE ixO^L?????
      call grad(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e(ixI^S,idir))
      normgrad2(ixI^S) = normgrad2(ixI^S) + grad_r_e(ixI^S,idir)**2
    end do
    fld_R(ixI^S) = dsqrt(normgrad2(ixI^S))/(fld_kappa*w(ixI^S,iw_rho)*w(ixI^S,r_e))

    !> Calculate the flux limiter, lambda
    !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
    fld_lambda(ixI^S) = (2+fld_R(ixI^S))/(6+3*fld_R(ixI^S)+fld_R(ixI^S)**2)
    ! use two instead of 2??? 3 three 6 six

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndir
      rad_flux(ixI^S, idir) = -const_c/unit_velocity*fld_lambda(ixI^S)/(fld_kappa*w(ixI^S,iw_rho)) *grad_r_e(ixI^S,idir)
    end do

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    rad_pressure(ixI^S) = (fld_lambda(ixI^S) + fld_lambda(ixI^S)**2 * fld_R(ixI^S)**2) * w(ixI^S, iw_r_e)

  end subroutine fld_get_radflux

  subroutine grad(q,ixI^L,ix^L,idir,x,gradq)
    ! Compute the true gradient of a scalar q within ixL in direction idir ie :
    !  - in cylindrical : d_/dr , (1/r)*d_/dth , d_/dz
    !  - in spherical   : d_/dr , (1/r)*d_/dth , (1/rsinth)*d_/dphi
    use mod_global_parameters

    integer :: ixI^L, ix^L, idir
    double precision :: q(ixI^S), gradq(ixI^S)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    integer :: jx^L, hx^L
    !-----------------------------------------------------------------------------
    jx^L=ix^L+kr(idir,^D);
    hx^L=ix^L-kr(idir,^D);
    gradq(ix^S)=(q(jx^S)-q(hx^S))/(x(jx^S,idir)-x(hx^S,idir))
    select case (typeaxial)
    case('slab') ! nothing to do
    case('cylindrical')
      if (idir==phi_) gradq(ix^S)=gradq(ix^S)/ x(ix^S,r_)
    case('spherical')
      if (idir==2   ) gradq(ix^S)=gradq(ix^S)/ x(ix^S,r_)
      if (idir==phi_) gradq(ix^S)=gradq(ix^S)/(x(ix^S,r_)*dsin(x(ix^S,2)))
    case default
      call mpistop('Unknown geometry')
    end select
  end subroutine grad

  subroutine fld_get_flux_cons(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: rad_flux(ixI^S, 1:ndim), rad_pressure(ixI^S)
    double precision                :: v(ixI^S)
    integer                         :: idir, itr

    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, rad_pressure)
    v(ixO^S) = w(ixO^S, iw_mom(idim)) / w(ixO^S, iw_rho)

    !> Radiation energy is v_i*r_e   (??+ F_i??)
    f(ixO^S, r_e) = v * w(ixO^S, r_e) + rad_flux(ixO^S,idim)

  end subroutine fld_get_flux_cons


  ! ! Calculate flux f_idim[iw]
  ! subroutine fld_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
  !   use mod_global_parameters
  !
  !   integer, intent(in)             :: ixI^L, ixO^L, idim
  !   ! conservative w
  !   double precision, intent(in)    :: wC(ixI^S, 1:nw)
  !   ! primitive w
  !   double precision, intent(in)    :: w(ixI^S, 1:nw)
  !   double precision, intent(in)    :: x(ixI^S, 1:ndim)
  !   double precision, intent(out)   :: f(ixI^S, nwflux)
  !   integer                         :: idir
  !
  !   double precision :: rad_flux(ixI^S, 1:ndim), rad_pressure(ixI^S)
  !   call fld_get_radflux(wC, x, ixI^L, ixO^L, rad_flux, rad_pressure)
  !
  !   !> Radiation energy is v_i*r_e   (??+ F_i??)
  !   f(ixO^S, r_e) = w(ixO^S,iw_mom(idim)) * wC(ixO^S, r_e) + rad_flux(ixO^S,idim)*w(ixO^S,iw_rho) !Do I need rho, if so w or wC
  !   !print*, f
  !
  ! end subroutine fld_get_flux


  subroutine fld_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_hd_phys, only: hd_get_pthermal  !needed to get temp
    !use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(out) :: dtnew
    double precision :: dxinv(ixI^S,1:ndim)

    double precision :: temperature(ixI^S)
    double precision :: rad_flux(ixI^S,1:ndim)
    double precision :: rad_pressure(ixI^S)
    double precision :: div_v(ixI^S)

    double precision :: radiation_force(ixI^S,1:ndim)
    double precision :: radiation_cooling(ixI^S)
    double precision :: radiation_heating(ixI^S)
    double precision :: photon_tiring(ixI^S)

    double precision :: fld_boltzman_cgs = 5.67036713d-8 !*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density) !fix units !AREN'T THESE csts KNOWN ANYWHERE ELSE?
    integer :: idir

    !> Calculate the radiative flux using the FLD Approximation
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, rad_pressure)

    do idir = 1,ndir
      !> Radiation force = kappa*rho/c *Flux
      radiation_force(ixI^S,idir) = fld_kappa*w(ixI^S,iw_rho)*rad_flux(ixI^S, idir) *unit_velocity/const_c
    end do

    !> calc Temperature as p/rho * (m_h*mu)/k
    call hd_get_pthermal(w,x,ixI^L,ixO^L,temperature)
    temperature(ixO^S)=(temperature(ixO^S)/w(ixO^S,iw_rho))*mp_cgs*fld_mu/kb_cgs&
    *unit_length**2/(unit_time**2 * unit_temperature)

    !> Cooling = -4pi kappa rho B = 4 kappa rho sigma T**4
    radiation_cooling(ixI^S) = 4*fld_kappa*w(ixI^S,iw_rho)*temperature(ixI^S)&
    *fld_boltzman_cgs*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density)
    !> Heating = c kappa rho E_rad
    radiation_heating(ixI^S) = fld_kappa*w(ixI^S,iw_rho)*w(ixI^S,r_e) *const_c/unit_velocity

    !> Photon tiring
    call divvector(w(ixI^S,iw_mom(:)),ixI^L,ixO^L,div_v)
    photon_tiring(ixI^S) = div_v(ixI^S)/rad_pressure(ixI^S)

    !> New dt based on radiation_force
    do idir = 1,ndir
      if (minval(sqrt(w(ixI^S,iw_rho)/radiation_force(ixI^S,idir)*dxinv(ixI^S,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixI^S,iw_rho)/radiation_force(ixI^S,idir)*dxinv(ixI^S,idir)))) !chuck in the density somwhere?
        print*, 'radiation_force dt', &
        minval(sqrt(w(ixI^S,iw_rho)/&
        radiation_force(ixI^S,idir)*&
        dxinv(ixI^S,idir)))
      end if
    end do

    !> New dt based on cooling term
    do idir = 1,ndir
      if (minval(sqrt(w(ixI^S,iw_mom(idir))/radiation_cooling*dxinv(:,:,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixI^S,iw_mom(idir))/radiation_cooling*dxinv(:,:,idir)))) !chuck in the density somwhere?
        print*, 'cooling term dt',  minval(sqrt(w(ixI^S,iw_mom(idir))/radiation_cooling*dxinv(:,:,idir)))
      end if
    end do

    !> New dt based on heating term
    do idir = 1,ndir
      if (minval(sqrt(w(ixI^S,iw_mom(idir))/radiation_heating*dxinv(:,:,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixI^S,iw_mom(idir))/radiation_heating*dxinv(:,:,idir)))) !chuck in the density somwhere?
        print*, 'heating term dt', minval(sqrt(w(ixI^S,iw_mom(idir))/radiation_heating*dxinv(:,:,idir)))
      end if
    end do

    !> New dt based on photon tiring
    do idir = 1,ndir
      if (minval(sqrt(w(ixI^S,iw_mom(idir))/photon_tiring*dxinv(:,:,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixI^S,iw_mom(idir))/photon_tiring*dxinv(:,:,idir)))) !chuck in the density somwhere?
        print*, 'photon tiring dt', minval(sqrt(w(ixI^S,iw_mom(idir))/photon_tiring*dxinv(:,:,idir)))
      end if
    end do

  end subroutine fld_get_dt


  ! subroutine alternating_direction(w,ixI^L,ixO^L,dtnew,dx^D,x)
  !   use mod_global_parameters
  !
  !   integer, intent(in) :: ixI^L, ixO^L
  !   double precision, intent(in) :: dx^D, x(ixI^S,1:ndim), dtnew
  !   double precision, intent(inout) :: w(ixI^S,1:nw)
  !
  !   double precision :: radiation_energy(ixI^S) = w(ixI^S,r_e)
  !   double precision :: matrix
  !   !SOlVE IN ONE DIRECTION
  !   ! call in routine to set matrix, beta1, B1 and D1
  !   ! call in routine solve system using some kind of decomposition,
  !   ! output should be placed in radiation_energy
  !
  !   !SOLVE IN OTHER DIRECTION
  !   ! call in routine to set matrix, beta2, B2 and D2
  !   ! call in routine to solve system
  !   ! output should be placed in radiation_energy
  !
  ! end subroutine alternating_direction

end module mod_fld
