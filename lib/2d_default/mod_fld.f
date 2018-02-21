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
  subroutine fld_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x,energy,qsourcesplit,active)

    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_hd_phys, only: hd_get_pthermal  !needed to get temp

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision :: temperature(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: rad_flux(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: div_v(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision :: radiation_cooling(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: radiation_heating(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: photon_tiring(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: fld_boltzman_cgs = 5.67036713d-8 !*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density) !fix units !AREN'T THESE csts KNOWN ANYWHERE ELSE?
    integer :: idir

    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      !> Calculate the radiative flux using the FLD Approximation
      call fld_get_radflux(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, rad_flux, rad_pressure)

      do idir = 1,ndir
        !> Radiation force = kappa*rho/c *Flux
        radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,&
           idir) = fld_kappa*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_rho)*rad_flux(ixImin1:ixImax1,ixImin2:ixImax2,&
            idir) *unit_velocity/const_c

        !> Momentum equation source term
        w(ixImin1:ixImax1,ixImin2:ixImax2,iw_mom(idir)) = w(ixImin1:ixImax1,&
           ixImin2:ixImax2,iw_mom(idir)) + qdt * &
           radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,idir)
      end do

      if(energy .and. .not.block%e_is_internal) then

        !> calc Temperature as p/rho * (m_h*mu)/k
        call hd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,temperature)
        temperature(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(temperature(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2)/wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw_rho))*mp_cgs*fld_mu/kb_cgs*unit_length**&
           2/(unit_time**2 * unit_temperature)

        !> Cooling = -4pi kappa rho B = 4 kappa rho sigma T**4
        radiation_cooling(ixImin1:ixImax1,&
           ixImin2:ixImax2) = 4*fld_kappa*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_rho)*temperature(ixImin1:ixImax1,&
           ixImin2:ixImax2)*fld_boltzman_cgs*unit_time**3 * &
           unit_temperature**4 /(unit_length**3 *unit_density)
        !> Heating = c kappa rho E_rad
        radiation_heating(ixImin1:ixImax1,&
           ixImin2:ixImax2) = fld_kappa*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_rho)*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           r_e) *const_c/unit_velocity

        !> Energy equation source terms
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw_e) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw_e) + qdt * radiation_heating(ixImin1:ixImax1,&
           ixImin2:ixImax2) - qdt * radiation_cooling(ixImin1:ixImax1,&
           ixImin2:ixImax2)

      end if

      !> Photon tiring
      call divvector(wCT(ixImin1:ixImax1,ixImin2:ixImax2,iw_mom(:)),ixImin1,&
         ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,div_v)
      photon_tiring(ixImin1:ixImax1,ixImin2:ixImax2) = div_v(ixImin1:ixImax1,&
         ixImin2:ixImax2)/rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2)

      !> Radiation Energy source term
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw_e) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw_e) - qdt * photon_tiring(ixImin1:ixImax1,&
         ixImin2:ixImax2) - qdt * radiation_heating(ixImin1:ixImax1,&
         ixImin2:ixImax2) + qdt * radiation_cooling(ixImin1:ixImax1,&
         ixImin2:ixImax2)

    end if

  end subroutine fld_add_source

  !> Calculate Radiation Flux
  !> Returns Radiation flux and radiation pressure
  subroutine fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, rad_flux, rad_pressure)
    use mod_global_parameters
    !use geometry, only: gradient

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: rad_flux(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim), rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: fld_lambda(ixImin1:ixImax1,ixImin2:ixImax2),&
        fld_R(ixImin1:ixImax1,ixImin2:ixImax2), normgrad2(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    integer :: idir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixImin1:ixImax1,ixImin2:ixImax2) = zero
    do idir = 1,ndir
      call gradient(w(ixImin1:ixImax1,ixImin2:ixImax2, iw_r_e),ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idir,&
         grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2,idir)) !!! IS IT WRONG TO USE ixOmin1?????,ixOmin2?????,ixOmax1?????,ixOmax2?????
      normgrad2(ixImin1:ixImax1,ixImin2:ixImax2) = normgrad2(ixImin1:ixImax1,&
         ixImin2:ixImax2) + grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2,idir)**2
    end do
    fld_R(ixImin1:ixImax1,ixImin2:ixImax2) = dsqrt(normgrad2(ixImin1:ixImax1,&
       ixImin2:ixImax2))/(fld_kappa*w(ixImin1:ixImax1,ixImin2:ixImax2,&
       iw_rho)*w(ixImin1:ixImax1,ixImin2:ixImax2,r_e))

    !> Calculate the flux limiter, lambda
    !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
    fld_lambda(ixImin1:ixImax1,ixImin2:ixImax2) = (2+fld_R(ixImin1:ixImax1,&
       ixImin2:ixImax2))/(6+3*fld_R(ixImin1:ixImax1,&
       ixImin2:ixImax2)+fld_R(ixImin1:ixImax1,ixImin2:ixImax2)**2)
    ! use two instead of 2??? 3 three 6 six

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndir
      rad_flux(ixImin1:ixImax1,ixImin2:ixImax2,&
          idir) = -const_c/unit_velocity*fld_lambda(ixImin1:ixImax1,&
         ixImin2:ixImax2)/(fld_kappa*w(ixImin1:ixImax1,ixImin2:ixImax2,&
         iw_rho)) *grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2,idir)
    end do

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2) = &
       (fld_lambda(ixImin1:ixImax1,ixImin2:ixImax2) + &
       fld_lambda(ixImin1:ixImax1,ixImin2:ixImax2)**2 * fld_R(ixImin1:ixImax1,&
       ixImin2:ixImax2)**2) * w(ixImin1:ixImax1,ixImin2:ixImax2, iw_r_e)

  end subroutine fld_get_radflux


  subroutine fld_get_flux_cons(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    double precision                :: rad_flux(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:ndim), rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir, itr

    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, rad_flux, rad_pressure)
    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        iw_mom(idim)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, iw_rho)

    !> Radiation energy is v_i*r_e   (??+ F_i??)
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_e) = v * w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, r_e) + rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)

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


  subroutine fld_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_hd_phys, only: hd_get_pthermal  !needed to get temp

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: dx1,dx2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(out) :: dtnew
    double precision :: dxinv(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

    double precision :: temperature(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: rad_flux(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: div_v(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision :: radiation_cooling(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: radiation_heating(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: photon_tiring(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: fld_boltzman_cgs = 5.67036713d-8 !*unit_time**3 * unit_temperature**4 /(unit_length**3 *unit_density) !fix units !AREN'T THESE csts KNOWN ANYWHERE ELSE?
    integer :: idir

    !> Calculate the radiative flux using the FLD Approximation
    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, rad_flux, rad_pressure)

    do idir = 1,ndir
      !> Radiation force = kappa*rho/c *Flux
      radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,&
         idir) = fld_kappa*w(ixImin1:ixImax1,ixImin2:ixImax2,&
         iw_rho)*rad_flux(ixImin1:ixImax1,ixImin2:ixImax2,&
          idir) *unit_velocity/const_c
    end do

    !> calc Temperature as p/rho * (m_h*mu)/k
    call hd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,temperature)
    temperature(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(temperature(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       iw_rho))*mp_cgs*fld_mu/kb_cgs*unit_length**2/(unit_time**2 * &
       unit_temperature)

    !> Cooling = -4pi kappa rho B = 4 kappa rho sigma T**4
    radiation_cooling(ixImin1:ixImax1,ixImin2:ixImax2) = &
       4*fld_kappa*w(ixImin1:ixImax1,ixImin2:ixImax2,&
       iw_rho)*temperature(ixImin1:ixImax1,&
       ixImin2:ixImax2)*fld_boltzman_cgs*unit_time**3 * unit_temperature**4 &
       /(unit_length**3 *unit_density)
    !> Heating = c kappa rho E_rad
    radiation_heating(ixImin1:ixImax1,ixImin2:ixImax2) = &
       fld_kappa*w(ixImin1:ixImax1,ixImin2:ixImax2,iw_rho)*w(ixImin1:ixImax1,&
       ixImin2:ixImax2,r_e) *const_c/unit_velocity

    !> Photon tiring
    call divvector(w(ixImin1:ixImax1,ixImin2:ixImax2,iw_mom(:)),ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,div_v)
    photon_tiring(ixImin1:ixImax1,ixImin2:ixImax2) = div_v(ixImin1:ixImax1,&
       ixImin2:ixImax2)/rad_pressure(ixImin1:ixImax1,ixImin2:ixImax2)

    !> New dt based on radiation_force
    do idir = 1,ndir
      if (minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
         iw_rho)/radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,&
         idir)*dxinv(ixImin1:ixImax1,ixImin2:ixImax2,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_rho)/radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,&
           idir)*dxinv(ixImin1:ixImax1,ixImin2:ixImax2,idir)))) !chuck in the density somwhere?
        print*, 'radiation_force dt', &
        minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_rho)/radiation_force(ixImin1:ixImax1,ixImin2:ixImax2,&
           idir)*dxinv(ixImin1:ixImax1,ixImin2:ixImax2,idir)))
      end if
    end do

    !> New dt based on cooling term
    do idir = 1,ndir
      if (minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
         iw_mom(idir))/radiation_cooling*dxinv(:,:,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_mom(idir))/radiation_cooling*dxinv(:,:,idir)))) !chuck in the density somwhere?
        print*, 'cooling term dt',  minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,iw_mom(idir))/radiation_cooling*dxinv(:,:,idir)))
      end if
    end do

    !> New dt based on heating term
    do idir = 1,ndir
      if (minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
         iw_mom(idir))/radiation_heating*dxinv(:,:,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_mom(idir))/radiation_heating*dxinv(:,:,idir)))) !chuck in the density somwhere?
        print*, 'heating term dt', minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,iw_mom(idir))/radiation_heating*dxinv(:,:,idir)))
      end if
    end do

    !> New dt based on photon tiring
    do idir = 1,ndir
      if (minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
         iw_mom(idir))/photon_tiring*dxinv(:,:,idir))) > zero) then
        dtnew = min( dtnew, minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,&
           iw_mom(idir))/photon_tiring*dxinv(:,:,idir)))) !chuck in the density somwhere?
        print*, 'photon tiring dt', minval(sqrt(w(ixImin1:ixImax1,ixImin2:ixImax2,iw_mom(idir))/photon_tiring*dxinv(:,:,idir)))
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
