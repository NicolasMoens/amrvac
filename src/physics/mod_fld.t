!> Module for including flux limited diffusion in hydrodynamics simulations
module mod_fld
    implicit none

    !> source split or not
    logical :: fld_split = .false.

    !> Opacity per unit of unit_density
    double precision, public :: fld_kappa = 0.4d0

    !> mean particle mass
    double precision, public :: fld_mu = 0.6d0

    !> Index of the radiation energy density
    integer, protected :: iw_r_e = -1

    !> Index of the radiation energy
    integer, public, protected :: r_e

    !> Dimensionless Boltzman constante sigma
    double precision, public :: fld_sigma_0

    !> Dimensionless speed of light
    double precision, public :: fld_speedofligt_0

    !> public methods
    public :: fld_add_source
    public :: fld_get_radflux
    public :: fld_get_fluxlimiter
    public :: fld_get_flux
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
  !> In  Erg/cm^3
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

    !> Make kappa dimensionless !!!STILL NEED TO MULTIPLY W RHO
    fld_kappa = fld_kappa*unit_time*unit_velocity

    !> set radiation energy mod_variables
    r_e = var_set_radiation_energy()

    !> Dimensionless speed of light
    fld_speedofligt_0 = const_c/unit_velocity

    !> Dimensionless Boltzman constante sigma
    fld_sigma_0 = (5.67051d-5*unit_temperature**4)/(unit_velocity*unit_pressure)

  end subroutine fld_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine fld_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)

    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_physics, only: phys_get_pthermal  !needed to get temp

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

    integer :: idir, i

    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      !> Calculate the radiative flux using the FLD Approximation
      call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux, rad_pressure)

      do idir = 1,ndir
        !> Radiation force = kappa*rho/c *Flux
        radiation_force(ixI^S,idir) = fld_kappa*wCT(ixI^S,iw_rho)/fld_speedofligt_0*rad_flux(ixI^S, idir)

        !> Momentum equation source term
        w(ixI^S,iw_mom(idir)) = w(ixI^S,iw_mom(idir)) &
            + qdt * radiation_force(ixI^S,idir)
      end do

      if(energy .and. .not.block%e_is_internal) then

        !> Get pressure
        call phys_get_pthermal(wCT,x,ixI^L,ixO^L,temperature)

        !> calc Temperature as p/rho
        temperature(ixI^S)=(temperature(ixO^S)/wCT(ixO^S,iw_rho))

        !> Cooling = 4 pi kappa B = 4 kappa sigma T**4
        radiation_cooling(ixO^S) = 4*fld_kappa*wCT(ixO^S,iw_rho)*fld_sigma_0*temperature(ixO^S)**4
        !> Heating = c kappa E_rad
        radiation_heating(ixO^S) = fld_speedofligt_0*fld_kappa*wCT(ixO^S,iw_rho)*wCT(ixO^S,iw_r_e)

        !> Write energy to file
        if (it == 0) open(1,file='energy_out6')
        write(1,222) it,global_time,w(5,5,iw_e)
        if (it == it_max) close(1)
        222 format(i8,2e15.5E3)

        ! print*, 'dt', dt

        !> Energy equation source terms
        w(ixO^S,iw_e) = w(ixO^S,iw_e) &
           + qdt * radiation_heating(ixO^S) &
           - qdt * radiation_cooling(ixO^S)

        ! print*, 'fld_kappa*wCT(ixO^S,iw_rho)',fld_kappa*wCT(5,5,iw_rho)
        !
        ! print*, 'fld_speedofligt_0*fld_kappa*wCT(ixO^S,iw_rho)', fld_speedofligt_0*fld_kappa*wCT(5,5,iw_rho)
        ! print*, 'qdt * radiation_heating(ixO^S)', qdt * radiation_heating(5,5)
        !
        ! print*, '4*fld_kappa*wCT(ixO^S,iw_rho)*fld_sigma_0', 4*fld_kappa*wCT(5,5,iw_rho)*fld_sigma_0
        ! print*, 'qdt * radiation_cooling(ixO^S)', qdt * radiation_cooling(5,5)

        print*,  it,global_time,w(5,5,iw_e)

      end if

      !> Photon tiring
      call divvector(wCT(ixI^S,iw_mom(:)),ixI^L,ixO^L,div_v)
      photon_tiring(ixO^S) = div_v(ixO^S)/rad_pressure(ixO^S)

      !> Radiation Energy source term
      w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) &
         - qdt * photon_tiring(ixO^S) &
         - qdt * radiation_heating(ixO^S) &
         + qdt * radiation_cooling(ixO^S)

    end if

  end subroutine fld_add_source


  !> Calculate fld flux limiter
  subroutine fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
    use mod_global_parameters
    !use geometry, only: gradient

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: fld_R(ixI^S), fld_lambda(ixI^S)
    double precision ::  normgrad2(ixI^S)
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

  end subroutine fld_get_fluxlimiter


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

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndir
      rad_flux(ixI^S, idir) = -const_c/unit_velocity*fld_lambda(ixI^S)/(fld_kappa*w(ixI^S,iw_rho)) *grad_r_e(ixI^S,idir)
    end do

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    rad_pressure(ixI^S) = (fld_lambda(ixI^S) + fld_lambda(ixI^S)**2 * fld_R(ixI^S)**2) * w(ixI^S, iw_r_e)

  end subroutine fld_get_radflux


  ! Calculate flux f_idim[iw]
  subroutine fld_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    ! primitive w
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    integer                         :: idir

    double precision :: rad_flux(ixI^S, 1:ndim), rad_pressure(ixI^S)
    call fld_get_radflux(wC, x, ixI^L, ixO^L, rad_flux, rad_pressure)

    !> Radiation energy is v_i*r_e   (??+ F_i??)
    f(ixO^S, r_e) = w(ixO^S,iw_mom(idim)) * wC(ixO^S, r_e) + rad_flux(ixO^S,idim)*w(ixO^S,iw_rho) !Do I need rho, if so w or wC
    !print*, f

  end subroutine fld_get_flux


  subroutine fld_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_physics, only: phys_get_pthermal  !needed to get temp

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

    integer :: idir

    !> Calculate the radiative flux using the FLD Approximation
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, rad_pressure)

    do idir = 1,ndir
      !> Radiation force = kappa*rho/c *Flux
      radiation_force(ixI^S,idir) = fld_kappa*w(ixI^S,iw_rho)/fld_speedofligt_0*rad_flux(ixI^S, idir)
    end do

    !> Get pressure
    call phys_get_pthermal(w,x,ixI^L,ixO^L,temperature)
    !> calc Temperature as p/rho
    temperature(ixI^S)=(temperature(ixI^S)/w(ixI^S,iw_rho))

    !> Cooling = 4 pi kappa B = 4 kappa sigma T**4
    radiation_cooling(ixI^S) = 4*fld_kappa*w(ixI^S,iw_rho)*fld_sigma_0*temperature(ixI^S)**4
    !> Heating = c kappa E_rad
    radiation_heating(ixI^S) = fld_speedofligt_0*fld_kappa*w(ixI^S,iw_rho)*w(ixI^S,iw_r_e)

    !> Photon tiring
    call divvector(w(ixI^S,iw_mom(:)),ixI^L,ixO^L,div_v)
    photon_tiring(ixI^S) = div_v(ixI^S)/rad_pressure(ixI^S)

    !> New dt based on radiation_force
    do idir = 1,ndir
      if (minval(dsqrt(w(ixI^S,iw_rho)/radiation_force(ixI^S,idir)*dxinv(ixI^S,idir))) > zero) then
        dtnew = min( dtnew, minval(dsqrt(w(ixI^S,iw_rho)/radiation_force(ixI^S,idir)*dxinv(ixI^S,idir)))) !chuck in the density somwhere?
        print*, 'radiation_force dt', &
        minval(dsqrt(w(ixI^S,iw_rho)/&
        radiation_force(ixI^S,idir)*&
        dxinv(ixI^S,idir)))
      end if
    end do

    ! !> New dt based on cooling term
    ! do idir = 1,ndir
    !   if (minval(dsqrt(w(ixI^S,iw_mom(idir))/radiation_cooling*dxinv(:,:,idir))) > zero) then
    !     dtnew = min( dtnew, minval(dsqrt(w(ixI^S,iw_mom(idir))/radiation_cooling*dxinv(:,:,idir)))) !chuck in the density somwhere?
    !     print*, 'cooling term dt',  minval(ddsqrt(w(ixI^S,iw_mom(idir))/radiation_cooling*dxinv(:,:,idir)))
    !   end if
    ! end do
    !
    ! !> New dt based on heating term
    ! do idir = 1,ndir
    !   if (minval(dsqrt(w(ixI^S,iw_mom(idir))/radiation_heating*dxinv(:,:,idir))) > zero) then
    !     dtnew = min( dtnew, minval(dsqrt(w(ixI^S,iw_mom(idir))/radiation_heating*dxinv(:,:,idir)))) !chuck in the density somwhere?
    !     print*, 'heating term dt', minval(dsqrt(w(ixI^S,iw_mom(idir))/radiation_heating*dxinv(:,:,idir)))
    !   end if
    ! end do
    !
    ! !> New dt based on photon tiring
    ! do idir = 1,ndir
    !   if (minval(dsqrt(w(ixI^S,iw_mom(idir))/photon_tiring*dxinv(:,:,idir))) > zero) then
    !     dtnew = min( dtnew, minval(dsqrt(w(ixI^S,iw_mom(idir))/photon_tiring*dxinv(:,:,idir)))) !chuck in the density somwhere?
    !     print*, 'photon tiring dt', minval(dsqrt(w(ixI^S,iw_mom(idir))/photon_tiring*dxinv(:,:,idir)))
    !   end if
    ! end do

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
  !
  !    !SET PSEUDO-TIMESTEPS
  !
  !   !SOlVE IN ONE DIRECTION
  !   ! call in routine to set matrix, beta1, B1 and D1
  !   call set_system_matrix
  !   ! call in routine solve system using some kind of decomposition,
  !   ! output should be placed in radiation_energy
  !
  !   !SOLVE IN OTHER DIRECTION
  !   ! call in routine to set matrix, beta2, B2 and D2
  !   ! call in routine to solve system
  !   ! output should be placed in radiation_energy
  !
  ! end subroutine alternating_direction


  ! subroutine set_system_matrix(w, x, ixI^L, ixO^L, idir, dt, ps_dt matrix, b_vec)
  !   use mod_global_parameters
  !
  !   integer, intent(in)          :: ixI^L, ixO^L, idir
  !   double precision, intent(in) :: dt, ps_dt
  !   double precision, intent(in) :: w(ixI^S, nw)
  !   double precision, intent(in) :: x(ixI^S, 1:ndim)
  !   double precision, intent(out):: matrix
  !   double precision :: fld_lambda(ixI^S), fld_R(ixI^S), normgrad2(ixI^S)
  !   double precision :: grad_r_e(ixI^S)
  !   double precision :: Diff_coef(ixI^S)
  !   double precision :: sys_h(ixI^S), sys_beta(ixI^S)
  !   integer :: jx^L
  !   integer :: idir1, i
  !
  !   !> Calculate R everywhere
  !   !> |grad E|/(rho kappa E)
  !   normgrad2(ixI^S) = zero
  !   do idir1 = 1,ndir
  !     !call gradient(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir)) !!! IS IT WRONG TO USE ixO^L?????
  !     call grad_face(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir1,idir,x,grad_r_e(ixI^S,idir))
  !     normgrad2(ixI^S) = normgrad2(ixI^S) + grad_r_e(ixI^S,idir)**2
  !   end do
  !   fld_R(ixI^S) = dsqrt(normgrad2(ixI^S))/(fld_kappa*w(ixI^S,iw_rho)*w(ixI^S,r_e))
  !
  !   !> Calculate the flux limiter, lambda
  !   !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
  !   fld_lambda(ixI^S) = (2+fld_R(ixI^S))/(6+3*fld_R(ixI^S)+fld_R(ixI^S)**2)
  !
  !   !> Diffusion coeficient c*lamda/ (rho*kappa)
  !   Diff_coef(ixI^S) = const_c/unit_velocity*fld_lambda(ixI^S)/(fld_kappa*w(ixI^S,iw_rho))
  !
  !   !> Calculate h
  !   !> This should be a single float !?!?!?!
  !   sys_h(ixI^S) = ps_dt/(2*dx(ixI^S,idir)**2)
  !
  !   !> Calculate beta
  !   jxI^S = ixI^S+kr(idir,^D)
  !   sys_beta(ixI^S) = one + ps_dt/(2*dt) + sys_h(ixI^S)*(Diff_coef(jxI^S) + Diff_coef(ixI^S))
  !
  !   !> Assemble matrix
  !   !> What to do with boundary conditions?!?!?
  !   do i =
  !     matrix(ixI^S,i,i) = sys_beta(ixI^S)
  !     matrix(ixI^S,i+1,i) = -sys_h*Diff_coef
  !     matrix(ixI^S,i,i+1) = -sys_h*Diff_coef
  !   end do
  !
  ! end subroutine set_system_matrix



  ! subroutine set_system_b_vec(w, x, ixI^L, ixO^L, idir, dt, ps_dt matrix, b_vec)
  !   use mod_global_parameters
  !
  !   integer, intent(in)          :: ixI^L, ixO^L, idir
  !   double precision, intent(in) :: dt, ps_dt
  !   double precision, intent(in) :: w(ixI^S, nw), Rad_e_n(ixI^S)
  !   double precision, intent(in) :: x(ixI^S, 1:ndim)
  !   double precision, intent(out):: b_vec
  !   double precision :: fld_lambda(ixI^S), fld_R(ixI^S), normgrad2(ixI^S)
  !   double precision :: grad_r_e(ixI^S)
  !   double precision :: Diff_coef(ixI^S)
  !   double precision :: sys_h(ixI^S)
  !   integer :: jx^L
  !   integer :: idir1, i
  !
  !   !> Calculate R everywhere
  !   !> |grad E|/(rho kappa E)
  !   normgrad2(ixI^S) = zero
  !   do idir1 = 1,ndir
  !     !call gradient(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir)) !!! IS IT WRONG TO USE ixO^L?????
  !     call grad_face(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir1,idir,x,grad_r_e(ixI^S,idir))
  !     normgrad2(ixI^S) = normgrad2(ixI^S) + grad_r_e(ixI^S,idir)**2
  !   end do
  !   fld_R(ixI^S) = dsqrt(normgrad2(ixI^S))/(fld_kappa*w(ixI^S,iw_rho)*w(ixI^S,r_e))
  !
  !   !> Calculate the flux limiter, lambda
  !   !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
  !   fld_lambda(ixI^S) = (2+fld_R(ixI^S))/(6+3*fld_R(ixI^S)+fld_R(ixI^S)**2)
  !
  !   !> Diffusion coeficient c*lamda/ (rho*kappa)
  !   Diff_coef(ixI^S) = const_c/unit_velocity*fld_lambda(ixI^S)/(fld_kappa*w(ixI^S,iw_rho))
  !
  !   !> Calculate h
  !   !> This should be a single float !?!?!?!
  !   sys_h(ixI^S) = ps_dt/(2*dx(ixI^S,idir)**2)
  !
  !   !> Assemble vector
  !   do i =
  !     b(i) = (one - sys_h*(Diff_coef(i,j+1)+Diff_coef(i,j)))*w(i,j,r_e) + &
  !     sys_h*Diff_coef(i,j+1)*w(i,j+1,r_e) + sys_h*Diff_coef(i,j)*w(i,j-1,r_e) + ps_dt/(2*dt)*Rad_e_n(i,j)
  !   end do
  !
  ! end subroutine set_system_matrix


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

  subroutine grad_face(q,ixI^L,ix^L,idir,fdir,x,gradq)
    ! Compute the true gradient of a scalar q within ixL in direction idir ie :

    !AT CELL FACE

    !  - in cylindrical : d_/dr , (1/r)*d_/dth , d_/dz
    !  - in spherical   : d_/dr , (1/r)*d_/dth , (1/rsinth)*d_/dphi
    use mod_global_parameters

    integer :: ixI^L, ix^L, idir, fdir
    double precision :: q(ixI^S), gradq(ixI^S)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    integer :: jx^L, hx^L
    !-----------------------------------------------------------------------------
    !Only move to face in direction fdir, stay centered in the other
    hx^L=ix^L-kr(fdir,^D);
    gradq(ix^S)=(q(ix^S)-q(hx^S))/(x(ix^S,idir)-x(hx^S,idir))
    select case (typeaxial)
    case('slab') ! nothing to do
    case('cylindrical')
      if (idir==phi_) gradq(ix^S)=gradq(ix^S)*two/ (x(ix^S,r_)+x(hx^S,r_))
    case('spherical')
      if (idir==2   ) gradq(ix^S)=gradq(ix^S)*two/ (x(ix^S,r_)+x(hx^S,r_))
      if (idir==phi_) gradq(ix^S)=gradq(ix^S)*two/((x(ix^S,r_)+x(hx^S,r_))*dsin((x(ix^S,2)+x(hx^S,2))/two))
    case default
      call mpistop('Unknown geometry')
    end select
  end subroutine grad_face

end module mod_fld
