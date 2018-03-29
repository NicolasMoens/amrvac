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

    !> Timestep for supertimestepping on the radiation energy equation
    !> the supertimestep == fld_numdt * dt, so this is a ratio!!!
    integer, public :: fld_numdt = 10

    !> Switch different terms on/off
    !> Solve parabolic system using ADI (Diffusion)
    logical :: fld_Diffusion = .true.

    !> Use radiation force sourceterm
    logical :: fld_Rad_force = .true.

    !> Use Heating and Cooling sourceterms in e_
    logical :: fld_HeatCool = .true.

    !> Use Heating and Cooling sourceterms in E_rad
    logical :: fld_HeatCool_rad = .true.

    !> Use Photon Tiring term in E_rad
    logical :: fld_Phot_Tiring = .true.

    !> public methods
    public :: fld_add_source
    public :: fld_get_radflux
    public :: fld_get_fluxlimiter
    public :: fld_get_flux
    public :: fld_get_csound2

  contains

  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa, fld_mu, fld_split, fld_numdt, fld_Diffusion,&
    fld_Rad_force, fld_HeatCool, fld_HeatCool_rad, fld_Phot_Tiring

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

    !> Check if fld_numdt is not 1
    if (fld_numdt .lt. 2) call mpistop("fld_numdt should be an integer larger than 1")

    !> Make kappa dimensionless !!!STILL NEED TO MULTIPLY W RHO
    fld_kappa = fld_kappa*unit_time*unit_velocity*unit_density

    !> set radiation energy mod_variables
    r_e = var_set_radiation_energy()

    !> Dimensionless speed of light
    fld_speedofligt_0 = const_c/unit_velocity

    !> Dimensionless Boltzman constante sigma
    fld_sigma_0 = 5.67051d-5*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)

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

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      !> Begin by evolving the radiation energy field
      if (fld_Diffusion) then
        call Evolve_ADI(w, x, fld_numdt, ixI^L, ixO^L)
      endif

      !> Calculate the radiative flux using the FLD Approximation
      call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux, rad_pressure)

      if (fld_Rad_force) then
        do idir = 1,ndir
          !> Radiation force = kappa*rho/c *Flux
          radiation_force(ixI^S,idir) = fld_kappa*wCT(ixI^S,iw_rho)/fld_speedofligt_0*rad_flux(ixI^S, idir)

          !> Momentum equation source term
          w(ixI^S,iw_mom(idir)) = w(ixI^S,iw_mom(idir)) &
              + qdt * radiation_force(ixI^S,idir)
        enddo
      endif

      !> Get pressure
      call phys_get_pthermal(wCT,x,ixI^L,ixO^L,temperature)

      !> calc Temperature as p/rho
      !temperature(ixI^S)=(temperature(ixO^S)/wCT(ixO^S,iw_rho))
      temperature(ixO^S)=(temperature(ixO^S)/wCT(ixO^S,iw_rho)) !????

      !> Cooling = 4 pi kappa B = 4 kappa sigma T**4
      radiation_cooling(ixO^S) = 4.d0*fld_kappa*wCT(ixO^S,iw_rho)*fld_sigma_0*temperature(ixO^S)**4
      !> Heating = c kappa E_rad
      radiation_heating(ixO^S) = fld_speedofligt_0*fld_kappa*wCT(ixO^S,iw_rho)*wCT(ixO^S,iw_r_e)

      !> Write energy to file
      if (it == 0) open(1,file='energy_out5')
      write(1,222) it,global_time,w(5,5,iw_e)
      if (it == it_max) close(1)
      222 format(i8,2e15.5E3)

      !> Energy equation source terms
      if (fld_HeatCool) then
        w(ixO^S,iw_e) = w(ixO^S,iw_e) &
           + qdt * radiation_heating(ixO^S) &
           - qdt * radiation_cooling(ixO^S)
      endif

      !> Radiation Energy source term
      if (fld_HeatCool_rad) then
        w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) &
           - qdt * radiation_heating(ixO^S) &
           + qdt * radiation_cooling(ixO^S)
      endif

      !> Photon tiring
      if (fld_Phot_Tiring) then
        call divvector(wCT(ixI^S,iw_mom(:)),ixI^L,ixO^L,div_v)
        photon_tiring(ixO^S) = div_v(ixO^S)*rad_pressure(ixO^S)

        w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) &
           - qdt * photon_tiring(ixO^S)
      endif
    end if

  end subroutine fld_add_source


  !> Calculate fld flux limiter
  subroutine fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
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
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_flux(ixI^S, 1:ndim), rad_pressure(ixI^S)
    double precision :: fld_lambda(ixI^S), fld_R(ixI^S), normgrad2(ixI^S), f(ixI^S)
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
      rad_flux(ixI^S, idir) = -fld_speedofligt_0*fld_lambda(ixI^S)/(fld_kappa*w(ixI^S,iw_rho)) *grad_r_e(ixI^S,idir)
    end do

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    f(ixI^S) = fld_lambda(ixI^S) + fld_lambda(ixI^S)**2 * fld_R(ixI^S)**2
    f(ixI^S) = one/two*(one-f(ixI^S)) + one/two*(3*f(ixI^S) - one)
    rad_pressure(ixI^S) = f(ixI^S) * w(ixI^S, iw_r_e)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  ^   THIS IS NOT YET CORRECT   ^  !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine fld_get_radflux


  subroutine fld_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: f(ixI^S, nwflux)

    f(ixO^S, iw_r_e) = w(ixO^S,iw_mom(idim)) * w(ixO^S, iw_r_e)

  end subroutine fld_get_flux



  subroutine Evolve_ADI(w, x, w_max, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, w_max
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: E_m(ixI^S), E_n(ixI^S)
    double precision :: diag1(ixImax1,ixImax2),sub1(ixImax1,ixImax2),sup1(ixImax1,ixImax2),bvec1(ixImax1,ixImax2)
    double precision :: diag2(ixImax2,ixImax1),sub2(ixImax2,ixImax1),sup2(ixImax2,ixImax1),bvec2(ixImax2,ixImax1)
    double precision :: Evec1(ixImax1), Evec2(ixImax2)
    double precision :: dw, delta_x
    integer g, h, m, j

    E_n(ixI^S) = w(ixI^S,iw_r_e)
    E_m(ixI^S) = w(ixI^S,iw_r_e)

    !> WHY CAN'T I USE dx ?!?!?!?!?
    delta_x = min( (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)), (x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)) )

    do m = 1,w_max
      !> Set pseudotimestep
      dw = delta_x/4.d0*(max(x(ixOmax1,ixOmin2,1)-x(ixOmin1,ixOmin2,1),x(ixOmin1,ixOmax2,2)-x(ixOmin1,ixOmin2,2))/delta_x)**((m-one)/(w_max-one))

      !> Setup matrix and vector for sweeping in direction 1
      call make_matrix(x,w,dw,E_m,E_n,1,ixImax1,ixI^L, ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)

      do j = ixImin2,ixImax2
        Evec1 = E_m(:,j)
        call solve_tridiag(ixOmin1,ixOmax1,ixImin1,ixImax1,diag1(:,j),sub1(:,j),sup1(:,j),bvec1(:,j),Evec1)
        E_m(:,j) = Evec1
      enddo

      call ADI_boundary_conditions(ixI^L,E_m,w)

      !> Setup matrix and vector for sweeping in direction 2
      call make_matrix(x,w,dw,E_m,E_n,2,ixImax2,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
      do j = ixImin1,ixImax1
        Evec2 = E_m(j,:)
        call solve_tridiag(ixOmin2,ixOmax2,ixImin2,ixImax2,diag2(:,j),sub2(:,j),sup2(:,j),bvec2(:,j),Evec2)
        E_m(j,:) = Evec2
      enddo
      call ADI_boundary_conditions(ixI^L,E_m,w)

    enddo

    w(ixO^S,iw_r_e) = E_m(ixO^S)

  end subroutine Evolve_ADI


  subroutine make_matrix(x,w,dw,E_m,E_n,sweepdir,ixImax,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
    use mod_global_parameters

    integer, intent(in) :: sweepdir, ixImax
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw), dw
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(in) :: E_n(ixI^S), E_m(ixI^S)
    double precision, intent(out):: diag1(ixImax1,ixImax2),sub1(ixImax1,ixImax2),sup1(ixImax1,ixImax2),bvec1(ixImax1,ixImax2)
    double precision, intent(out):: diag2(ixImax2,ixImax1),sub2(ixImax2,ixImax1),sup2(ixImax2,ixImax1),bvec2(ixImax2,ixImax1)
    double precision :: fld_lambda(ixI^S), fld_R(ixI^S)
    double precision :: D_center(ixI^S), D(ixI^S,1:ndim), h, beta(ixImax), delta_x
    double precision :: grad_r_e(ixI^S, 1:ndim)
    integer :: idir,i,j

    !calculate lambda
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)

    !calculate diffusion coefficient
    D_center(ixI^S) = fld_speedofligt_0*fld_lambda(ixI^S)/(fld_kappa*w(ixI^S,iw_rho))

    !> Go from cell center to cell face
    do i = ixImin1+1, ixImax1
    do j = ixImin2+1, ixImax2
      D(i,j,1) = (D_center(i,j) + D_center(i-1,j))/two
      D(i,j,2) = (D_center(i,j) + D_center(i,j-1))/two
    enddo
    enddo
    D(ixImin1,ixImin2,:) = D_center(ixImin1,ixImin2)

    !calculate h
    delta_x = min( (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)), (x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)) )
    h = dw/(two*delta_x**two)

    !> Matrix depends on sweepingdirection
    if (sweepdir == 1) then
      !calculate matrix for sweeping in 1-direction
      do j = ixImin2,ixImax2
       !calculate beta
       do i = ixImin1,ixImax1-1
         beta(i) = one + dw/(two*dt) + h*(D(i+1,j,1)+D(i,j,1))
       enddo

       do i = ixImin1,ixImax1-1
         diag1(i,j) = beta(i)
         sub1(i+1,j) = -h*D(i+1,j,1)
         sup1(i,j) = -h*D(i+1,j,1)
         bvec1(i,j) = (1 - h*(D(i,j+1,2)+D(i,j,2)))*E_m(i,j) &
         + h*D(i,j+1,2)*E_m(i,j+1) + h*D(i,j,2)*E_m(i,j-1) + dw/(two*dt)*E_n(i,j)
       enddo

       !> Boundary conditions on matrix
       sub1(ixOmin1,j) = zero
       sup1(ixOmax1,j) = zero
       diag1(ixOmin1,j) = beta(ixOmin1) - h*D(ixOmin1,j,1)
       diag1(ixOmax1,j) = beta(ixOmax1) - h*D(ixOmax1+1,j,1)
      enddo

    elseif ( sweepdir == 2 ) then
      !calculate matrix for sweeping in 2-direction
      do j = ixImin1,ixImax1
       !calculate beta
       do i = ixImin2,ixImax2-1
         beta(i) = one + dw/(two*dt) + h*(D(j,i+1,2)+D(j,i,2))
       enddo

       do i = ixImin2,ixImax2-1
         diag2(i,j) = beta(i)
         sub2(i+1,j) = -h*D(j,i+1,2)
         sup2(i,j) = -h*D(j,i+1,2)
         bvec2(i,j) = (1 - h*(D(j+1,i,1)+D(j,i,1)))*E_m(j,i) &
         + h*D(j+1,i,1)*E_m(j+1,i) + h*D(j,i,1)*E_m(j-1,i) + dw/(two*dt)*E_n(j,i)
       enddo

       !> Boundary conditions on matrix
       sub2(ixOmin2,j) = zero
       sup2(ixOmax2,j) = zero
       diag2(ixOmin2,j) = beta(ixOmin2) - h*D(j,ixOmin2,2)
       diag2(ixOmax2,j) = beta(ixOmax2) - h*D(j,ixOmax2+1,2)

      enddo
    else
      call mpistop("sweepdirection unknown")
    endif

  end subroutine make_matrix


  subroutine solve_tridiag(ixOmin,ixOmax,ixImin,ixImax,diag,sub,sup,bvec,Evec)
    use mod_global_parameters
    implicit none

    integer, intent(in) :: ixOmin,ixOmax,ixImin,ixImax
    double precision, intent(in) :: diag(ixImax), bvec(ixImax)
    double precision, intent(in) :: sub(ixImax), sup(ixImax)
    double precision, intent(out) :: Evec(ixImax)
    double precision :: cp(ixImax), dp(ixImax)
    integer :: i

    ! initialize c-prime and d-prime
    cp(ixOmin) = sup(ixOmin)/diag(ixOmin)
    dp(ixOmin) = bvec(ixOmin)/diag(ixOmin)

    ! solve for vectors c-prime and d-prime
    do i = ixOmin ,ixOmax
      cp(i) = sup(i)/( diag(i)-cp(i-1)*sub(i))
      dp(i) = (bvec(i)-dp(i-1)*sub(i))/(diag(i)-cp(i-1)*sub(i))
    enddo

    ! initialize x
    Evec(ixOmax) = dp(ixOmax)

    ! solve for x from the vectors c-prime and d-prime
    do i = ixOmax-1, ixOmin, -1
      Evec(i) = dp(i)-cp(i)*Evec(i+1)
    end do

  end subroutine solve_tridiag



  subroutine ADI_boundary_conditions(ixI^L,E_m,w)
    use mod_global_parameters

    integer, intent(in) :: ixI^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: E_m(ixI^S)
    integer g, h

    !Edges
    !do g = 0,nghostcells-1
    do g = 0,nghostcells !> THIS IS VERRRYYYYY SJOEMEL-Y
      E_m(ixImin1+g,:) = w(ixImin1+nghostcells,:,iw_r_e)
      E_m(ixImax1-g,:) = w(ixImax1-nghostcells,:,iw_r_e)
      E_m(:,ixImin2+g) = w(:,ixImin2+nghostcells,iw_r_e)
      E_m(:,ixImax2-g) = w(:,ixImax2-nghostcells,iw_r_e)
    end do

    !Corners
    !do g = 0,nghostcells-1
    do g = 0,nghostcells !> THIS IS VERRRYYYYY SJOEMEL-Y
      !do h = 0, nghostcells-1
      do h = 0,nghostcells !> THIS IS VERRRYYYYY SJOEMEL-Y
        E_m(ixImin1+g,ixImax2-h) = w(ixImin1+nghostcells,ixImax2-nghostcells,iw_r_e)
        E_m(ixImax1-g,ixImax2-h) = w(ixImax1-nghostcells,ixImax2-nghostcells,iw_r_e)
        E_m(ixImin1+g,ixImin2+h) = w(ixImin1+nghostcells,ixImin2+nghostcells,iw_r_e)
        E_m(ixImax1-g,ixImin2+h) = w(ixImax1-nghostcells,ixImin2+nghostcells,iw_r_e)
      end do
    end do

  end subroutine ADI_boundary_conditions

  subroutine fld_get_csound2(w,x,ixI^L,ixO^L,hd_gamma,csound2)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw), hd_gamma
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)

    double precision :: pth(ixI^S), prad(ixI^S), rad_flux(ixI^S,1:ndim)

    call fld_get_radflux(w,x,ixI^L,ixO^L,rad_flux,prad)
    call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)
    csound2(ixO^S) = prad(ixO^S) + pth(ixO^S)
    csound2(ixO^S)=max(hd_gamma,4.d0/3.d0)*csound2(ixO^S)/w(ixO^S,iw_rho)

  end subroutine fld_get_csound2



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

end module mod_fld
