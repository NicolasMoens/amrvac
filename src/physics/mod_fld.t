!> Nicolas Moens
!> Module for including flux limited diffusion in hydrodynamics simulations
!> Based on Turner and stone 2001
module mod_fld
    implicit none

    !> source split or not
    logical :: fld_split = .false.

    !> Opacity per unit of unit_density
    double precision, public :: fld_kappa0 = 0.4d0

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

    !> Maximum amount of pseudotimesteps before trying something else
    integer, public :: fld_maxdw = 100

    !> Tolerance for bisection method for Energy sourceterms
    !> This is a percentage of the minimum of gas- and radiation energy
    double precision, public :: fld_bisect_tol = 1.d-3

    !> Tolerance for adi method for radiative Energy diffusion
    double precision, public :: fld_adi_tol = 1.d-2

    double precision :: fld_max_fracdt = 50.d0

    !> Switch different terms on/off
    !> Solve parabolic system using ADI (Diffusion)
    logical :: fld_Diffusion = .true.

    !> Use radiation force sourceterm
    logical :: fld_Rad_force = .true.

    !> Use Heating and Cooling sourceterms in e_
    logical :: fld_Energy_interact = .true.

    !> Let Vac advect radiative energy
    logical :: fld_Energy_advect = .true.

    !> Use constant Opacity
    logical :: fld_const_opacity = .false.

    logical :: fld_complete_diffusion_limit = .false.

    !> Boundary conditions for radiative Energy in ADI.
    character(len=8) :: fld_bound_min1 = 'periodic'
    character(len=8) :: fld_bound_max1 = 'periodic'
    character(len=8) :: fld_bound_min2 = 'periodic'
    character(len=8) :: fld_bound_max2 = 'periodic'

    !> Set Diffusion coefficient to unity
    logical :: fld_diff_testcase = .false.

    !> public methods
    !> these are called in mod_hd_phys
    public :: fld_add_source
    public :: fld_get_flux
    public :: fld_get_csound2
    !> these are used in specialvar_output
    public :: fld_get_radflux
    public :: fld_get_radpress
    public :: fld_get_fluxlimiter
    public :: fld_get_opacity

  contains

  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_mu, fld_split, fld_maxdw, fld_Diffusion,&
    fld_Rad_force, fld_Energy_interact, fld_Energy_advect, fld_bisect_tol, fld_diff_testcase,&
    fld_bound_min1, fld_bound_max1, fld_bound_min2, fld_bound_max2, fld_adi_tol, fld_max_fracdt,&
    fld_const_opacity, fld_complete_diffusion_limit

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
    if (fld_maxdw .lt. 2) call mpistop("fld_maxdw should be an integer larger than 1")

    !> Make kappa dimensionless !!!STILL NEED TO MULTIPLY W RHO
    fld_kappa0 = fld_kappa0*unit_time*unit_velocity*unit_density

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

    double precision :: fld_kappa(ixO^S)

    double precision :: rad_flux(ixO^S,1:ndim)
    double precision :: radiation_force(ixO^S,1:ndim)

    integer :: idir, i

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      !> Begin by evolving the radiation energy field
      if (fld_Diffusion) then
        call Evolve_E_rad(w, x, ixI^L, ixO^L)
        print*, it, " ######################"
      endif

      !> Add momentum sourceterms
      if (fld_Rad_force) then
        !> Calculate the radiative flux using the FLD Approximation
        call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
        call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux)

        do idir = 1,ndir
          !> Radiation force = kappa*rho/c *Flux
          radiation_force(ixO^S,idir) = fld_kappa(ixO^S)*wCT(ixO^S,iw_rho)/fld_speedofligt_0*rad_flux(ixO^S, idir)

          !> Momentum equation source term
          w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
              + qdt * radiation_force(ixO^S,idir)
        enddo
      endif

      !> Add energy sourceterms
      if (fld_Energy_interact) then
        call Energy_interaction(w, x, ixI^L, ixO^L)
      endif

      !> Write energy to file
      ! if (it == 0) open(1,file='energy_out0')
      ! write(1,222) it,global_time,w(3,3,iw_e),w(3,3,iw_r_e)
      ! if (it == it_max) close(1)
      ! 222 format(i8,3e15.5E3)
      ! print*, it, w(3,3,iw_e),w(3,3,r_e)

    end if
  end subroutine fld_add_source


  subroutine fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: fld_kappa(ixO^S)
    double precision :: temperature(ixI^S)
    double precision :: rho0,n

    integer :: i

    if (fld_const_opacity) then
       fld_kappa = fld_kappa0
    else

      !> CNST growth factor
      rho0 = 0.5d0 !> Take lower value of rho in domain
      n = 1.d0 !~0.5
      fld_kappa(ixO^S) = fld_kappa0*(one + 0.33d0*(w(ixO^S,iw_rho)/rho0)**n )

      ! !> Opacity bump
      ! rho0 = 0.5d0 !> Take central value of rho in domain
      ! n = 2.d0
      ! fld_kappa(ixO^S) = fld_kappa0*(one + n*dexp(-dlog(rho0/w(ixO^S,iw_rho))**two))

    endif
  end subroutine fld_get_opacity


  !> Calculate fld flux limiter
  subroutine fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: fld_R(ixO^S), fld_lambda(ixO^S)
    double precision :: fld_kappa(ixO^S)
    double precision ::  normgrad2(ixO^S)
    double precision :: grad_r_e(ixO^S)
    integer :: idir

    if (fld_complete_diffusion_limit) then
      fld_lambda = one/3.d0
    else
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixO^S) = zero
      do idir = 1,ndir
        call grad(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixO^S) = (2+fld_R(ixO^S))/(6+3*fld_R(ixO^S)+fld_R(ixO^S)**2)
    endif
  end subroutine fld_get_fluxlimiter


  !> Calculate Radiation Flux
  !> Returns Radiation flux and radiation pressure
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_flux(ixO^S, 1:ndim)
    double precision :: fld_kappa(ixO^S)
    double precision :: fld_lambda(ixO^S), fld_R(ixO^S), normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixO^S, 1:ndim)
    integer :: idir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixO^S) = zero
    do idir = 1,ndir
      call grad(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e(ixO^S,idir))
      normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S,idir)**2
    end do

    call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
    fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,r_e))

    !> Calculate the flux limiter, lambda
    !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
    if (fld_complete_diffusion_limit) then
      fld_lambda = one/3.d0
    else
      fld_lambda(ixO^S) = (two+fld_R(ixO^S))/(6.d0+3.d0*fld_R(ixO^S)+fld_R(ixO^S)**two)
    endif

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndir
      rad_flux(ixO^S, idir) = -fld_speedofligt_0*fld_lambda(ixO^S)/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)) *grad_r_e(ixO^S,idir)
    end do

    !rad_flux(:,ixOmax2-1, :) = rad_flux(:,ixOmax2-2, :)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! THIS IS SJOEMELY !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine fld_get_radflux


  !> Calculate Radiation Pressure
  !> Returns Radiation Pressure
  subroutine fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_pressure(ixO^S)
    double precision :: fld_kappa(ixO^S)
    double precision :: fld_lambda(ixO^S), fld_R(ixO^S), normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixO^S, 1:ndim)
    integer :: idir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixO^S) = zero

    do idir = 1,ndir
      call grad(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e(ixO^S,idir))
      normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S,idir)**two
    end do

    call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
    fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,r_e))

    !> Calculate the flux limiter, lambda
    !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
    if (fld_complete_diffusion_limit) then
      fld_lambda = one/3.d0
    else
      fld_lambda(ixO^S) = (two+fld_R(ixO^S))/(6.d0+3.d0*fld_R(ixO^S)+fld_R(ixO^S)**two)
    endif

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    f(ixO^S) = fld_lambda(ixO^S) + fld_lambda(ixO^S)**two * fld_R(ixO^S)**two
    f(ixO^S) = one/two*(one-f(ixO^S)) + one/two*(3.d0*f(ixO^S) - one)

    rad_pressure(ixO^S) = f(ixO^S) * w(ixO^S, iw_r_e)
  end subroutine fld_get_radpress


  subroutine fld_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: f(ixI^S, nwflux)

    if (fld_Energy_advect) then
      f(ixO^S, iw_r_e) = w(ixO^S,iw_mom(idim)) * w(ixO^S, iw_r_e)
    else
      f(ixO^S, iw_r_e) = zero
    endif
  end subroutine fld_get_flux


  subroutine Evolve_E_rad(w, x, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: E_new(ixI^S), E_old(ixI^S), ADI_Error
    double precision :: frac_grid
    integer :: w_max, frac_dt, i
    logical :: converged

    !> Reset E_new
    E_old(ixI^S) = w(ixI^S,iw_r_e)
    E_new(ixI^S) = w(ixI^S,iw_r_e)

    converged = .false.
    ADI_Error = bigdouble
    w_max = 1
    frac_grid = two
    frac_dt = 1

    do while (converged .eqv. .false.)

      !> Check if solution converged
      if (ADI_Error .lt. fld_adi_tol) then
        !> If converged in former loop, break loop
        converged = .true.
        goto 3000
      else
        !> If no convergence, adapt pseudostepping
        w_max = 2*w_max
        frac_grid = 2*frac_grid
        !> Reset E_new
        E_old(ixI^S) = w(ixI^S,iw_r_e)
        E_new(ixI^S) = w(ixI^S,iw_r_e)
      endif

      !> Evolve using ADI
      call Evolve_ADI(w, x, E_new, E_old, w_max, frac_grid, ixI^L, ixO^L)
      call Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error) !> SHOULD THIS BE DONE EVERY ITERATION???
      print*, w_max, ADI_Error

      !> If adjusting pseudostep doesn't work, divide the actual timestep in smaller parts
      if (w_max .gt. fld_maxdw) then
        !> use a smaller timestep than the hydrodynamical one
        call half_timestep_ADI(w, x, E_new, E_old, ixI^L, ixO^L, converged)
      endif
    enddo

    3000 w(ixO^S,iw_r_e) = E_new(ixO^S)
  end subroutine Evolve_E_rad


  subroutine half_timestep_ADI(w, x, E_new, E_old, ixI^L, ixO^L, converged)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(in) :: E_old(ixI^S)
    double precision, intent(out) :: E_new(ixI^S)
    logical, intent(inout) :: converged
    double precision :: frac_grid
    double precision :: E_loc(ixI^S)
    double precision :: saved_dt, ADI_Error
    integer :: i,  w_max, frac_dt

    saved_dt = dt
    ADI_Error = bigdouble
    frac_dt = 1
    5231 frac_dt = 2*frac_dt
    w_max = 1
    frac_grid = two

    if (frac_dt .gt. fld_max_fracdt) call mpistop("No convergence after halving timestep N times")
    dt = dt/frac_dt

    E_loc = E_old

    print*, "halving time"

    do i = 1,frac_dt
      !---------------------------------------------------------------
      do while (converged .eqv. .false.)
        !> Check if solution converged
        if (ADI_Error .lt. fld_adi_tol) then
          !> If converged in former loop, break loop
          converged = .true.
          goto 7895
        else
          !> If no convergence, adapt pseudostepping
          w_max = 2*w_max
          frac_grid = 2*frac_grid
        endif

        !> Evolve using ADI
        call Evolve_ADI(w, x, E_new, E_loc, w_max, frac_grid, ixI^L, ixO^L)
        call Error_check_ADI(w, x, E_new, E_loc, ixI^L, ixO^L, ADI_Error) !> SHOULD THIS BE DONE EVERY ITERATION???

        !> If adjusting pseudostep doesn't work, divide the actual timestep in smaller parts
        if (w_max .gt. fld_maxdw) goto 5231

      enddo
      !---------------------------------------------------------------
      E_loc = E_new
    enddo

    7895 dt = saved_dt
  end subroutine half_timestep_ADI


  subroutine Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S, 1:ndim), w(ixI^S, 1:nw)
    double precision, intent(in) :: E_new(ixI^S), E_old(ixI^S)
    double precision, intent(out) :: ADI_Error
    double precision :: LHS(ixO^S), RHS(ixO^S), D(ixI^S,1:ndim)
    integer :: jx1^L, hx1^L,jx2^L, hx2^L

    jx1^L=ixO^L+kr(1,^D);
    hx1^L=ixO^L-kr(1,^D);
    jx2^L=ixO^L+kr(2,^D);
    hx2^L=ixO^L-kr(2,^D);

    call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)

    !> LHS = dx^2/dt * (E_new - E_old)
    LHS(ixO^S) = (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))*&
    (x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/dt*&
    (E_new(ixO^S) - E_old(ixO^S))

    !> RHS = D1(E_+ - E) - D1(E - E_-) + D2(E_+ - E) - D2(E - E_-)
    RHS(ixO^S) = &
    D(jx1^S,1)*(E_new(jx1^S) - E_new(ixO^S)) &
    - D(ixO^S,1)*(E_new(ixO^S) - E_new(hx1^S)) &
    + D(jx2^S,2)*(E_new(jx2^S) - E_new(ixO^S)) &
    - D(ixO^S,2)*(E_new(ixO^S) - E_new(hx2^S))

    !ADI_Error = maxval(abs((RHS-LHS)/(E_old/dt)))!> Try mean value or smtn
    ADI_Error = sum(abs((RHS-LHS)/(E_old/dt)))/((ixOmax1-ixOmin1)*(ixOmax2-ixOmin2))
  end subroutine Error_check_ADI


  subroutine Evolve_ADI(w, x, E_new, E_old, w_max, frac_grid, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, w_max
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim), frac_grid
    double precision, intent(in) :: E_old(ixI^S)
    double precision, intent(out):: E_new(ixI^S)
    double precision :: E_m(ixI^S), E_n(ixI^S)
    double precision :: diag1(ixImax1,ixImax2),sub1(ixImax1,ixImax2),sup1(ixImax1,ixImax2),bvec1(ixImax1,ixImax2)
    double precision :: diag2(ixImax2,ixImax1),sub2(ixImax2,ixImax1),sup2(ixImax2,ixImax1),bvec2(ixImax2,ixImax1)
    double precision :: Evec1(ixImin1:ixImax1), Evec2(ixImin2:ixImax2)

    double precision :: a_1(ixOmin1:ixOmax1),b_1(ixOmin1:ixOmax1),c_1(ixOmin1:ixOmax1), d_1(ixOmin1:ixOmax1), x_1(ixOmin1:ixOmax1)
    double precision :: a_2(ixOmin2:ixOmax2),b_2(ixOmin2:ixOmax2),c_2(ixOmin2:ixOmax2), d_2(ixOmin2:ixOmax2), x_2(ixOmin2:ixOmax2)

    double precision :: dw, w0, w1
    integer :: m, j

    w0 = (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))*(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/frac_grid
    w1 = (x(ixOmax1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))*(x(ixOmin1,ixOmax2,2)-x(ixOmin1,ixOmin2,2))/frac_grid !4.d0

    E_m = E_new

    do m = 1,w_max
      E_n = E_old

      !> Set pseudotimestep
      dw = w0*(w1/w0)**((m-one)/(w_max-one))

      !> Setup matrix and vector for sweeping in direction 1
      call make_matrix(x,w,dw,E_m,E_n,1,ixImax1,ixI^L, ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
      do j = ixImin2,ixImax2
        ! Evec1(ixImin1:ixImax1) = E_m(ixImin1:ixImax1,j)
        ! call solve_tridiag(ixOmin1,ixOmax1,ixImin1,ixImax1,diag1(:,j),sub1(:,j),sup1(:,j),bvec1(:,j),Evec1)
        ! E_m(ixOmin1:ixOmax1,j) = Evec1(ixOmin1:ixOmax1)
        a_1(ixOmin1:ixOmax1) = sub1(ixOmin1:ixOmax1,j)
        b_1(ixOmin1:ixOmax1) = diag1(ixOmin1:ixOmax1,j)
        c_1(ixOmin1:ixOmax1) = sup1(ixOmin1:ixOmax1,j)
        d_1(ixOmin1:ixOmax1) = bvec1(ixOmin1:ixOmax1,j)
        call solve_tridiag_per(ixOmin1,ixOmax1,a_1,b_1,c_1,d_1,x_1)
        E_m(ixOmin1:ixOmax1,j) = x_1(ixOmin1:ixOmax1)
      enddo

      !call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)

      !> Setup matrix and vector for sweeping in direction 2
      call make_matrix(x,w,dw,E_m,E_n,2,ixImax2,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
      do j = ixImin1,ixImax1
        ! Evec2(ixImin2:ixImax2) = E_m(j,ixImin2:ixImax2)
        ! call solve_tridiag(ixOmin2,ixOmax2,ixImin2,ixImax2,diag2(:,j),sub2(:,j),sup2(:,j),bvec2(:,j),Evec2)
        ! E_m(j,ixOmin2:ixOmax2) = Evec2(ixOmin2:ixOmax2)
        a_2(ixOmin2:ixOmax2) = sub2(ixOmin2:ixOmax2,j)
        b_2(ixOmin2:ixOmax2) = diag2(ixOmin2:ixOmax2,j)
        c_2(ixOmin2:ixOmax2) = sup2(ixOmin2:ixOmax2,j)
        d_2(ixOmin2:ixOmax2) = bvec2(ixOmin2:ixOmax2,j)
        call solve_tridiag_per(ixOmin2,ixOmax2,a_2,b_2,c_2,d_2,x_2)
        E_m(j,ixOmin2:ixOmax2) = x_2(ixOmin2:ixOmax2)
      enddo

      !call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)
    enddo
    E_new = E_m
  end subroutine Evolve_ADI


  subroutine fld_get_diffcoef(w, x, ixI^L, ixO^L, D)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: D(ixI^S,1:ndim)
    double precision :: fld_kappa(ixO^S)
    double precision :: fld_lambda(ixO^S), fld_R(ixO^S)
    double precision :: D_center(ixI^S)
    integer :: idir,i,j

    if (fld_diff_testcase) then
      D = unit_length/unit_velocity
    else
      !> calculate lambda
      call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)

      !> set Opacity
      call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)

      !> calculate diffusion coefficient
      D_center(ixO^S) = fld_speedofligt_0*fld_lambda(ixO^S)/(fld_kappa(ixO^S)*w(ixO^S,iw_rho))

      !> Extrapolate lambda to ghostcells
      !> Edges
      !> To calculate the diffusion coefficient at the ghostcells, copy lambda from grid, but use correct kappa and rho
      do i = 0,nghostcells-1
        D_center(ixImin1+i,:) = D_center(ixImin1+nghostcells,:)
        D_center(ixImax1-i,:) = D_center(ixImax1-nghostcells,:)
        D_center(:,ixImin2+i) = D_center(:,ixImin2+nghostcells)
        D_center(:,ixImax2-i) = D_center(:,ixImax2-nghostcells)
      end do

      ! call Diff_boundary_conditions(ixI^L,ixO^L,D)

      !> Corners
      do i = 0,nghostcells-1
        do j = 0, nghostcells-1
          D_center(ixImin1+i,ixImax2-j) = D_center(ixImin1+nghostcells,ixImax2-nghostcells)
          D_center(ixImax1-i,ixImax2-j) = D_center(ixImax1-nghostcells,ixImax2-nghostcells)
          D_center(ixImin1+i,ixImin2+j) = D_center(ixImin1+nghostcells,ixImin2+nghostcells)
          D_center(ixImax1-i,ixImin2+j) = D_center(ixImax1-nghostcells,ixImin2+nghostcells)
        end do
      end do

      !> Go from cell center to cell face
      do i = ixImin1+1, ixImax1
      do j = ixImin2+1, ixImax2
        !D(i,j,1) = (D_center(i,j) + D_center(i-1,j))/two
        !D(i,j,2) = (D_center(i,j) + D_center(i,j-1))/two
        D(i,j,1) = (D_center(i,j) + D_center(i-1,j) + D_center(i,j+1) + D_center(i-1,j+1) + D_center(i,j-1) + D_center(i-1,j-1))/6.d0
        D(i,j,2) = (D_center(i,j) + D_center(i,j-1) + D_center(i+1,j) + D_center(i+1,j-1) + D_center(i-1,j) + D_center(i-1,j-1))/6.d0
      enddo
      enddo
      D(ixImin1,:,1) = D_center(ixImin1,:)
      D(:,ixImin2,1) = D_center(:,ixImin2)
      D(ixImin1,:,2) = D_center(ixImin1,:)
      D(:,ixImin2,2) = D_center(:,ixImin2)

      !D(:,ixImax2-2,:) = D(:,ixImax2-3,:)

    endif
  end subroutine fld_get_diffcoef


  subroutine make_matrix(x,w,dw,E_m,E_n,sweepdir,ixImax,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
    use mod_global_parameters

    integer, intent(in) :: sweepdir, ixImax
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw), dw
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(in) :: E_n(ixI^S), E_m(ixI^S)
    double precision, intent(out):: diag1(ixImin1:ixImax1,ixImin2:ixImax2),sub1(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: sup1(ixImin1:ixImax1,ixImin2:ixImax2),bvec1(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: diag2(ixImin2:ixImax2,ixImin1:ixImax1),sub2(ixImin2:ixImax2,ixImin1:ixImax1)
    double precision, intent(out):: sup2(ixImin2:ixImax2,ixImin1:ixImax1),bvec2(ixImin2:ixImax2,ixImin1:ixImax1)
    double precision :: D(ixI^S,1:ndim), h, beta(ixImax), delta_x
    integer :: idir,i,j

    call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)

    !calculate h
    if (sweepdir == 1) then
      delta_x = x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
      !delta_x = x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    elseif (sweepdir == 2) then
      !delta_x = x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
      delta_x = x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    endif
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
         bvec1(i,j) = (one - h*(D(i,j+1,2)+D(i,j,2)))*E_m(i,j) &
         + h*D(i,j+1,2)*E_m(i,j+1) + h*D(i,j,2)*E_m(i,j-1) + dw/(two*dt)*E_n(i,j)
       enddo

       !> Boundary conditions on matrix
       sub1(ixImin1,j) = zero
       sup1(ixImax1,j) = zero
       diag1(ixImin1,j) = beta(ixImin1) - h*D(ixImin1,j,1)
       diag1(ixImax1,j) = beta(ixImax1-1) - h*D(ixImax1,j,1)
       bvec1(ixImax1,j) = (one - h*(D(ixImax1,j+1,2)+D(ixImax1,j,2)))*E_m(ixImax1,j) &
       + h*D(ixImax1,j+1,2)*E_m(ixImax1,j+1) + h*D(ixImax1,j,2)*E_m(ixImax1,j-1) + dw/(two*dt)*E_n(ixImax1,j)

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
         bvec2(i,j) = (one - h*(D(j+1,i,1)+D(j,i,1)))*E_m(j,i) &
         + h*D(j+1,i,1)*E_m(j+1,i) + h*D(j,i,1)*E_m(j-1,i) + dw/(two*dt)*E_n(j,i)
       enddo

       !> Boundary conditions on matrix
       sub2(ixImin2,j) = zero
       sup2(ixImax2,j) = zero
       diag2(ixImin2,j) = beta(ixImin2) - h*D(j,ixImin2,2)
       diag2(ixImax2,j) = beta(ixImax2-1) - h*D(j,ixImax2,2)
       bvec2(ixImax2,j) = (one - h*(D(j+1,ixImax2,1)+D(j,ixImax2,1)))*E_m(j,ixImax2) &
       + h*D(j+1,ixImax2,1)*E_m(j+1,ixImax2) + h*D(j,ixImax2,1)*E_m(j-1,ixImax2) + dw/(two*dt)*E_n(j,ixImax2)
      enddo

    else
      call mpistop("sweepdirection unknown")
    endif
  end subroutine make_matrix


  subroutine solve_tridiag(ixOmin,ixOmax,ixImin,ixImax,diag,sub,sup,bvec,Evec)
    use mod_global_parameters

    integer, intent(in) :: ixOmin,ixOmax,ixImin,ixImax
    double precision, intent(in) :: diag(ixImin:ixImax), bvec(ixImin:ixImax)
    double precision, intent(in) :: sub(ixImin:ixImax), sup(ixImin:ixImax)
    double precision, intent(out) :: Evec(ixImin:ixImax)
    double precision :: cp(ixImin:ixImax), dp(ixImin:ixImax)
    integer :: i

    ! initialize c-prime and d-prime
    cp(ixImin) = sup(ixImin)/diag(ixImin)
    dp(ixImin) = bvec(ixImin)/diag(ixImin)

    ! solve for vectors c-prime and d-prime
    do i = ixImin+1 ,ixImax-1
      cp(i) = sup(i)/(diag(i)-cp(i-1)*sub(i))
      dp(i) = (bvec(i)-dp(i-1)*sub(i))/(diag(i)-cp(i-1)*sub(i))
    enddo
    dp(ixImax) = (bvec(ixImax)-dp(ixImax-1)*sub(ixImax))/(diag(ixImax)-cp(ixImax-1)*sub(ixImax))

    ! initialize x
    Evec(ixImax) = dp(ixImax)

    ! solve for x from the vectors c-prime and d-prime
    do i = ixImax-1, ixImin, -1
      Evec(i) = dp(i)-cp(i)*Evec(i+1)
    end do
  end subroutine solve_tridiag

  subroutine solve_tridiag_per(ixOmin,ixOmax,a_s,b_s,c_s,d_s,x_s)
    use mod_global_parameters
    integer, intent(in) :: ixOmin,ixOmax
    double precision, intent(in) :: a_s(ixOmin:ixOmax),b_s(ixOmin:ixOmax),c_s(ixOmin:ixOmax), d_s(ixOmin:ixOmax)
    double precision, intent(out) :: x_s(ixOmin:ixOmax)
    double precision ::  u_s(ixOmin:ixOmax),v_s(ixOmin:ixOmax),y_s(ixOmin:ixOmax), q_s(ixOmin:ixOmax)
    integer :: i

    u_s = zero
    v_s = zero

    u_s(ixOmin) = -b_s(ixOmin)
    u_s(ixOmax) = c_s(ixOmax)

    v_s(ixOmin) = one
    v_s(ixOmax) = -a_s(ixOmin)/b_s(ixOmin)

    call solve_tridiag2(ixOmin,ixOmax,a_s,b_s,c_s,d_s,y_s)
    call solve_tridiag2(ixOmin,ixOmax,a_s,b_s,c_s,q_s,u_s)

    do i = ixOmin,ixOmax
      x_s(i) = y_s(i) - ((v_s(i)*y_s(i))/(one + (v_s(i)*q_s(i))))*q_s(i)
    enddo
  end subroutine solve_tridiag_per

  subroutine solve_tridiag2(ixOmin,ixOmax,a_s,b_s,c_s,d_s,x_s)
    use mod_global_parameters
    integer, intent(in) :: ixOmin,ixOmax
    double precision, intent(in) :: a_s(ixOmin:ixOmax),b_s(ixOmin:ixOmax),c_s(ixOmin:ixOmax), d_s(ixOmin:ixOmax)
    double precision, intent(out) :: x_s(ixOmin:ixOmax)
    double precision :: m_s, bb_s(ixOmin:ixOmax), dd_s(ixOmin:ixOmax)
    integer :: i

    bb_s = b_s
    dd_s = d_s

    do i = ixOmin+1,ixOmax
      m_s = a_s(i)/b_s(i-1)
      bb_s = b_s(i) - m_s*c_s(i-1)
      dd_s = d_s(i) - m_s*d_s(i-1)
    enddo

    x_s(ixOmax) = d_s(ixOmax)/b_s(ixOmax)

    do i = ixOmax - 1, ixOmin, -1
      x_s(i) = (d_s(i)-c_s(i)*x_s(i+1))/b_s(i)
    enddo
  end subroutine solve_tridiag2

  subroutine ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: E_m(ixI^S)
    double precision :: fld_R(ixO^S), fld_lambda(ixO^S), fld_kappa(ixO^S)
    integer g, h, i, j

    select case (fld_bound_min1)
    case('periodic')
      E_m(ixImin1:ixOmin1-1,ixOmin2:ixOmax2) = E_m(ixOmax1-1:ixOmax1,ixOmin2:ixOmax2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      E_m(ixOmin1-1,:) = 2.d0*E_m(ixOmin1,:) - E_m(ixOmin1+1,:)
      E_m(ixImin1,:) = 2.d0*E_m(ixOmin1-1,:) - E_m(ixOmin1,:)
    case('fixed')
      E_m(ixImin1:ixOmin1-1,:) = w(ixImin1:ixOmin1-1,:,iw_r_e)
    case('mean')
      do i = ixImin2, ixImax2
        E_m(ixImin1:ixOmin1-1,i) = sum(E_m(ixOmin1:ixOmax1,i))/size(E_m(ixOmin1:ixOmax1,i))
      enddo
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max1)
    case('periodic')
      E_m(ixOmax1+1:ixImax1,ixOmin2:ixOmax2) = E_m(ixOmin1:ixOmin1+1,ixOmin2:ixOmax2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      E_m(ixOmax1+1,:) = 2.d0*E_m(ixOmax1,:) - E_m(ixOmax1-1,:)
      E_m(ixImax1,:) = 2.d0*E_m(ixOmax1+1,:) - E_m(ixOmax1,:)
    case('fixed')
      E_m(ixImax1:ixOmax1+1,:) = w(ixImax1:ixOmax1+1,:,iw_r_e)
    case('mean')
      do i = ixImin2, ixImax2
        E_m(ixOmax1+1:ixImax1,i) = sum(E_m(ixOmin1:ixOmax1,i))/size(E_m(ixOmin1:ixOmax1,i))
        print*, E_m(ixOmax1+1:ixImax1,i)
      enddo
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_min2)
    case('periodic')
      E_m(ixOmin1:ixOmax1,ixImin2:ixOmin2-1) = E_m(ixOmin1:ixOmax1,ixOmax2-1:ixOmax2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      E_m(:,ixOmin2-1) = 2.d0*E_m(:,ixOmin2) - E_m(:,ixOmin2+1)
      E_m(:,ixImin2) = 2.d0*E_m(:,ixOmin2-1) - E_m(:,ixOmin2)
    case('fixed')
      E_m(ixOmin1:ixOmax1,ixImin2:ixOmin2-1) = w(ixOmin1:ixOmax1,ixImin2:ixOmin2-1,iw_r_e)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max2)
    case('periodic')
      E_m(ixOmin1:ixOmax1,ixImax2:ixOmax2+1) = E_m(ixOmin1:ixOmax1,ixOmin2+1:ixOmin2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      E_m(:,ixOmax2+1) = 2.d0*E_m(:,ixOmax2) - E_m(:,ixOmax2-1)
      E_m(:,ixImax2) = 2.d0*E_m(:,ixOmax2+1) - E_m(:,ixOmax2)
    case('fixed')
      E_m(:,ixImax2:ixOmax2+1) = w(:,ixImax2:ixOmax2+1,iw_r_e)
    case('cons')
      !> Conservation law
      call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
      call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
      do i = ixImax2,ixOmax2+1
        E_m(:,i+1) = ((w(:,i-1,iw_mom(2))*w(:,i-1,r_e)/w(:,i+1,iw_rho)&
         - fld_lambda(:,ixOmax2)*fld_speedofligt_0/(w(:,i-1,iw_rho)*fld_kappa(:,ixOmax2)) &
         * (w(:,i-2,r_e) - w(:,i,r_e))/abs(x(:,i+1,2)- x(:,i-1,2))&
         - w(:,i,iw_mom(2))*w(:,i,r_e)/w(:,i,iw_rho))&
         * w(:,i,iw_rho)*fld_kappa(:,ixOmax2)/(fld_lambda(:,ixOmax2)*fld_speedofligt_0)*abs(x(:,i+1,2)- x(:,i-1,2)))&
         + w(:,i-1,r_e)
         do j = ixImin1,ixImax1
           E_m(j,i+1) = min(w(j,i+1,r_e), w(j,i,r_e))
         enddo
      enddo
    case default
      call mpistop("ADI boundary not defined")
    end select

    ! !Corners
    ! do g = 0,nghostcells-1
    !   do h = 0, nghostcells-1
    !     E_m(ixImin1+g,ixImax2-h) = w(ixImin1+nghostcells,ixImax2-nghostcells,iw_r_e)
    !     E_m(ixImax1-g,ixImax2-h) = w(ixImax1-nghostcells,ixImax2-nghostcells,iw_r_e)
    !     E_m(ixImin1+g,ixImin2+h) = w(ixImin1+nghostcells,ixImin2+nghostcells,iw_r_e)
    !     E_m(ixImax1-g,ixImin2+h) = w(ixImax1-nghostcells,ixImin2+nghostcells,iw_r_e)
    !   end do
    ! end do
  end subroutine ADI_boundary_conditions


  subroutine Diff_boundary_conditions(ixI^L,ixO^L,D)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(inout) :: D(ixI^S)
    double precision :: Dmn1(ixI^S),Dmx1(ixI^S),Dmn2(ixI^S),Dmx2(ixI^S)

    select case (fld_bound_min1)
    case('periodic')
      D(ixImin1:ixOmin1-1,:) = D(ixOmax1-1:ixOmax1,:)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(ixOmin1-1,:) = 2.d0*D(ixOmin1,:) - D(ixOmin1+1,:)
      D(ixImin1,:) = 2.d0*D(ixOmin1-1,:) - D(ixOmin1,:)
    case('fixed')
      if (it==0) then
        Dmn1 = D
      endif
      D(ixImin1:ixOmin1-1,:) = Dmn1(ixImin1:ixOmin1-1,:)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max1)
    case('periodic')
      D(ixImax1:ixOmax1+1,:) = D(ixOmin1+1:ixOmin1,:)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(ixOmax1+1,:) = 2.d0*D(ixOmax1,:) - D(ixOmax1-1,:)
      D(ixImax1,:) = 2.d0*D(ixOmax1+1,:) - D(ixOmax1,:)
    case('fixed')
      if (it==0) then
        Dmx1 = D
      endif
      D(ixImax1:ixOmax1+1,:) = Dmx1(ixImax1:ixOmax1+1,:)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_min2)
    case('periodic')
      D(:,ixImin2:ixOmin2-1) = D(:,ixOmax2-1:ixOmax2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(:,ixOmin2-1) = 2.d0*D(:,ixOmin2) - D(:,ixOmin2+1)
      D(:,ixImin2) = 2.d0*D(:,ixOmin2-1) - D(:,ixOmin2)
    case('fixed')
      if (it==0) then
        Dmn2 = D
      endif
      D(:,ixImin2:ixOmin2-1) = Dmn2(:,ixImin2:ixOmin2-1)
      ! D(:,ixImin2+1) = Dmn2(:,ixOmin2)
      ! D(:,ixImin2) = Dmn2(:,ixOmin2)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max2)
    case('periodic')
      D(:,ixImax2:ixOmax2+1) = D(:,ixOmin2+1:ixOmin2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(:,ixOmax2+1) = 2.d0*D(:,ixOmax2) - D(:,ixOmax2-1)
      D(:,ixImax2) = 2.d0*D(:,ixOmax2+1) - D(:,ixOmax2)
    case('fixed')
      if (it==0) then
        Dmx2 = D
      endif
      D(:,ixImax2:ixOmax2+1) = Dmx2(:,ixImax2:ixOmax2+1)
    case default
      call mpistop("ADI boundary not defined")
    end select
  end subroutine Diff_boundary_conditions


  subroutine Energy_interaction(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: fld_kappa(ixO^S)
    double precision :: rad_pressure(ixO^S)
    double precision :: temperature(ixI^S), div_v(ixI^S), vel(ixI^S,1:ndim)
    double precision :: a1(ixO^S), a2(ixO^S), a3(ixO^S)
    double precision :: c0(ixO^S), c1(ixO^S)
    double precision :: e_gas(ixO^S), E_rad(ixO^S)

    integer :: i,j,idir

    !> Calculate the radiative flux using the FLD Approximation
    call fld_get_radpress(w,x,ixI^L,ixO^L,rad_pressure)

    !> Get pressure
    call phys_get_pthermal(w,x,ixI^L,ixO^L,temperature)

    !> calc Temperature as p/rho
    temperature(ixO^S)=temperature(ixO^S)/w(ixO^S,iw_rho)

    !> set Opacity
    call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)

    !> calc photon tiring term
    do idir=1,ndim
      vel(ixI^S,idir)= w(ixI^S,iw_mom(idir))/w(ixI^S,iw_rho)
    enddo

    call divvector(vel,ixI^L,ixO^L,div_v)

    e_gas(ixO^S) = w(ixO^S,iw_e)
    E_rad(ixO^S) = w(ixO^S,iw_r_e)

    !> Calculate coefficients for polynomial
    a1(ixO^S) = 4*fld_kappa(ixO^S)*w(ixO^S,iw_rho)*fld_sigma_0*(temperature(ixO^S)/e_gas(ixO^S))**4.d0*dt
    a2(ixO^S) = fld_speedofligt_0*fld_kappa(ixO^S)*w(ixO^S,iw_rho)*dt
    a3(ixO^S) = div_v(ixO^S)*rad_pressure(ixO^S)/E_rad(ixO^S)*dt

    c0(ixO^S) = ((one + a1(ixO^S) + a3(ixO^S))*e_gas(ixO^S) + a2(ixO^S)*E_rad(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))
    c1(ixO^S) = (one + a1(ixO^S) + a3(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))

    !> Loop over every cell for bisection method
    do i = ixOmin1,ixOmax1
    do j =  ixOmin2,ixOmax2
      call Bisection_method(e_gas(i,j), E_rad(i,j), c0(i,j), c1(i,j))
    enddo
    enddo

    !> Update gas-energy in w
    w(ixO^S,iw_e) = e_gas(ixO^S)

    !> Calculate new radiation energy
    !> Get new pressure
    call phys_get_pthermal(w,x,ixI^L,ixO^L,temperature)

    !> calc new Temperature as p/rho
    temperature(ixO^S)=(temperature(ixO^S)/w(ixO^S,iw_rho))

    !> Update a1
    a1(ixO^S) = 4*fld_kappa(ixO^S)*w(ixO^S,iw_rho)*fld_sigma_0*(temperature(ixO^S)/e_gas(ixO^S))**4.d0*dt

    !> advance E_rad
    E_rad(ixO^S) = (a1*e_gas(ixO^S)**4.d0 + E_rad(ixO^S))/(one + a2 + a3)

    !> Update rad-energy in w
    w(ixO^S,iw_r_e) = E_rad(ixO^S)
  end subroutine Energy_interaction


  subroutine Bisection_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: bisect_a, bisect_b, bisect_c

    bisect_a = zero
    bisect_b = min(abs(c0/c1),abs(c0)**(1.d0/4.d0))

    do while (abs(Polynomial_Bisection(bisect_b, c0, c1)-Polynomial_Bisection(bisect_a, c0, c1))&
       .ge. fld_bisect_tol*min(e_gas,E_rad))
      bisect_c = (bisect_a + bisect_b)/two

      if (Polynomial_Bisection(bisect_a, c0, c1)*&
      Polynomial_Bisection(bisect_b, c0, c1) .le. zero) then

        if (Polynomial_Bisection(bisect_a, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .le. zero) then
          bisect_b = bisect_c
        elseif (Polynomial_Bisection(bisect_b, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .le. zero) then
          bisect_a = bisect_c
        else
          call mpistop("Problem with fld bisection method")
        endif
      else
        bisect_a = e_gas
        bisect_b = e_gas
        print*, "IGNORING ENERGY GAS-RAD EXCHANGE "
      endif
    enddo

    e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method


  function Polynomial_Bisection(e_gas, c0, c1) result(pol_result)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: pol_result

    pol_result = e_gas**4.d0 + c1*e_gas - c0
  end function Polynomial_Bisection


  subroutine fld_get_csound2(w,x,ixI^L,ixO^L,hd_gamma,csound2)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw), hd_gamma
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision :: pth(ixI^S), prad(ixO^S)

    call fld_get_radpress(w,x,ixI^L,ixO^L,prad)
    call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)

    csound2(ixO^S) = pth(ixO^S) + prad(ixO^S)
    csound2(ixO^S) = max(hd_gamma,4.d0/3.d0)*csound2(ixO^S)/w(ixO^S,iw_rho)
  end subroutine fld_get_csound2


  subroutine grad(q,ixI^L,ixO^L,idir,x,gradq)
    ! Compute the true gradient of a scalar q within ixL in direction idir ie :
    !  - in cylindrical : d_/dr , (1/r)*d_/dth , d_/dz
    !  - in spherical   : d_/dr , (1/r)*d_/dth , (1/rsinth)*d_/dphi
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idir
    double precision, intent(in) :: q(ixI^S), x(ixI^S,1:ndim)
    double precision, intent(out) ::  gradq(ixO^S)

    integer :: jx^L, hx^L
    !-----------------------------------------------------------------------------
    jx^L=ixO^L+kr(idir,^D);
    hx^L=ixO^L-kr(idir,^D);

    gradq(ixO^S)=(q(jx^S)-q(hx^S))/(x(jx^S,idir)-x(hx^S,idir))

    select case (typeaxial)
    case('slab') ! nothing to do
    case('cylindrical')
      if (idir==phi_) gradq(ixO^S)=gradq(ixO^S)/ x(ixO^S,r_)
    case('spherical')
      if (idir==2   ) gradq(ixO^S)=gradq(ixO^S)/ x(ixO^S,r_)
      if (idir==phi_) gradq(ixO^S)=gradq(ixO^S)/(x(ixO^S,r_)*dsin(x(ixO^S,2)))
    case default
      call mpistop('Unknown geometry')
    end select
  end subroutine grad

end module mod_fld
