!> Module for including dust species, which interact with the gas through a drag
!> force
module mod_dust
  use mod_global_parameters, only: std_len
  use mod_physics

  implicit none
  private

  !> The number of dust species
  integer, public, protected      :: dust_n_species = 0

  integer, protected              :: gas_rho_ = -1
  integer, allocatable, protected :: gas_mom(:)
  integer, protected              :: gas_e_   = -1

  !> Mean molecular weight of gas molecules
  double precision, protected, public :: gas_mu = -huge(1.0d0)

  !> Indices of the dust densities
  integer, allocatable, public, protected :: dust_rho(:)

  !> Indices of the dust momentum densities
  integer, allocatable, public, protected :: dust_mom(:, :)

  !> Size of each dust species
  double precision, allocatable, public :: dust_size(:)

  !> Internal density of each dust species
  double precision, allocatable, public :: dust_density(:)

  !> Dust temperature (if dust_temperature_type is constant)
  double precision :: dust_temperature = -1.0d0

  !> If dust_temperature_type is stellar, it will be calculated according to Tielens (2005),
  !> eqn. 5.44 using a stellar luminosity in solar luminosities
  double precision :: dust_stellar_luminosity = -1.0d0

  !> Set small dust densities to zero to avoid numerical problems
  logical :: dust_small_to_zero = .false.

  !> Minimum dust density
  double precision :: dust_min_rho = -1.0d0

  !> TODO: 1. Introduce this generically in advance, 2: document
  logical :: dust_source_split = .false.

  !> What type of dust drag force to use. Can be 'Kwok', 'sticking', 'linear',or 'none'.
  character(len=std_len) :: dust_method = 'Kwok'

  !> Can be 'graphite' or 'silicate', affects the dust temperature
  character(len=std_len) :: dust_species = 'graphite'

  !> Determines the dust temperature, can be 'constant', 'ism', or 'stellar'
  character(len=std_len) :: dust_temperature_type = 'constant'

  ! Public methods
  public :: dust_init
  public :: dust_get_dt
  public :: dust_get_flux
  public :: dust_get_cmax
  public :: dust_get_flux_prim
  public :: dust_get_cmax_prim
  public :: dust_add_source
  public :: dust_to_conserved
  public :: dust_to_primitive
  public :: dust_check_params

contains

  subroutine dust_init(g_rho, g_mom, g_energy)
    use mod_global_parameters
    integer, intent(in) :: g_rho
    integer, intent(in) :: g_mom(ndir)
    integer, intent(in) :: g_energy ! Negative value if not present
    integer             :: n, idir
    character(len=2)    :: dim

    call dust_read_params(par_files)

    allocate(gas_mom(ndir))
    gas_rho_ = g_rho
    gas_mom  = g_mom
    gas_e_   = g_energy

    allocate(dust_size(dust_n_species))
    allocate(dust_density(dust_n_species))
    dust_size(:) = -1.0d0
    dust_density(:) = -1.0d0

    allocate(dust_rho(dust_n_species))
    allocate(dust_mom(ndir, dust_n_species))

    ! Set index of dust densities
    do n = 1, dust_n_species
      dust_rho(n) = var_set_fluxvar("rhod", "rhod", n)
    end do

    ! Dust momentum
    do idir = 1, ndir
      write(dim, "(I0,A)") idir, "d"
      do n = 1, dust_n_species
        dust_mom(idir, n) = var_set_fluxvar("m"//dim, "v"//dim, n)
      end do
    end do

  end subroutine dust_init

  !> Read this module"s parameters from a file
  subroutine dust_read_params(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /dust_list/ dust_n_species, dust_min_rho, gas_mu, dust_method,&
        dust_small_to_zero, dust_source_split, dust_temperature,&
        dust_temperature_type

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, dust_list, end=111)
111   close(unitpar)
    end do

  end subroutine dust_read_params

  subroutine dust_check_params()
    if (gas_mu <= 0.0d0) call mpistop ("Dust error: gas_mu (molecular weight)"//&
         "negative or not set")

    if (dust_temperature_type == "constant") then
       if (dust_temperature < 0.0d0) then
          call mpistop("Dust error: dust_temperature < 0 or not set")
       end if
    else if (dust_temperature_type == "stellar") then
       if (dust_stellar_luminosity < 0.0d0) then
          call mpistop("Dust error: dust_stellar_luminosity < 0 or not set")
       end if
    end if

    if (any(dust_size < 0.0d0)) call mpistop(&
       "Dust error: any(dust_size < 0) or not set")

    if (any(dust_density < 0.0d0)) call mpistop(&
       "Dust error: any(dust_density < 0) or not set")
  end subroutine dust_check_params

  subroutine dust_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    integer                         :: n, idir

    do n = 1, dust_n_species
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
            n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            dust_rho(n)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
            n))
      end do
    end do
  end subroutine dust_to_conserved

  subroutine dust_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    integer                         :: n, idir

    do n = 1, dust_n_species
      ! Convert momentum to velocity
      do idir = 1, ndir
        where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) > dust_min_rho)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
              n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
              n)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n))
        elsewhere
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir, n)) = 0
        end where
      end do
    end do
  end subroutine dust_to_primitive

  subroutine dust_get_flux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
    double precision, intent(inout) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    integer                         :: n, idir

    do n = 1, dust_n_species
      where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) > dust_min_rho)
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, dust_mom(idim, n))
      elsewhere             ! TODO: remove?
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) = 0.0d0
      end where

      do idir = 1, ndir
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
            n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
            n)) * get_vdust(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, idim, n)
      end do
    end do
  end subroutine dust_get_flux

  subroutine dust_get_flux_prim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
    double precision, intent(inout) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    integer                         :: n, idir

    do n = 1, dust_n_species
      where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) > dust_min_rho)
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, dust_mom(idim, n))*w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, dust_rho(n))
      elsewhere             ! TODO: remove?
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) = 0.0d0
      end where

      do idir = 1, ndir
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
            n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
            n)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            dust_rho(n)) * get_vdust_prim(w, ixImin1,ixImin2,ixImax1,ixImax2,&
            ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, n)
      end do
    end do
  end subroutine dust_get_flux_prim

  function get_vdust(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim, n
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) > dust_min_rho)
      vdust = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idim,&
          n)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n))
    elsewhere
      vdust = 0.0d0
    end where
  end function get_vdust

  function get_vdust_prim(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim, n
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) > dust_min_rho)
      vdust = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idim, n))
    elsewhere
      vdust = 0.0d0
    end where
  end function get_vdust_prim

  ! Force dust density to zero if dust_rho <= dust_min_rho
  subroutine set_dusttozero(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,  wCT,  w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    integer                         :: n, idir

    do n = 1, dust_n_species
      where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) <= dust_min_rho)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n)) = 0.0d0
      end where

      do idir = 1, ndir
        where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            dust_rho(n)) <= dust_min_rho)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir, n)) = 0.0d0
        end where
      end do
    end do
  end subroutine set_dusttozero

  ! Calculate drag force based on Epstein's or Stokes' law
  ! From Kwok 1975, page 584 (between eqn 8 and 9)
  subroutine get_3d_dragforce(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x, fdrag, ptherm, vgas)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(out)   :: fdrag(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndir, 1:dust_n_species)
    double precision, intent(in)    :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2),&
        vgas(ixImin1:ixImax1,ixImin2:ixImax2, ndir)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: vt2,&
        deltav, fd, vdust
    double precision                   :: alpha_T(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:dust_n_species)
    integer                            :: n, idir
    double precision                   :: K

    vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 3.0d0*ptherm(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, gas_rho_)

    select case( TRIM(dust_method) )
    case ('Kwok') ! assume sticking coefficient equals 0.25

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              dust_rho(n)) > dust_min_rho)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, dust_mom(idir, n)) / w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (vgas(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

            ! 0.75 from sticking coefficient
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = 0.75d0*w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, dust_rho(n))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                gas_rho_)*deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) / (dust_density(n) * dust_size(n))

            ! 0.75 from spherical grainvolume
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = -fd(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)*0.75d0*dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) + deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2)
          elsewhere
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
          end where
          fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir, n) = fd(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
        end do
      end do

    case ('sticking') !Calculate sticking coefficient based on the gas and dust temperatures
      !  Equation from Decin et al. 2006
      if (gas_e_ < 0) call mpistop("dust sticking requires gas energy")

      call get_sticking(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, alpha_T, ptherm)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n))>dust_min_rho)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,dust_mom(idir, n)) / w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (vgas(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = &
               (one-alpha_T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               n)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                dust_rho(n))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                gas_rho_) * deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) / (dust_density(n)*dust_size(n))
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = -fd(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) * 0.75d0 * dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) + deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2)
          else where
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
          end where
          fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir,n) = fd(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
        end do
      end do
    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      K = 3.4d5 / dust_n_species
      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n))>dust_min_rho)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,dust_mom(idir, n))/w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (vgas(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = &
               -K*deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          else where
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
          end where
          fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir,n) = fd(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
        end do
      end do
    case('none')
      fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2, :, :) = 0.0d0
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

  end subroutine get_3d_dragforce

  !> Get sticking coefficient
  !>
  !> Assume cgs units, and use of convert factors for conversion
  !> Equation from Decin et al. 2006
  subroutine get_sticking(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, alpha_T, ptherm)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(out) :: alpha_T(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:dust_n_species)
    double precision, intent(in)  :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision              :: Tgas(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                       :: n

    call get_tdust(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, alpha_T)

    Tgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (ptherm(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w_convert_factor(gas_e_)*mH_cgs) / (w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, gas_rho_) * w_convert_factor(gas_rho_) * kB_cgs)

    do n = 1, dust_n_species
      alpha_T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         n) =  max(0.35d0 * exp(-sqrt((Tgas(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) + alpha_T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         n))/5.0d2))+0.1d0, smalldouble)
    end do
  end subroutine get_sticking

  !> Returns dust temperature (in K), either as constant or based on equ. 5.41,
  !> 5.42 and 5.44 from Tielens (2005)
  !>
  !> Note that this calculation assumes cgs!!!! with conversion between physical
  !> and scaled quantities done through the convert factors!!!!
  !>
  !> It takes as input the stellar luminosoity in solar units and/or a fixed
  !> dust temperature in Kelvin
  subroutine get_tdust(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, Td)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(out) :: Td(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:dust_n_species)
    double precision              :: G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer                       :: n

    select case( trim(dust_temperature_type) )
    case( 'constant' )
      Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2, :) = dust_temperature
    case( 'ism' )
      select case( trim(dust_species) )
      case( 'graphite' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              n) = 15.8d0*((0.0001d0/(dust_size(n)*length_convert_factor))**&
             0.06d0)
        end do
      case( 'silicate' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              n) = 13.6d0*((0.0001d0/(dust_size(n)*length_convert_factor))**&
             0.06d0)
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select
    case( 'stellar' )
      select case( trim(typeaxial) )
      case( 'spherical' )
        G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, 1)*length_convert_factor, smalldouble)
      case( 'cylindrical' )
        G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(dsqrt(sum(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,:)**2,dim=ndim+1))*length_convert_factor,&
            smalldouble)
      end select

      G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
         2.1d4*(dust_stellar_luminosity/1.0d8)*((3.0857d17/G0(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))**2)

      select case( trim(dust_species) )
      case( 'graphite' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              n) = 61.0d0*((0.0001d0/(dust_size(n)*length_convert_factor))**&
             0.06d0) *(G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**(one/5.8d0))
        end do
      case( 'silicate' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              n) = 50.0d0*((0.0001d0/(dust_size(n)*length_convert_factor))**&
             0.06d0) *(G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**(one/6.0d0))
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select
    case default
      call mpistop( "=== Dust temperature undetermined===" )
    end select

  end subroutine get_tdust

  !> w[iw]= w[iw]+qdt*S[wCT,  x] where S is the source based on wCT within ixO
  subroutine dust_add_source(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, wCT,w, x, qsourcesplit, active)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    double precision :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2),&
        vgas(ixImin1:ixImax1,ixImin2:ixImax2, ndir)
    double precision :: fdrag(ixImin1:ixImax1,ixImin2:ixImax2, ndir,&
        dust_n_species)
    integer          :: n, idir

    select case( TRIM(dust_method) )
    case( 'none' )
      !do nothing here
    case default !all regular dust methods here
      if (qsourcesplit .eqv. dust_source_split) then
        active = .true.

        call phys_get_pthermal(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2,&
            ixOmin1,ixOmin2,ixOmax1,ixOmax2, ptherm)
        do idir=1,ndir
          vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,gas_mom(idir))/wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,gas_rho_)
        end do

        call get_3d_dragforce(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2, wCT, x, fdrag, ptherm, vgas)
        fdrag = fdrag * qdt

        do idir = 1, ndir

          do n = 1, dust_n_species
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                gas_mom(idir))  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                gas_mom(idir)) + fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir,&
                n)

            if (gas_e_ > 0) then
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, gas_e_) = w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2, gas_e_) + (wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2, gas_mom(idir)) / wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2, gas_rho_)) * fdrag(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2, idir, n)
            end if

            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
                n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idir,&
                n)) - fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir, n)
          end do
        end do

        if (dust_small_to_zero) then
          call set_dusttozero(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,  wCT,  w, x)
        end if
      endif
    end select

  end subroutine dust_add_source

  !> Get dt related to dust and gas stopping time (Laibe 2011)
  subroutine dust_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:2)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2),&
        vgas(ixImin1:ixImax1,ixImin2:ixImax2, ndir)
    double precision, dimension(1:dust_n_species)        :: dtdust
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)         :: vt2, deltav, tstop, vdust
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:dust_n_species) :: alpha_T
    double precision                           :: K
    integer                                    :: n, idir

    call phys_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, ptherm)
    do idir = 1, ndir
      vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,gas_mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         gas_rho_)
    end do
    select case( TRIM(dust_method) )

    case( 'Kwok' ) ! assume sticking coefficient equals 0.25
      dtdust(:) = bigdouble

      vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 3.0d0*ptherm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, gas_rho_)

      ! Tgas, mu = mean molecular weight
      ptherm(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = ( ptherm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * w_convert_factor(gas_e_) * mH_cgs*gas_mu) / &
         (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          gas_rho_) * w_convert_factor(gas_rho_)*kB_cgs)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n))>dust_min_rho)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,dust_mom(idir, n))/w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (vgas(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = &
               4.0d0*(dust_density(n)*dust_size(n))/ &
               (3.0d0*(0.75d0)*dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) + deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)**2)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                dust_rho(n)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, gas_rho_)))
          else where
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
              dtdust(n))
        end do
      end do

      dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case( 'sticking' ) !Calculate sticking coefficient based on the gas temperature
      dtdust(:) = bigdouble

      vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 3.0d0*ptherm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, gas_rho_)

      ! Sticking coefficient
      call get_sticking(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, alpha_T, ptherm)

      ! Tgas, mu = mean molecular weight
      ptherm(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = ( ptherm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*w_convert_factor(gas_e_) * mH_cgs*gas_mu) / &
         (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          gas_rho_) * w_convert_factor(gas_rho_)*kB_cgs)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n))>dust_min_rho)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,dust_mom(idir, n))/w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (vgas(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = &
               4.0d0*(dust_density(n)*dust_size(n))/ &
               (3.0d0*(one-alpha_T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               n))*dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) + deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)**2)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                dust_rho(n)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, gas_rho_)))
          else where
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
              dtdust(n))
        end do
      end do

      dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      K = 3.4d5/dust_n_species
      dtdust(:) = bigdouble

      do n = 1, dust_n_species
        where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_rho(n))>dust_min_rho)
          tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = (w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, dust_rho(n))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              gas_rho_))/ (K*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              dust_rho(n)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, gas_rho_)))
        else where
          tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = bigdouble
        end where

        dtdust(n) = min(minval(tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
            dtdust(n))
      end do

      dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)
    case('none')
      ! no dust timestep
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

    if (dtnew < dtmin) then
      write(unitterm,*)"-------------------------------------"
      write(unitterm,*&
         )"Warning: found DUST related time step too small! dtnew=", dtnew
      write(unitterm,*)"on grid with index:", saveigrid," grid level=",&
          node(plevel_, saveigrid)
      write(unitterm,*)"grid corners are=",rnode(rpxmin1_, saveigrid),&
          rnode(rpxmax1_, saveigrid),rnode(rpxmin2_, saveigrid),&
          rnode(rpxmax2_, saveigrid)
      write(unitterm,*)" dtdust =", dtdust(:)
      write(unitterm,*)"on processor:", mype
      write(unitterm,*)"-------------------------------------"
    endif

  end subroutine dust_get_dt

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                          :: vdust(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                                   :: n

    do n = 1, dust_n_species
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = get_vdust(w, ixImin1,ixImin2,&
         ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, n)

      if (present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = min(cmin(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2), vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2), vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2), abs(vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
      end if
    end do
  end subroutine dust_get_cmax

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax_prim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                          :: vdust(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                                   :: n

    do n = 1, dust_n_species
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = get_vdust_prim(w, ixImin1,&
         ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, n)

      if (present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = min(cmin(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2), vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2), vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2), abs(vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
      end if
    end do
  end subroutine dust_get_cmax_prim

end module mod_dust
