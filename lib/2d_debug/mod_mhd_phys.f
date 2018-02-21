!> Magneto-hydrodynamics module
module mod_mhd_phys
  use mod_global_parameters, only: std_len
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: mhd_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: mhd_radiative_cooling = .false.

  !> Whether viscosity is added
  logical, public, protected              :: mhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: mhd_gravity = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_Hall = .false.

  !> Whether particles module is added
  logical, public, protected              :: mhd_particles = .false.

  !> Whether magnetofriction is added
  logical, public, protected              :: mhd_magnetofriction = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: mhd_glm = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b 
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0

  !> MHD fourth order
  logical, public, protected              :: mhd_4th_order = .false.

  !> Number of tracer species
  integer, public, protected              :: mhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> The adiabatic index
  double precision, public                :: mhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: mhd_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: mhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: mhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public                :: mhd_etah = 0.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len) :: typedivbfix  = 'linde'

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*2)=.true.

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  ! Public methods
  public :: mhd_phys_init
  public :: mhd_kin_en
  public :: mhd_get_pthermal
  public :: mhd_get_v
  public :: mhd_to_conserved
  public :: mhd_to_primitive
  public :: mhd_get_csound2
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_gamma, mhd_adiab,mhd_eta,&
        mhd_eta_hyper, mhd_etah, mhd_glm_alpha, mhd_magnetofriction,&
       mhd_thermal_conduction, mhd_radiative_cooling, mhd_Hall, mhd_gravity,&
       mhd_viscosity, mhd_4th_order, typedivbfix, divbdiff, typedivbdiff,&
        compactres, divbwave, He_abundance, SI_unit, B0field,&
        B0field_forcefree, Bdip, Bquad, Boct, Busr,mhd_particles,&
        boundary_divbfix

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine mhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = mhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mhd_write_info

  subroutine mhd_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout)    :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim),  wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    integer, intent(in)                :: idim
    integer                            :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, iw
    double precision                   :: inv_volume(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call mpistop("to do")
    ! ! shifted indexes
    ! hxO^L=ixO^L-kr(idim,^D);
    ! ! all the indexes
    ! kxCmin^D=hxOmin^D;
    ! kxCmax^D=ixOmax^D;
    !
    ! inv_volume = 1.0d0/block%dvolume(ixO^S)
    !
    ! select case(typeaxial)
    ! case ("cylindrical")
    !   do iw=1,nwflux
    !     if (idim==r_ .and. iw==iw_mom(phi_)) then
    !       fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/x(ixO^S,r_))
    !     else
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !             inv_volume
    !     endif
    !   enddo
    ! case ("spherical")
    !   do iw=1,nwflux
    !     if     (idim==r_ .and. (iw==iw_mom(2) .or. iw==iw_mom(phi_))) then
    !       fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/x(ixO^S,r_))
    !     elseif (idim==2  .and. iw==iw_mom(phi_)) then
    !       fC(kxC^S,iw,idim)=fC(kxC^S,iw,idim)*dsin(x(kxC^S,2)+half*block%dx(kxC^S,2)) ! (x(4,3,1)-x(3,3,1)))
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !            (inv_volume/dsin(x(ixO^S,2)))
    !     else
    !       wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !             inv_volume
    !     endif
    !   enddo
    !
    !   ! if (idim==r_) then
    !   !   fC(kxC^S,iw_mom(phi_),idim)= fC(kxC^S,iw_mom(phi_),idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !   !   fC(kxC^S,iw_mom(phi_),idim)= fC(kxC^S,iw_mom(phi_),idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,r_))
    !   !   wnew(ixO^S,iw_mom(phi_))=wnew(ixO^S,iw_mom(phi_)) + (fC(ixO^S,iw_mom(phi_),idim)-fC(hxO^S,iw_mom(phi_),idim)) * &
    !   !        (inv_volume/x(ixO^S,r_))
    !   !
    !   ! elseif (idim==2) then
    !   !   fC(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,iw,idim)=fC(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,iw,idim)*dsin(x(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,2)+half*block%dx(hxOmin1:hxOmax1,hxOmin2,hxOmin3:hxOmax3,2)) ! (x(4,3,1)-x(3,3,1)))
    !   !   fC(ixO^S,iw,idim)=fC(ixO^S,iw,idim)*dsin(x(ixO^S,2)+half*block%dx(ixO^S,2)) ! (x(4,3,1)-x(3,3,1)))
    !   !   wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
    !   !        (inv_volume/dsin(x(ixO^S,2)))
    !   ! endif
    !
    ! end select

  end subroutine mhd_angmomfix

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_magnetofriction, only: magnetofriction_init
    use mod_physics

    integer :: itr, idir

    call mhd_read_params(par_files)

    physics_type = "mhd"
    phys_energy=mhd_energy
    use_particles=mhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
      type_divb=0
    case ('glm1')
      mhd_glm=.true.
      need_global_cmax=.true.
      type_divb=1
    case ('glm2')
      mhd_glm=.true.
      need_global_cmax=.true.
      type_divb=2
    case ('glm3')
      mhd_glm=.true.
      need_global_cmax=.true.
      type_divb=3
    case ('powel', 'powell')
      type_divb=4
    case ('janhunen')
      type_divb=5
    case ('linde')
      type_divb=6
    case ('lindejanhunen')
      type_divb=7
    case ('lindepowel')
      type_divb=8
    case default
      call mpistop('Unknown divB fix')
    end select

    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (mhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    if (mhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    allocate(tracer(mhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, mhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    do idir=1,ndir
      if(ndim>1) flux_type(idir,mag(idir))=flux_tvdlf
    end do
    if(mhd_glm .and. ndim>1) flux_type(:,psi_)=flux_tvdlf

    phys_get_dt              => mhd_get_dt
    phys_get_cmax            => mhd_get_cmax
    phys_get_cbounds         => mhd_get_cbounds
    phys_get_flux            => mhd_get_flux
    phys_add_source_geom     => mhd_add_source_geom
    phys_add_source          => mhd_add_source
    phys_to_conserved        => mhd_to_conserved
    phys_to_primitive        => mhd_to_primitive
    phys_check_params        => mhd_check_params
    phys_check_w             => mhd_check_w
    phys_get_pthermal        => mhd_get_pthermal
    phys_boundary_adjust     => mhd_boundary_adjust
    phys_write_info          => mhd_write_info
    phys_angmomfix           => mhd_angmomfix
    phys_handle_small_values => mhd_handle_small_values

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb <6) phys_req_diagonal = .true.

    ! derive units from basic units
    call mhd_physical_units()

    if(.not. mhd_energy .and. mhd_thermal_conduction) then
      call mpistop("thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_radiative_cooling) then
      call mpistop("radiative cooling needs mhd_energy=T")
    end if

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
      phys_req_diagonal = .true.
      call thermal_conduction_init(mhd_gamma)
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init(mhd_gamma,He_abundance)
    end if

    ! Initialize viscosity module
    if(mhd_viscosity) then
      call viscosity_init()
    end if

    ! Initialize gravity module
    if(mhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(mhd_particles) then
      call particles_init()
      phys_req_diagonal = .true.
    end if

    ! initialize magnetofriction module
    if (mhd_magnetofriction) then
      phys_req_diagonal = .true.
      call magnetofriction_init()
    end if

    if (mhd_glm) then
       ! Solve the Riemann problem for the linear 2x2 system for normal
       ! B-field and GLM_Psi according to Dedner 2002:
       phys_modify_wLR => glmSolve
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if (mhd_hall) then
       if (mhd_4th_order) then
          phys_wider_stencil = 2
       else
          phys_wider_stencil = 1
       end if
    end if

  end subroutine mhd_phys_init

  subroutine mhd_check_params
    use mod_global_parameters

    ! after user parameter setting
    gamma_1=mhd_gamma-1.d0
    inv_gamma_1=1.d0/gamma_1

    if (.not. mhd_energy) then
       if (mhd_gamma <= 0.0d0) call mpistop ("Error: mhd_gamma <= 0")
       if (mhd_adiab < 0.0d0) call mpistop ("Error: mhd_adiab < 0")
       small_pressure = mhd_adiab*small_density**mhd_gamma
    else
       if (mhd_gamma <= 0.0d0 .or. mhd_gamma == 1.0d0) call mpistop &
          ("Error: mhd_gamma <= 0 or mhd_gamma == 1")
       small_e = small_pressure * inv_gamma_1
    end if

  end subroutine mhd_check_params

  subroutine mhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
    end if
    if(unit_velocity==0) then
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=(2.d0+3.d0*He_abundance)&
         *unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    end if

  end subroutine mhd_physical_units

  subroutine mhd_check_w(primitive,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    integer, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)

    flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_

    if (mhd_energy) then
       if (block%e_is_internal) then
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              e_) < small_pressure*inv_gamma_1) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) = e_
       else
         if (primitive)then
           where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               e_) < small_pressure) flag(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) = e_
         else
        ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gamma_1*(w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,e_) - mhd_kin_en(w,ixImin1,ixImin2,ixImax1,&
              ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2)-mhd_mag_en(w,ixImin1,&
              ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2))
           where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < small_pressure) &
              flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = e_
         end if
       end if
    end if
  end subroutine mhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    integer                         :: idir, itr

    if (mhd_energy) then
       ! Calculate total energy from pressure, kinetic and magnetic energy
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          p_)*inv_gamma_1
      if(.not.block%e_is_internal) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_) + 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
          dim=ndim+1) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          rho_) + mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))
    end do

    if (check_small_values) call mhd_handle_small_values(.false., w, x,&
        ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'mhd_to_conserved')
  end subroutine mhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision                :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    integer                         :: itr, idir

    inv_rho = 1.0d0 / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)

    if (mhd_energy) then
      if(block%e_is_internal) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, p_) = gamma_1*w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, e_)
      else
       ! Calculate pressure = (gamma-1) * (e-0.5*(2ek+2eb))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, p_) = gamma_1*(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, e_) - mhd_kin_en(w, ixImin1,ixImin2,ixImax1,&
           ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, inv_rho) - mhd_mag_en(w,&
            ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2))
      end if
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(idir))*inv_rho
    end do

    if (check_small_values) call mhd_handle_small_values(.true., w, x, ixImin1,&
       ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'mhd_to_primitive')
  end subroutine mhd_to_primitive

  subroutine mhd_handle_small_values(primitive, w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: smallone
    integer :: idir, flag(ixImin1:ixImax1,ixImin2:ixImax2)

    if (small_values_method == "ignore") return

    call mhd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, w, flag)

    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,rho_) = small_density

          do idir = 1, ndir
             where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = 0.0d0
          end do

          if (mhd_energy) then
             if(primitive) then
               smallone = small_pressure
             else
               smallone = small_e
             end if
             where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) &
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = smallone
          end if
       case ("average")
          call small_values_average(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2, w, x, flag)
       case default
          call small_values_error(w, x, ixImin1,ixImin2,ixImax1,ixImax2,&
              ixOmin1,ixOmin2,ixOmax1,ixOmax2, flag, subname)
       end select
    end if
  end subroutine mhd_handle_small_values

  !> Convert energy to entropy
  subroutine e_to_rhos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision,intent(inout)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    if (mhd_energy) then
      if(.not.block%e_is_internal) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) - mhd_kin_en(w, ixImin1,&
         ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
         ixOmax2) - mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = gamma_1* w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, rho_)**(1.0d0 - mhd_gamma) * w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, e_)
    else
      call mpistop("e_to_rhos can not be used without energy equation!")
    end if
  end subroutine e_to_rhos

  !> Convert entropy to energy
  subroutine rhos_to_e(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    if (mhd_energy) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)**gamma_1 * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) * inv_gamma_1
       if(.not.block%e_is_internal) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) + mhd_kin_en(w, ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
          ixOmax2) + mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

  !> Calculate v vector
  subroutine mhd_get_v(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,ndir)

    integer :: idir

    do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mom(idir)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    end do

  end subroutine mhd_get_v

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)

    call mhd_get_csound(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idim,cmax)

    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(idim))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_))+cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

  end subroutine mhd_get_cmax

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine mhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: umean,&
        dmean, csoundL, csoundR, tmp1,tmp2,tmp3

    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=1.d0/(sqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))+sqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      call mhd_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,csoundR)
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundL(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2+tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+0.5d0*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)))**2
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(dmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
    else
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=wmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)
      call mhd_get_csound(wmean,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
    end if

  end subroutine mhd_get_cbounds

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2),&
        AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2), b2(ixImin1:ixImax1,&
       ixImin2:ixImax2), kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    inv_rho=1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    call mhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,csound)
    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)        = mhd_mag_en_all(w,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = b2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cfast2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2-4.0d0*csound(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * mhd_mag_i_all(w,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(AvMinCs2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    if (.not. MHD_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))),&
           mhd_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound_prim(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2),&
        AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2), b2(ixImin1:ixImax1,&
       ixImin2:ixImax2), kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    inv_rho=1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    if(mhd_energy) then
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,p_)*inv_rho
    else
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*mhd_adiab*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**gamma_1
    end if
    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)        = mhd_mag_en_all(w,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = b2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cfast2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2-4.0d0*csound(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * mhd_mag_i_all(w,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(AvMinCs2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    if (.not. MHD_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))),&
           mhd_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound_prim

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,pth)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2)

    if(mhd_energy) then
      if(block%e_is_internal) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gamma_1*w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)
      else
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gamma_1*(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)- mhd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2)- mhd_mag_en(w,ixImin1,ixImin2,&
           ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2))
      end if
    else
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_adiab*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)**mhd_gamma
    end if
  end subroutine mhd_get_pthermal

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine mhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    if(mhd_energy) then
      call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*csound2(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)
    else
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*mhd_adiab*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**gamma_1
    end if
  end subroutine mhd_get_csound2

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine mhd_get_p_total(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,p)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: p(ixImin1:ixImax1,ixImin2:ixImax2)

    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,p)

    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mag(:))**2, dim=ndim+1)

  end subroutine mhd_get_p_total

  !> Calculate fluxes within ixO^L.
  subroutine mhd_get_flux(wC,w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,f)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,nwflux)

    double precision             :: ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, allocatable:: vHall(:,:,:)
    integer                      :: idirmin, iw, idir

    if (mhd_Hall) then
      allocate(vHall(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir))
      call mhd_getv_Hall(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,vHall)
    end if

    if(B0field) tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum(block%B0(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,idim)*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(:)),dim=ndim+1)

    if(mhd_energy) then
      ptotal=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)+0.5d0*sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mag(:))**2, dim=ndim+1)
    else
      ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_adiab*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)**mhd_gamma+0.5d0*sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mag(:))**2, dim=ndim+1)
    end if

    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tracer(iw))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=ptotal(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))
        if(B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))=f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))= -w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))
      end if
      if (B0field) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=f(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idir))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           idim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,idim)
      end if
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=f(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if (mhd_energy) then
       if (block%e_is_internal) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)
          if (mhd_Hall) then
             call mpistop("solve pthermal not designed for Hall MHD")
          endif
       else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mom(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             e_) + ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2))- w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mag(idim))*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(:)),dim=ndim+1)

          if (B0field) then
             f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = f(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,e_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                mom(idim)) * tmp(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2) - sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                mom(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),&
                dim=ndim+1) * block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
                idim)
          end if

          if (mhd_Hall) then
          ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
             if (mhd_etah>zero) then
                f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = f(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,e_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   idim) * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2,&
                   dim=ndim+1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   mag(idim)) * sum(vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   :)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),dim=ndim+1)
                if (B0field) then
                   f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = f(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,e_) + vHall(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,idim) * tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2) - sum(vHall(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,:)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      mag(:)),dim=ndim+1) * block%B0(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,idim,idim)
                end if
             end if
          end if
       end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))

        if (B0field) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=f(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mag(idir))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idim))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             idim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idir))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,idim)
        end if

        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (mhd_etah>zero) then
            if (B0field) then
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = f(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,mag(idir)) - vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,idir)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(idim))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
                 idim)) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idim)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
                 idim))
            else
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = f(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,mag(idir)) - vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,idir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(idim)) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))
            end if
          end if
        end if

      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)  = &
         cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim))
    end if

  end subroutine mhd_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mhd_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if (mhd_energy .and. block%e_is_internal) then
        active = .true.
        call internal_energy_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      endif

      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        call add_source_B0split(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(mhd_eta)>smalldouble)then
        if (.not.slab) call mpistop("no resistivity in non-slab geometry")
        active = .true.
        if (compactres)then
          call add_source_res1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        else
          call add_source_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        end if
      end if

      if (mhd_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      end if
      
      if(istep==nstep) then
        ! Sources related to div B
        select case (type_divb)
        case (0)
          ! Do nothing
        case (1)
          active = .true.
          call add_source_glm1(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case (2)
          active = .true.
          call add_source_glm2(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case (3)
          active = .true.
          call add_source_glm3(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case (4)
          active = .true.
          call add_source_powel(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case (5)
          active = .true.
          call add_source_janhunen(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case (6)
          active = .true.
          call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case (7)
          active = .true.
          call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
          call add_source_janhunen(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case (8)
          active = .true.
          call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
          call add_source_powel(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pw(saveigrid)%wold,w,x)
        case default
          call mpistop('Unknown divB fix')
        end select
      end if
     
    end if

    if(mhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    end if

    if(mhd_viscosity) then
      call viscosity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,mhd_energy,qsourcesplit,active)
    end if

    if(mhd_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,mhd_energy,qsourcesplit,active)
    end if

  end subroutine mhd_add_source

  subroutine internal_energy_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),divv(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,v)
    call divvector(v,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divv)
    call mhd_get_pthermal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,pth)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)-qdt*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine internal_energy_add_source

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: a(ixImin1:ixImax1,ixImin2:ixImax2,3),&
        b(ixImin1:ixImax1,ixImin2:ixImax2,3), axb(ixImin1:ixImax1,&
       ixImin2:ixImax2,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir)=block%B0(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:ndir,0)

      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=block%J0(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,idir)
      end do
      call cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,a,b,axb)
      axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=axb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,:)*qdt
      ! add J0xB0 source term in momentum equations
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1:ndir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(1:ndir))+axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:ndir)
    end if

    if(mhd_energy) then
      if(.not.block%e_is_internal) then
        a=0.d0
        ! for free-free field -(vxB0) dot J0 =0
        b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mag(:))
        ! store full magnetic field B0+B1 in b
        if(.not.B0field_forcefree) b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           :)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)+block%B0(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,:,0)
        ! store velocity in a
        do idir=1,ndir
          a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mom(idir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             rho_)
        end do
        call cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,a,b,axb)
        axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=axb(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,:)*qdt
        ! add -(vxB) dot J0 source term in energy equation
        do idir=7-2*ndir,3
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_)-axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idir)*block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
        end do
      end if
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_B0')

  end subroutine add_source_B0split

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idir,jdir,kdir,idirmin,idim,&
       jxOmin1,jxOmin2,jxOmax1,jxOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,ix
    integer :: lxOmin1,lxOmin2,lxOmax1,lxOmax2, kxOmin1,kxOmin2,kxOmax1,&
       kxOmax2

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3),&
       eta(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: gradeta(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (mhd_4th_order) then
      ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmax1=ixOmax1+2;ixAmax2=ixOmax2+2;
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2) call mpistop&
       ("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=mhd_eta
       gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixAmin1,ixAmin2,ixAmax1,ixAmax2,idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,idim,tmp)
          gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)=tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
       end do
    end if

    if(B0field) then
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=wCT(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndir,0)
    else
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=wCT(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))
    end if

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (mhd_4th_order) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=Bf(ixImin1:ixImax1,&
            ixImin2:ixImax2,idir)
         do idim=1,ndim
            lxOmin1=ixOmin1+2*kr(idim,1);lxOmin2=ixOmin2+2*kr(idim,2)
            lxOmax1=ixOmax1+2*kr(idim,1);lxOmax2=ixOmax2+2*kr(idim,2);
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmax1=ixOmax1+kr(idim,1);jxOmax2=ixOmax2+kr(idim,2);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
            kxOmin1=ixOmin1-2*kr(idim,1);kxOmin2=ixOmin2-2*kr(idim,2)
            kxOmax1=ixOmax1-2*kr(idim,1);kxOmax2=ixOmax2-2*kr(idim,2);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+(-tmp2(lxOmin1:lxOmax1,&
               lxOmin2:lxOmax2)+16.0d0*tmp2(jxOmin1:jxOmax1,&
               jxOmin2:jxOmax2)-30.0d0*tmp2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+16.0d0*tmp2(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2)-tmp2(kxOmin1:kxOmax1,&
               kxOmin2:kxOmax2)) /(12.0d0 * dxlevel(idim)**2)
         end do
       else
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=Bf(ixImin1:ixImax1,&
            ixImin2:ixImax2,idir)
         do idim=1,ndim
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmax1=ixOmax1+kr(idim,1);jxOmax2=ixOmax2+kr(idim,2);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+(tmp2(jxOmin1:jxOmax1,&
               jxOmin2:jxOmax2)-2.0d0*tmp2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+tmp2(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2))/dxlevel(idim)**2
         end do
       end if

       ! Multiply by eta
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (mhd_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2)-gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      jdir)*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir)
                else
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2)+gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      jdir)*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir)
                end if
             end if
          end do; end do
       end if

       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(idir))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if (mhd_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_)+qdt*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
       end if

    end do ! idir

    if (mhd_energy) then
       ! de/dt+=eta*J**2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
       do idir=idirmin,3
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)+current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)**2
       end do
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_)+qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_res1')

  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idir,jdir,kdir,idirmin,iw,idim,&
       idirmin1

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3),&
       eta(ixImin1:ixImax1,ixImin2:ixImax2),curlj(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3)
    double precision :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,1:3)

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmax1=ixOmax1+2;ixAmax2=ixOmax2+2;

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2) call mpistop&
       ("Error in add_source_res2: Non-conforming input limits")

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=mhd_eta
    else
       call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixAmin1,ixAmin2,ixAmax1,ixAmax2,idirmin,x,current,eta)
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,jdir)=current(ixAmin1:ixAmax1,&
          ixAmin2:ixAmax2,jdir)*eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2)
    end do
    call curlvector(tmpvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,curlj,idirmin1,1,3)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)
    end do

    if(mhd_energy) then
      ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
      ! de1/dt= eta J^2 - B1 dot curl(eta J)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+qdt*(sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)**2,&
         dim=ndim+1)*eta(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(1:ndir))*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir),&
         dim=ndim+1))
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_res2')

  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    !.. local ..
    double precision                :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndir:3)
    double precision                :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3),tmpvec2(ixImin1:ixImax1,ixImin2:ixImax2,1:3),tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2),ehyper(ixImin1:ixImax1,ixImin2:ixImax2,1:3)
    integer                         :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idir,&
       jdir,kdir,idirmin,idirmin1

    ixAmin1=ixOmin1-3;ixAmin2=ixOmin2-3;ixAmax1=ixOmax1+3;ixAmax2=ixOmax2+3;
    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2) call mpistop&
       ("Error in add_source_hyperres: Non-conforming input limits")

    call get_current(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,idirmin,current)
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,jdir)=current(ixAmin1:ixAmax1,&
          ixAmin2:ixAmax2,jdir)
    end do

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmax1=ixOmax1+2;ixAmax2=ixOmax2+2;
    call curlvector(tmpvec,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,tmpvec2,idirmin1,1,3)

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    call curlvector(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,tmpvec,idirmin1,1,3)
    ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir) = - tmpvec(ixAmin1:ixAmax1,&
       ixAmin2:ixAmax2,1:ndir)*mhd_eta_hyper

    ixAmin1=ixOmin1;ixAmin2=ixOmin2;ixAmax1=ixOmax1;ixAmax2=ixOmax2;
    tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    call curlvector(ehyper,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-tmpvec2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*qdt
    end do

    if (mhd_energy) then
      ! de/dt= +div(B x Ehyper)
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
      tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idir) = tmpvec(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,idir)+ lvc(idir,jdir,kdir)*wCT(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,mag(jdir))*ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           kdir)
      end do; end do; end do
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
      call divvector(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,tmp)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*qdt
    end if

    if (check_small_values)  call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision:: divb(ixImin1:ixImax1,ixImin2:ixImax2)
    integer          :: idim,idir
    double precision :: gradPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         abs(mhd_glm_alpha)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
             ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
             gradPsi)
       case("limited")
          call gradientS(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
             ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
             gradPsi)
       end select
       if (mhd_energy .and. .not.block%e_is_internal) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,e_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(idim))*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*mhd_mag_i_all(w,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_glm1')

  end subroutine add_source_glm1

  subroutine add_source_glm2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 38_
    ! giving the non conservative EGLM-MHD scheme.
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qdt,   wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision:: divb(ixImin1:ixImax1,ixImin2:ixImax2),v(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndir)
    integer          :: idim,idir
    double precision :: gradPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    ! calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,v)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         abs(mhd_glm_alpha)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
         call gradient(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
            ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
            gradPsi)
       case("limited")
         call gradientS(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
            ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
            gradPsi)
       end select

      ! Psi=Psi - qdt (v . grad(Psi))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,psi_)-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idim)*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

      if (mhd_energy .and. .not.block%e_is_internal) then
      ! e  = e  - qdt (b . grad(Psi))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
    end do

    if (mhd_energy .and. .not.block%e_is_internal) then
    ! e = e - qdt (v . b) * div b
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) - qdt * divb(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          :)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),dim=ndim+1)
    end if
    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*mhd_mag_i_all(w,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_glm2')

  end subroutine add_source_glm2

  subroutine add_source_glm3(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation (1a), (1b), (4), (1d), 19
    ! conservative hyperbolic mixed GLM-MHD with no additional source terms.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         abs(mhd_glm_alpha)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    end if

  end subroutine add_source_glm3

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2),&
       v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,v)

    if (mhd_energy .and. .not.block%e_is_internal) then
      ! e = e - qdt (v . b) * div b
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)-qdt*sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)*wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(:)),dim=ndim+1)*divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*mhd_mag_i_all(w,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: idim, idir, ixpmin1,ixpmin2,ixpmax1,ixpmax2, i1,i2, iside
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2),&
       graddivb(ixImin1:ixImax1,ixImin2:ixImax2)


    ! Calculate div B
    ixpmin1=ixOmin1-1;ixpmin2=ixOmin2-1;ixpmax1=ixOmax1+1;ixpmax2=ixOmax2+1;
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,ixpmax1,&
       ixpmax2,divb)
    ! for AMR stability, retreat one cell layer from the boarders of level jump
    ixpmin1=ixOmin1;ixpmin2=ixOmin2;ixpmax1=ixOmax1;ixpmax2=ixOmax2;
    do idim=1,ndim
      select case(idim)
       case(1)
          do iside=1,2
            i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
            if (leveljump(i1,i2)) then
              if (iside==1) then
                ixpmin1=ixOmin1-i1
              else
                ixpmax1=ixOmax1-i1
              end if
            end if
          end do
       
       case(2)
          do iside=1,2
            i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
            if (leveljump(i1,i2)) then
              if (iside==1) then
                ixpmin2=ixOmin2-i2
              else
                ixpmax2=ixOmax2-i2
              end if
            end if
          end do
       
      end select
    end do

    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
            ixpmax1,ixpmax2,idim,graddivb)
       case("limited")
         call gradientS(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
            ixpmax1,ixpmax2,idim,graddivb)
       end select

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab) then
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2)*divbdiff/(1.0d0/dxlevel(1)**2+&
             1.0d0/dxlevel(2)**2)
       else
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2)*divbdiff /(1.0d0/block%dx(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,1)**2+1.0d0/block%dx(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,2)**2)
       end if

       w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,mag(idim))=w(ixpmin1:ixpmax1,&
          ixpmin2:ixpmax2,mag(idim))+graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)

       if (mhd_energy .and. typedivbdiff=='all' .and. &
          .not.block%e_is_internal) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,e_)=w(ixpmin1:ixpmax1,&
            ixpmin2:ixpmax2,e_)+wCT(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
            mag(idim))*graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)
       end if
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_linde')

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision                   :: divb(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    double precision                   :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)

    bvec(ixImin1:ixImax1,ixImin2:ixImax2,:)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mag(:))

    select case(typediv)
    case("central")
      call divvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,divb)
    case("limited")
      call divvectorS(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,divb)
    end select
  end subroutine get_divb

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision                   :: divb(ixImin1:ixImax1,&
       ixImin2:ixImax2), dsurface(ixImin1:ixImax1,ixImin2:ixImax2)

    integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idims

    call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)
    if(slab) then
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0.5d0*abs(divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/sqrt(mhd_mag_en_all(w,ixImin1,ixImin2,ixImax1,&
         ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2))/sum(1.d0/dxlevel(:))
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;
      ixAmax1=ixOmax1-1;ixAmax2=ixOmax2-1;
      dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= &
         sum(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:),dim=ndim+1)
      do idims=1,ndim
        ixAmin1=ixOmin1-kr(idims,1);ixAmin2=ixOmin2-kr(idims,2)
        ixAmax1=ixOmax1-kr(idims,1);ixAmax2=ixOmax2-kr(idims,2);
        dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsurface(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+block%surfaceC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/sqrt(mhd_mag_en_all(w,ixImin1,ixImin2,ixImax1,&
         ixImax2,ixOmin1,ixOmin2,ixOmax1,&
         ixOmax2))*block%dvolume(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idirmin,current)
    use mod_global_parameters

    integer :: idirmin0
    integer :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, idirmin, ixImin1,ixImin2,&
       ixImax1,ixImax2
    double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: idir

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3),&
       bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    idirmin0 = 7-2*ndir

    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mag(1:ndir))

    call curlvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin0:3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin0:3)+block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idirmin0:3)

  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mhd_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx1,dx2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2)

    dtnew = bigdouble

    dxarr(1)=dx1;dxarr(2)=dx2;
    if (mhd_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mhd_eta
    else if (mhd_eta<zero)then
       call get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,idirmin,current)
       call usr_special_resistivity(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,idirmin,x,current,eta)
       dtnew=bigdouble
       do idim=1,ndim
          dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/dxarr(idim)**2)))
       end do
    end if

    if (mhd_eta_hyper>zero)then
       dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mhd_eta_hyper,dtnew)
    end if

    if(mhd_radiative_cooling) then
      call cooling_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

    if(mhd_viscosity) then
      call viscosity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

    if(mhd_gravity) then
      call gravity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

  end subroutine mhd_get_dt

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmin2,h1xmax1,h1xmax2, h2xmin1,&
       h2xmin2,h2xmax1,h2xmax2
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (typeaxial)
    case ('cylindrical')
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
       call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,tmp)
       if(phi_>0) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mr_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,bphi_)**2+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mphi_)**2/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mphi_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*(-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mr_)/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            rho_) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,br_))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bphi_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,bphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bphi_)*wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mphi_)) /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
       else
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mr_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
       if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,br_)=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,br_)+qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          psi_)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
    case ('spherical')
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmax1=ixOmax1-kr(1,1);h1xmax2=ixOmax2-kr(1,2)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmax1=ixOmax1-kr(2,1);h2xmax2=ixOmax2-kr(2,2);
       call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,tmp1)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp1(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       if(B0field) then
         tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum(block%B0(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,:,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),&
            dim=ndim+1)
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
       ! m1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)-block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,&
          1))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(idir))**2/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))**2
           if(B0field) tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              idir,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       ! b1
       if(mhd_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mag(1))+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
       end if

       
       ! m2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp1(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       if(B0field) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
       ! This will make hydrostatic p=const an exact solution
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(2))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-(wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mom(2))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2)))
       if (B0field) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
             0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(2)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,0)
       end if
       if(ndir==3) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(3))**2/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3))**2)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               3,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mag(3))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
         end if
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(2))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       ! b2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(1)))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
       if(B0field) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
            0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
            0))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
       end if
       if(mhd_glm) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,psi_)
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(2))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
      

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-(wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(1))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(1)))  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(3))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3))) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
           if (B0field) then
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
                 0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
                 0)  +(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
                 0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
                 0)) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
           end if
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(3))=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(3))+qdt*tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(1)))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            rho_)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            2)) /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            rho_)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
               0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(3))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
               0))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(3))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
               0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
               0))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               2)) /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mag(3))+qdt*tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       end if
    end select
  end subroutine mhd_add_source_geom

  !> Compute 2 times total magnetic energy
  function mhd_mag_en_all(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (B0field) then
      mge = sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,block%iw0))**2,&
          dim=ndim+1)
    else
      mge = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2, dim=ndim+1)
    end if
  end function mhd_mag_en_all

  !> Compute full magnetic field by direction
  function mhd_mag_i_all(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (B0field) then
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,block%iw0)
    else
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(idir))
    end if
  end function mhd_mag_i_all

  !> Compute evolving magnetic energy
  function mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    mge = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2,&
        dim=ndim+1)
  end function mhd_mag_en

  !> compute kinetic energy
  function mhd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    if (present(inv_rho)) then
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) * inv_rho
    else
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
    end if
  end function mhd_kin_en

  subroutine mhd_getv_Hall(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3)

    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3)

    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,idirmin,current)
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:3) = zero
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin:3) = - mhd_etah*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin:3)
    do idir = idirmin, 3
       vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = vHall(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,idir)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    end do

  end subroutine mhd_getv_Hall

  subroutine mhd_getdt_Hall(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,dx1,dx2,dthall)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: dthall
    !.. local ..
    double precision :: dxarr(ndim)
    double precision :: bmag(ixImin1:ixImax1,ixImin2:ixImax2)

    dthall=bigdouble

    ! because we have that in cmax now:
    return

    dxarr(1)=dx1;dxarr(2)=dx2;

    if (.not. B0field) then
       bmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(sum(w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(:))**2, dim=ndim+1))
       bmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(sum((w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(:)) + block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1:ndir,block%iw0))**2))
    end if

    dthall=dtdiffpar*minval(dxarr(1:ndim))**&
       2.0d0/(mhd_etah*maxval(bmag(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)))
  end subroutine mhd_getdt_Hall

  !> This implements eq. (42) in Dedner et al. 2002 JcP 175
  !> Gives the Riemann solution on the interface
  !> for the normal B component and Psi in the GLM-MHD system.
  !> 23/04/2013 Oliver Porth
  subroutine glmSolve(wLC,wRC,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idir)
    use mod_global_parameters
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    double precision                :: dB(ixImin1:ixImax1,ixImin2:ixImax2),&
        dPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = wRC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(idir)) - wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mag(idir))
    dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wRC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_) - wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)

    wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))   = 0.5d0 * &
       (wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) + wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(idir))) - 0.5d0/cmax_global * dPsi(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)       = 0.5d0 * &
       (wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) + wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_)) - 0.5d0*cmax_global * dB(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(idir))
    wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = wLC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,psi_)

  end subroutine glmSolve

  subroutine mhd_boundary_adjust
    use mod_global_parameters
    integer :: iB, idim, iside, iigrid, igrid
    integer :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2, i1,i2

    ixGmin1=ixGlo1;ixGmin2=ixGlo2;ixGmax1=ixGhi1;ixGmax2=ixGhi2;
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>pw(igrid)
        dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
        do idim=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
           do iside=1,2
              i1=kr(1,idim)*(2*iside-3);i2=kr(2,idim)*(2*iside-3);
              if (neighbor_type(i1,i2,igrid)/=1) cycle
              iB=(idim-1)*2+iside
              if(.not.boundary_divbfix(iB)) cycle
              if(any(typeboundary(:,iB)=="special")) then
                select case (idim)
                case (1)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin1=ixGmax1+1-nghostcells;ixOmin2=ixGmin2;
                      ixOmax1=ixGmax1;ixOmax2=ixGmax2;
                   else
                      ! minimal boundary
                      ixOmin1=ixGmin1;ixOmin2=ixGmin2;
                      ixOmax1=ixGmin1-1+nghostcells;ixOmax2=ixGmax2;
                   end if 
                case (2)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin1=ixGmin1;ixOmin2=ixGmax2+1-nghostcells;
                      ixOmax1=ixGmax1;ixOmax2=ixGmax2;
                   else
                      ! minimal boundary
                      ixOmin1=ixGmin1;ixOmin2=ixGmin2;
                      ixOmax1=ixGmax1;ixOmax2=ixGmin2-1+nghostcells;
                   end if 
                end select
                call fixdivB_boundary(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,pw(igrid)%wb,pw(igrid)%x,iB)
              end if
           end do
        end do
     end do

  end subroutine mhd_boundary_adjust

  subroutine fixdivB_boundary(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,iB
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)

    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix1,ix2,ixFmin1,ixFmin2,ixFmax1,ixFmax2

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=w(ix1+1,ixFmin2:ixFmax2,&
              mag(1)) +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=( (w(ix1+1,ixFmin2:ixFmax2,&
              mag(1))+w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1,&
              ixFmin2:ixFmax2,1)+(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,&
              ixFmin2:ixFmax2,mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              2)-(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,&
              mag(2)))*block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,&
              2) )/block%surfaceC(ix1-1,ixFmin2:ixFmax2,1)-w(ix1,&
              ixFmin2:ixFmax2,mag(1))
         end do
       end if
      
       
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     case(2)
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=w(ix1-1,ixFmin2:ixFmax2,&
              mag(1)) -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=( (w(ix1-1,ixFmin2:ixFmax2,&
              mag(1))+w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1-1,&
              ixFmin2:ixFmax2,1)-(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,&
              ixFmin2:ixFmax2,mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              2)+(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,&
              mag(2)))*block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,&
              2) )/block%surfaceC(ix1,ixFmin2:ixFmax2,1)-w(ix1,ixFmin2:ixFmax2,&
              mag(1))
         end do
       end if
      
       
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     case(3)
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=w(ixFmin1:ixFmax1,ix2+1,&
              mag(2)) +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         enddo
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              2)+(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              1)-(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,&
              1) )/block%surfaceC(ixFmin1:ixFmax1,ix2-1,2)-w(ixFmin1:ixFmax1,&
              ix2,mag(2))
         end do
       end if
      
       
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     case(4)
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=w(ixFmin1:ixFmax1,ix2-1,&
              mag(2)) -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2-1,&
              2)-(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              1)+(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,&
              1) )/block%surfaceC(ixFmin1:ixFmax1,ix2,2)-w(ixFmin1:ixFmax1,ix2,&
              mag(2))
         end do
       end if
      
       
       if(mhd_energy.and..not.block%e_is_internal) call &
          mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x)
     
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

   end subroutine fixdivB_boundary

end module mod_mhd_phys
