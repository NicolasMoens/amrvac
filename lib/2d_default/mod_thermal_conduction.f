!> Thermal conduction for HD and MHD
!>
!> 10.07.2011 developed by Chun Xia and Rony Keppens
!> 01.09.2012 moved to modules folder by Oliver Porth
!> 13.10.2013 optimized further by Chun Xia
!> 12.03.2014 implemented RKL2 super timestepping scheme to reduce iterations 
!> and improve stability and accuracy up to second order in time by Chun Xia.
!> 23.08.2014 implemented saturation and perpendicular TC by Chun Xia
!> 12.01.2017 modulized by Chun Xia
!> 
!> PURPOSE: 
!> IN MHD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(KAPPA_i,j . GRAD_j T)
!> where KAPPA_i,j = tc_k_para b_i b_j + tc_k_perp (I - b_i b_j)
!> b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!> IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(tc_k_para . GRAD T)
!> USAGE:
!> 1. in mod_usr.t -> subroutine usr_init(), add 
!>        unit_length=<your length unit>
!>        unit_numberdensity=<your number density unit>
!>        unit_velocity=<your velocity unit>
!>        unit_temperature=<your temperature unit>
!>    before call (m)hd_activate()
!> 2. to switch on thermal conduction in the (m)hd_list of amrvac.par add:
!>    (m)hd_thermal_conduction=.true.
!> 3. in the tc_list of amrvac.par :
!>    tc_perpendicular=.true.  ! (default .false.) turn on thermal conduction perpendicular to magnetic field 
!>    tc_saturate=.false.  ! (default .true. ) turn off thermal conduction saturate effect
!>    tc_dtpar=0.9/0.45/0.3 ! stable time step coefficient for 1D/2D/3D, decrease it for more stable run
!>    tc_slope_limiter='MC' ! choose limiter for slope-limited anisotropic thermal conduction in MHD

module mod_thermal_conduction
  use mod_global_parameters, only: std_len
  implicit none
  !> Coefficient of thermal conductivity (parallel to magnetic field)
  double precision, public :: tc_k_para

  !> Coefficient of thermal conductivity perpendicular to magnetic field
  double precision, public :: tc_k_perp

  !> Time step of thermal conduction
  double precision :: dt_tc

  !> Number of sub-steps of supertime stepping
  integer, public :: s

  !> Index of the density (in the w array)
  integer, private :: rho_

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

  !> Indices of the magnetic field
  integer, allocatable, private, protected :: mag(:)

  !> The adiabatic index
  double precision, private :: tc_gamma

  !> The small_est allowed energy
  double precision, private :: small_e

  !> Time step coefficient
  double precision, private :: tc_dtpar=0.9d0

  !> Maximal substeps of TC within one fluid time step to limit fluid time step
  integer :: tc_ncycles=1000

  !> Calculate thermal conduction perpendicular to magnetic field (.true.) or not (.false.)
  logical, private :: tc_perpendicular=.false.

  !> Consider thermal conduction saturation effect (.true.) or not (.false.)
  logical, private :: tc_saturate=.true.

  !> Logical switch for prepare mpi datatype only once
  logical, private :: first=.true.

  !> Logical switch for test constant conductivity
  logical, private :: tc_constant=.false.

  !> Whether to conserve fluxes at the current partial step
  logical :: fix_conserve_at_step = .true.

  !> Name of slope limiter for transverse component of thermal flux 
  character(len=std_len), private  :: tc_slope_limiter

  procedure(thermal_conduction), pointer   :: phys_thermal_conduction => &
     null()
  procedure(get_heatconduct), pointer   :: phys_get_heatconduct => null()
  procedure(getdt_heatconduct), pointer :: phys_getdt_heatconduct => null()

  abstract interface
    subroutine thermal_conduction
    ! Meyer 2012 MNRAS 422,2102
      use mod_global_parameters
      use mod_ghostcells_update
    end subroutine thermal_conduction

    subroutine get_heatconduct(tmp,tmp1,tmp2,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x,qvec)
      use mod_global_parameters
      
      integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2
      double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
      !! tmp store the heat conduction energy changing rate
      double precision, intent(out) :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
         tmp1(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
         ixImin2:ixImax2)
      double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    end subroutine get_heatconduct

    subroutine getdt_heatconduct(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,&
       ixmin2,ixmax1,ixmax2,dtnew,dx1,dx2,x)
      use mod_global_parameters
      
      integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
         ixmax1,ixmax2
      double precision, intent(in) :: dx1,dx2, x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,1:ndim)
      ! note that depending on small_values_method=='error' etc, w values may change 
      ! through call to getpthermal
      double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
         1:nw), dtnew
    end subroutine getdt_heatconduct
  end interface

contains
  !> Read this module"s parameters from a file
  subroutine tc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_dtpar,&
        tc_slope_limiter, tc_k_para, tc_k_perp, tc_ncycles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, tc_list, end=111)
111    close(unitpar)
    end do

  end subroutine tc_params_read

  !> Initialize the module
  subroutine thermal_conduction_init(phys_gamma)
    use mod_global_parameters
    use mod_physics

    double precision, intent(in) :: phys_gamma

    tc_gamma=phys_gamma

    tc_dtpar=tc_dtpar/dble(ndim)

    tc_slope_limiter='MC'

    tc_k_para=0.d0

    tc_k_perp=0.d0

    call tc_params_read(par_files)
  
    phys_thermal_conduction => do_thermal_conduction
    if(physics_type=='hd') then
      phys_get_heatconduct   => hd_get_heatconduct
      phys_getdt_heatconduct => hd_getdt_heatconduct
    else if(physics_type=='mhd') then
      phys_get_heatconduct   => mhd_get_heatconduct
      phys_getdt_heatconduct => mhd_getdt_heatconduct
    end if

    allocate(mom(ndir), mag(ndir))
    rho_ = iw_rho
    e_ = iw_e
    mom(:) = iw_mom(:)
    mag(:) = iw_mag(:)

    small_e = small_pressure/(tc_gamma - 1.0d0)

    if(tc_k_para==0.d0 .and. tc_k_perp==0.d0) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        tc_k_para=8.d-12*unit_temperature**&
           3.5d0/unit_length/unit_density/unit_velocity**3 
        ! thermal conductivity perpendicular to magnetic field
        tc_k_perp=4.d-30*unit_numberdensity**2/unit_magneticfield**&
           2/unit_temperature**3*tc_k_para
      else
        ! Spitzer thermal conductivity with cgs units
        tc_k_para=8.d-7*unit_temperature**&
           3.5d0/unit_length/unit_density/unit_velocity**3 
        ! thermal conductivity perpendicular to magnetic field
        tc_k_perp=4.d-26*unit_numberdensity**2/unit_magneticfield**&
           2/unit_temperature**3*tc_k_para
      end if
    else
      tc_constant=.true.
    end if

  end subroutine thermal_conduction_init

  subroutine do_thermal_conduction
  ! Meyer 2012 MNRAS 422,2102
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    
    double precision :: omega1,cmu,cmut,cnu,cnut
    double precision, allocatable :: bj(:)
    integer:: iigrid, igrid,j
    logical :: evenstep

    ! point bc mpi datatype to partial type for thermalconduction
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1 
    ! create bc mpi datatype for ghostcells update
    if(first) then
      call create_bc_mpi_datatype(e_-1,1)
      first=.false.
    end if

    call init_comm_fix_conserve(1,ndim,1)
    fix_conserve_at_step = time_advance .and. levmax>levmin

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      if(.not. allocated(pw(igrid)%w2)) allocate(pw(igrid)%w2(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,1:nw))
      if(.not. allocated(pw(igrid)%w3)) allocate(pw(igrid)%w3(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,1:nw))
      pw(igrid)%w1=pw(igrid)%w
      pw(igrid)%w2=pw(igrid)%w
      pw(igrid)%w3=pw(igrid)%w
    end do
    
    allocate(bj(0:s))
    bj(0)=1.d0/3.d0
    bj(1)=bj(0)
    if(s>1) then
      omega1=4.d0/dble(s**2+s-2)
      cmut=omega1/3.d0
    else
      omega1=0.d0
      cmut=1.d0
    endif
    
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      block=>pw(igrid)
      typelimiter=type_limiter(node(plevel_,igrid))
      typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
      call evolve_step1(igrid,cmut,dt_tc,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,&
         ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w1,pw(igrid)%w,pw(igrid)%x,&
         pw(igrid)%w3)
      pw(igrid)%wb=>pw(igrid)%w1
    end do
    !$OMP END PARALLEL DO
    ! fix conservation of AMR grid by replacing flux from finer neighbors
    if (fix_conserve_at_step) then
      call recvflux(1,ndim)
      call sendflux(1,ndim)
      call fix_conserve(1,ndim,e_,1)
    end if
    bcphys=.false.
    call getbc(global_time,0.d0,e_-1,1)
    if(s==1) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,e_)=pw(igrid)%w1(ixGlo1:ixGhi1,&
           ixGlo2:ixGhi2,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
      bcphys=.true.
      deallocate(bj)
      return
    endif
    evenstep=.true.
    do j=2,s
      bj(j)=dble(j**2+j-2)/dble(2*j*(j+1))
      cmu=dble(2*j-1)/dble(j)*bj(j)/bj(j-1)
      cmut=omega1*cmu
      cnu=dble(1-j)/dble(j)*bj(j)/bj(j-2)
      cnut=(bj(j-1)-1.d0)*cmut
      if(evenstep) then
    !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          block=>pw(igrid)
          typelimiter=type_limiter(node(plevel_,igrid))
          typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
          dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
          call evolve_stepj(igrid,cmu,cmut,cnu,cnut,dt_tc,ixGlo1,ixGlo2,ixGhi1,&
             ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w1,pw(igrid)%w2,&
             pw(igrid)%w,pw(igrid)%x,pw(igrid)%w3)
          pw(igrid)%wb=>pw(igrid)%w2
        end do
    !$OMP END PARALLEL DO
        ! fix conservation of AMR grid by replacing flux from finer neighbors
        if (fix_conserve_at_step) then
          call recvflux(1,ndim)
          call sendflux(1,ndim)
          call fix_conserve(1,ndim,e_,1)
        end if
        call getbc(global_time,0.d0,e_-1,1)
        evenstep=.false.
      else
    !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          block=>pw(igrid)
          typelimiter=type_limiter(node(plevel_,igrid))
          typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
          dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
          call evolve_stepj(igrid,cmu,cmut,cnu,cnut,dt_tc,ixGlo1,ixGlo2,ixGhi1,&
             ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w2,pw(igrid)%w1,&
             pw(igrid)%w,pw(igrid)%x,pw(igrid)%w3)
          pw(igrid)%wb=>pw(igrid)%w1
        end do
    !$OMP END PARALLEL DO
        ! fix conservation of AMR grid by replacing flux from finer neighbors
        if (fix_conserve_at_step) then
          call recvflux(1,ndim)
          call sendflux(1,ndim)
          call fix_conserve(1,ndim,e_,1)
        end if
        call getbc(global_time,0.d0,e_-1,1)
        evenstep=.true.
      end if 
    end do
    if(evenstep) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,e_)=pw(igrid)%w1(ixGlo1:ixGhi1,&
           ixGlo2:ixGhi2,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do 
    else
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,e_)=pw(igrid)%w2(ixGlo1:ixGhi1,&
           ixGlo2:ixGhi2,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do 
    end if
    deallocate(bj)
    ! point bc mpi data type back to full type for (M)HD
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    bcphys=.true.
  
  end subroutine do_thermal_conduction

  subroutine evolve_stepj(igrid,qcmu,qcmut,qcnu,qcnut,qdt,ixImin1,ixImin2,&
     ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,w1,w2,w,x,wold)
    use mod_global_parameters
    use mod_fix_conserve
    
    integer, intent(in) :: igrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: qcmu,qcmut,qcnu,qcnut,qdt
    double precision, intent(in) :: w1(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),wold(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w2(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:1,1:ndim)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call phys_get_heatconduct(tmp,tmp1,tmp2,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,w1,x,fC)

    w2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=qcmu*w1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,e_)+qcnu*w2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+(1.d0-qcmu-qcnu)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+qcmut*qdt*tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+qcnut*wold(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)
    
    if (fix_conserve_at_step) then
      fC=qcmut*qdt*fC
      call storeflux(igrid,fC,1,ndim,1)
    end if

  end subroutine evolve_stepj

  subroutine evolve_step1(igrid,qcmut,qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,w1,w,x,wold)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: igrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: qcmut, qdt, w(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) ::w1(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       wold(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:1,1:ndim)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer :: lowindex(ndim), ix1,ix2

    call phys_get_heatconduct(tmp,tmp1,tmp2,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x,fC)

    wold(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=qdt*tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! update internal energy
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = tmp1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + qcmut*wold(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)
    
    ! ensure you never trigger negative pressure 
    ! hence code up energy change with respect to kinetic and magnetic
    ! part(nonthermal)
    if(small_values_method=='error') then
      if(any(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<small_e).and. .not.crash) &
         then
        lowindex=minloc(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        lowindex(1)=lowindex(1)+ixOmin1-1;lowindex(2)=lowindex(2)+ixOmin2-1;
        write(*,*)'too small internal energy = ',minval(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)),'at x=',x(lowindex(1),lowindex(2),1:ndim),&
           lowindex,' with limit=',small_e,' on time=',global_time,' step=',it,&
            'where w(1:nwflux)=',w(lowindex(1),lowindex(2),1:nwflux)
        crash=.true.
      else
        if(solve_internal_e) then
          w1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = tmp1(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
        else
          w1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = tmp2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        end if
      end if
    else
      where(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<small_e)
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=small_e
      endwhere
      if(solve_internal_e) then
        w1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)
      else
        w1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = tmp2(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
    end if

    if (fix_conserve_at_step) then
      fC=qcmut*qdt*fC
      call storeflux(igrid,fC,1,ndim,1)
    end if
  
  end subroutine evolve_step1

  !> anisotropic thermal conduction with slope limited symmetric scheme
  !> Sharma 2007 Journal of Computational Physics 227, 123
  subroutine mhd_get_heatconduct(qd,tmp1,tmp2,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x,qvec)
    use mod_global_parameters
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    !! qd store the heat conduction energy changing rate
    double precision, intent(out) :: qd(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim) :: qvec
    
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir) :: mf,&
       Bc,Bcf
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim) :: gradT
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: Te,ka,kaf,&
       ke,kef,qdd,qe,Binv,minq,maxq,Bnorm
    double precision :: alpha,dxinv(ndim)
    integer, dimension(ndim) :: lowindex
    integer :: idims,idir,ix1,ix2,ixmin1,ixmin2,ixmax1,ixmax2,ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2,ixAmin1,ixAmin2,ixAmax1,ixAmax2,ixBmin1,ixBmin2,ixBmax1,&
       ixBmax2

    ! coefficient of limiting on normal component
    if(ndim<3) then
      alpha=0.75d0
    else
      alpha=0.85d0
    end if
    ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

    dxinv=1.d0/dxlevel

    ! tmp1 store internal energy
    if(solve_internal_e) then
      tmp1(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)
    else
      ! tmp2 store kinetic+magnetic energy before addition of heat conduction source
      tmp2(ixImin1:ixImax1,ixImin2:ixImax2) = 0.5d0 * (sum(w(ixImin1:ixImax1,&
         ixImin2:ixImax2,mom(:))**2,dim=ndim+1)/w(ixImin1:ixImax1,&
         ixImin2:ixImax2,rho_) + sum(w(ixImin1:ixImax1,ixImin2:ixImax2,&
         mag(:))**2,dim=ndim+1))
      tmp1(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)-tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
    end if

    ! Clip off negative pressure if small_pressure is set
    if(small_values_method=='error') then
       if (any(tmp1(ixImin1:ixImax1,ixImin2:ixImax2)<small_e) .and. &
          .not.crash) then
         lowindex=minloc(tmp1(ixImin1:ixImax1,ixImin2:ixImax2))
         lowindex(1)=lowindex(1)+ixImin1-1;lowindex(2)=lowindex(2)+ixImin2-1;
         write(*,*)'too low internal energy = ',minval(tmp1(ixImin1:ixImax1,&
            ixImin2:ixImax2)),' at x=',x(lowindex(1),lowindex(2),1:ndim),&
            lowindex,' with limit=',small_e,' on time=',global_time, ' it=',it
         write(*,*) 'w',w(lowindex(1),lowindex(2),:)
         crash=.true.
       end if
    else
    do ix2=ixImin2,ixImax2
    do ix1=ixImin1,ixImax1
       if(tmp1(ix1,ix2)<small_e) then
          tmp1(ix1,ix2)=small_e
       end if
    end do
    end do
    end if
    ! compute the temperature
    Te(ixImin1:ixImax1,ixImin2:ixImax2)=tmp1(ixImin1:ixImax1,&
       ixImin2:ixImax2)*(tc_gamma-one)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    ! B vector
    if(B0field) then
      mf(ixImin1:ixImax1,ixImin2:ixImax2,:)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         mag(:))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,:,0)
    else
      mf(ixImin1:ixImax1,ixImin2:ixImax2,:)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         mag(:));
    end if
    ! |B|
    Binv(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(sum(mf(ixmin1:ixmax1,ixmin2:ixmax2,&
       :)**2,dim=ndim+1))
    where(Binv(ixmin1:ixmax1,ixmin2:ixmax2)/=0.d0)
      Binv(ixmin1:ixmax1,ixmin2:ixmax2)=1.d0/Binv(ixmin1:ixmax1,ixmin2:ixmax2)
    elsewhere
      Binv(ixmin1:ixmax1,ixmin2:ixmax2)=bigdouble
    end where
    ! b unit vector: magnetic field direction vector
    do idims=1,ndim
      mf(ixmin1:ixmax1,ixmin2:ixmax2,idims)=mf(ixmin1:ixmax1,ixmin2:ixmax2,&
         idims)*Binv(ixmin1:ixmax1,ixmin2:ixmax2)
    end do
    ! ixC is cell-corner index
    ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;
    ! b unit vector at cell corner
    Bc=0.d0
    do ix2=0,1
    do ix1=0,1
      ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;
      ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;
      Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:ndim)=Bc(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1:ndim)+mf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndim)
    end do
    end do
    Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:ndim)=Bc(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,1:ndim)*0.5d0**ndim
    if(tc_constant) then
      if(tc_perpendicular) then
        ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=tc_k_para-tc_k_perp
        ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=tc_k_perp
      else
        ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=tc_k_para
      end if
    else
      ! conductivity at cell center
      minq(ixmin1:ixmax1,ixmin2:ixmax2)=tc_k_para*Te(ixmin1:ixmax1,&
         ixmin2:ixmax2)**2.5d0
      ka=0.d0
      do ix2=0,1
      do ix1=0,1
        ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;
        ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;
        ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ka(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+minq(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
      end do
      end do
      ! cell corner conductivity
      ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0**ndim*ka(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      ! compensate with perpendicular conductivity
      if(tc_perpendicular) then
        minq(ixmin1:ixmax1,ixmin2:ixmax2)=tc_k_perp*w(ixmin1:ixmax1,&
           ixmin2:ixmax2,rho_)**2*Binv(ixmin1:ixmax1,&
           ixmin2:ixmax2)**2/dsqrt(Te(ixmin1:ixmax1,ixmin2:ixmax2))
        ke=0.d0
        do ix2=0,1
        do ix1=0,1
          ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;
          ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;
          ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ke(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)+minq(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
        end do
        end do
        ! cell corner conductivity: k_parallel-k_perpendicular
        ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0**ndim*ke(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
        where(ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2)<ka(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))
          ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ka(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        elsewhere
          ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.d0
          ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ka(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)
        end where
      end if
    end if
    ! T gradient at cell faces
    gradT=0.d0
    do idims=1,ndim
      ixBmin1=ixmin1;ixBmin2=ixmin2;
      ixBmax1=ixmax1-kr(idims,1);ixBmax2=ixmax2-kr(idims,2);
      ixAmin1=ixBmin1+kr(idims,1);ixAmin2=ixBmin2+kr(idims,2)
      ixAmax1=ixBmax1+kr(idims,1);ixAmax2=ixBmax2+kr(idims,2);
      gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idims)=(Te(ixAmin1:ixAmax1,&
         ixAmin2:ixAmax2)-Te(ixBmin1:ixBmax1,ixBmin2:ixBmax2))*dxinv(idims)
    end do
    if(tc_slope_limiter=='no') then
      ! calculate thermal conduction flux with symmetric scheme
      do idims=1,ndim
        qd=0.d0
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims ) then
             ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;
             ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;
             qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qd(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)+gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idims)
           end if
        end do
        end do
        ! temperature gradient at cell corner
        qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idims)=qd(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*0.5d0**(ndim-1)
      end do
      ! b grad T at cell corner
      qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sum(qvec(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1:ndim)*Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:ndim),&
         dim=ndim+1)
      do idims=1,ndim
        ! TC flux at cell corner
        gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idims)=ka(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idims)*qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        if(tc_perpendicular) gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idims)=gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idims)+kef(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*qvec(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idims)
      end do
      ! TC flux at cell face
      qvec=0.d0
      do idims=1,ndim
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmax1=ixOmax1-kr(idims,1);ixBmax2=ixOmax2-kr(idims,2);
        ixAmax1=ixOmax1;ixAmax2=ixOmax2; ixAmin1=ixBmin1;ixAmin2=ixBmin2;
        ixBmin1=ixAmin1;ixBmin2=ixAmin2;ixBmax1=ixAmax1;ixBmax2=ixAmax2;
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims ) then
             ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2; 
             ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2; 
             qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=qvec(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2,idims)+gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                idims)
           end if
        end do
        end do
        qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=qvec(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,idims)*0.5d0**(ndim-1)
        if(tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c is isothermal sound speed
          Bcf=0.d0
          do ix2=0,1 
          do ix1=0,1 
             if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims ) then
               ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;
               ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;
               Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=Bcf(ixAmin1:ixAmax1,&
                  ixAmin2:ixAmax2,idims)+Bc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                  idims)
             end if
          end do
          end do
          ! averaged b at face centers
          Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=Bcf(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2,idims)*0.5d0**(ndim-1)
          ixBmin1=ixAmin1+kr(idims,1);ixBmin2=ixAmin2+kr(idims,2)
          ixBmax1=ixAmax1+kr(idims,1);ixBmax2=ixAmax2+kr(idims,2);
          qd(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=2.75d0*(w(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2,rho_)+w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             rho_))*dsqrt(0.5d0*(Te(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2)+Te(ixBmin1:ixBmax1,&
             ixBmin2:ixBmax2)))**3*dabs(Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             idims))
         do ix2=ixAmin2,ixAmax2
         do ix1=ixAmin1,ixAmax1
            if(dabs(qvec(ix1,ix2,idims))>qd(ix1,ix2)) then
              qvec(ix1,ix2,idims)=sign(1.d0,qvec(ix1,ix2,idims))*qd(ix1,ix2)
            end if
         end do
         end do
        end if
      end do
    else
      ! calculate thermal conduction flux with slope-limited symmetric scheme
      qvec=0.d0
      do idims=1,ndim
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmax1=ixOmax1-kr(idims,1);ixBmax2=ixOmax2-kr(idims,2);
        ixAmax1=ixOmax1;ixAmax2=ixOmax2; ixAmin1=ixBmin1;ixAmin2=ixBmin2;
        ! calculate normal of magnetic field
        ixBmin1=ixAmin1+kr(idims,1);ixBmin2=ixAmin2+kr(idims,2)
        ixBmax1=ixAmax1+kr(idims,1);ixBmax2=ixAmax2+kr(idims,2);
        Bnorm(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=0.5d0*(mf(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,idims)+mf(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idims))
        Bcf=0.d0
        kaf=0.d0
        kef=0.d0
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims ) then
             ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;
             ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;
             Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndim)=Bcf(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2,1:ndim)+Bc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                1:ndim)
             kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=kaf(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2)+ka(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
             if(tc_perpendicular) kef(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2)=kef(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2)+ke(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
           end if
        end do
        end do
        ! averaged b at face centers
        Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndim)=Bcf(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,1:ndim)*0.5d0**(ndim-1)
        ! averaged thermal conductivity at face centers
        kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=kaf(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2)*0.5d0**(ndim-1)
        if(tc_perpendicular) kef(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2)=kef(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2)*0.5d0**(ndim-1)
        ! limited normal component
        minq(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=min(alpha*gradT(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,idims),gradT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)/alpha)
        maxq(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=max(alpha*gradT(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,idims),gradT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)/alpha)
        ! eq (19)
        qdd=0.d0
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims ) then
             ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;
             ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;
             qdd(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qdd(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)+gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idims)
           end if
        end do
        end do
        ! temperature gradient at cell corner
        qdd(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qdd(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*0.5d0**(ndim-1)
        ! eq (21)
        qe=0.d0
        do ix2=0,1 
        do ix1=0,1 
           qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qdd(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2)
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims ) then
             ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;
             ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;
             where(qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2)<=minq(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2))
               qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2)=minq(ixAmin1:ixAmax1,&
                  ixAmin2:ixAmax2)
             elsewhere(qd(ixBmin1:ixBmax1,&
                ixBmin2:ixBmax2)>=maxq(ixAmin1:ixAmax1,ixAmin2:ixAmax2))
               qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2)=maxq(ixAmin1:ixAmax1,&
                  ixAmin2:ixAmax2)
             end where
             qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=qvec(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2,idims)+Bc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                idims)**2*qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
             if(tc_perpendicular) qe(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2)=qe(ixAmin1:ixAmax1,&
                ixAmin2:ixAmax2)+qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2) 
           end if
        end do
        end do
        qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=kaf(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2)*qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)*0.5d0**(ndim-1)
        ! add normal flux from perpendicular conduction
        if(tc_perpendicular) qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)+kef(ixAmin1:ixAmax1,ixAmin2:ixAmax2)*qe(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2)*0.5d0**(ndim-1)
        ! limited transverse component, eq (17)
        ixBmin1=ixAmin1;ixBmin2=ixAmin2;
        ixBmax1=ixAmax1+kr(idims,1);ixBmax2=ixAmax2+kr(idims,2);
        do idir=1,ndim
          if(idir==idims) cycle
          qd(ixImin1:ixImax1,ixImin2:ixImax2)=slope_limiter(gradT(&
             ixImin1:ixImax1,ixImin2:ixImax2,idir),ixImin1,ixImin2,ixImax1,&
             ixImax2,ixBmin1,ixBmin2,ixBmax1,ixBmax2,idir,-1)
          qd(ixImin1:ixImax1,ixImin2:ixImax2)=slope_limiter(qd,ixImin1,ixImin2,&
             ixImax1,ixImax2,ixAmin1,ixAmin2,ixAmax1,ixAmax2,idims,1)
          qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=qvec(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2,idims)+kaf(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2)*Bnorm(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2)*Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             idir)*qd(ixAmin1:ixAmax1,ixAmin2:ixAmax2)
        end do

        ! consider magnetic null point
        !where(Binv(ixA^S)==0.d0)
        !  qvec(ixA^S,idims)=tc_k_para*(0.5d0*(Te(ixA^S)+Te(ixB^S)))**2.5d0*gradT(ixA^S,idims)
        !end where

        if(tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c is isothermal sound speed
          ixBmin1=ixAmin1+kr(idims,1);ixBmin2=ixAmin2+kr(idims,2)
          ixBmax1=ixAmax1+kr(idims,1);ixBmax2=ixAmax2+kr(idims,2);
          qd(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=2.75d0*(w(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2,rho_)+w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             rho_))*dsqrt(0.5d0*(Te(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2)+Te(ixBmin1:ixBmax1,&
             ixBmin2:ixBmax2)))**3*dabs(Bnorm(ixAmin1:ixAmax1,&
             ixAmin2:ixAmax2))
         do ix2=ixAmin2,ixAmax2
         do ix1=ixAmin1,ixAmax1
            if(dabs(qvec(ix1,ix2,idims))>qd(ix1,ix2)) then
        !      write(*,*) 'it',it,qvec(ix^D,idims),qd(ix^D),' TC saturated at ',&
        !      x(ix^D,:),' rho',w(ix^D,rho_),' Te',Te(ix^D)
              qvec(ix1,ix2,idims)=sign(1.d0,qvec(ix1,ix2,idims))*qd(ix1,ix2)
            end if
         end do
         end do
        end if
      end do
    end if

    qd=0.d0
    do idims=1,ndim
      qvec(ixmin1:ixmax1,ixmin2:ixmax2,idims)=dxinv(idims)*qvec(ixmin1:ixmax1,&
         ixmin2:ixmax2,idims)
      ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
      ixBmax1=ixOmax1-kr(idims,1);ixBmax2=ixOmax2-kr(idims,2);
      qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=qd(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idims)-qvec(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idims)
    end do
    
  end subroutine mhd_get_heatconduct

  function slope_limiter(f,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idims,pm) result(lf)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idims, pm
    double precision, intent(in) :: f(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: lf(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: signf(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: ixBmin1,ixBmin2,ixBmax1,ixBmax2

    ixBmin1=ixOmin1+pm*kr(idims,1);ixBmin2=ixOmin2+pm*kr(idims,2)
    ixBmax1=ixOmax1+pm*kr(idims,1);ixBmax2=ixOmax2+pm*kr(idims,2);
    signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sign(1.d0,f(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
    select case(tc_slope_limiter)
     case('minmod')
       ! minmod limiter
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*max(0.d0,min(abs(f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)),signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2)))
     case ('MC')
       ! montonized central limiter Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=two*signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)* max(zero,min(dabs(f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)),signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2),&
          signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*quarter*(f(ixBmin1:ixBmax1,&
          ixBmin2:ixBmax2)+f(ixOmin1:ixOmax1,ixOmin2:ixOmax2))))
     case ('superbee')
       ! Roes superbee limiter (eq.3.51i)
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)* max(zero,min(two*dabs(f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)),signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2)),&
          min(dabs(f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
          two*signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*f(ixBmin1:ixBmax1,&
          ixBmin2:ixBmax2)))
     case ('koren')
       ! Barry Koren Right variant
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)* max(zero,min(two*dabs(f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)),two*signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2),&
          (two*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2)*signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+dabs(f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))*third))
     case default
       call mpistop("Unknown slope limiter for thermal conduction")
    end select

  end function slope_limiter

  subroutine mhd_getdt_heatconduct(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    !Check diffusion time limit dt < tc_dtpar*dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_physics
    
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: dx1,dx2, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    ! note that depending on small_values_method=='error' etc, w values may change 
    ! through call to getpthermal
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        dtnew
    
    double precision :: dxinv(1:ndim),mf(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)
    double precision :: tmp2(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2),Te(ixImin1:ixImax1,&
       ixImin2:ixImax2),B2(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: dtdiff_tcond, dtdiff_tsat
    integer          :: idim,idir,ix1,ix2

    dxinv(1)=one/dx1;dxinv(2)=one/dx2;
    
    call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,tmp)

    !temperature
    Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    !tc_k_para_i
    if(tc_constant) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tc_k_para
    else
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tc_k_para*dsqrt(Te(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**5)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    end if
    
    ! B
    if(B0field) then
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,0)
    else
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(:))
    end if
    ! B^-2
    B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum(mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       :)**2,dim=ndim+1)
    ! B_i**2/B**2
    where(B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0.d0)
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)**2/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)**2/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2);
    elsewhere
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=1.d0
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=1.d0;
    end where

    if(tc_saturate) B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=22.d0*dsqrt(Te(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2))

    do idim=1,ndim
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)
      if(tc_saturate) then
        where(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>B2(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=B2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
        end where
      end if
      ! dt< tc_dtpar * dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
      dtdiff_tcond=tc_dtpar/(tc_gamma-1.d0)/maxval(tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*dxinv(idim)**2)
      ! limit the time step
      dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  
  end subroutine mhd_getdt_heatconduct

  subroutine hd_get_heatconduct(tmp,tmp1,tmp2,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x,qvec)
    use mod_global_parameters
    use mod_small_values, only: small_values_method
    
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    !! tmp store the heat conduction energy changing rate
    double precision, intent(out) :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2),qflux(ixImin1:ixImax1,&
       ixImin2:ixImax2),qsatflux(ixImin1:ixImax1,ixImin2:ixImax2)
    integer:: ixmin1,ixmin2,ixmax1,ixmax2,idims,ix1,ix2
    integer, dimension(ndim)       :: lowindex
    
    ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Need two extra layers for thermal conduction")
    
    ! store old kinetic energy
    tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=half*sum(w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(:))**2,dim=ndim+1)/w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)
    ! store old internal energy
    tmp1(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       e_)-tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
    ! Clip off negative pressure if small_pressure is set
    if(small_values_method=='error' .and. .not.crash) then
      if(any(tmp1(ixImin1:ixImax1,ixImin2:ixImax2)<small_e)) then
        lowindex=minloc(tmp1(ixImin1:ixImax1,ixImin2:ixImax2))
        write(*,*)'low internal energy= ',minval(tmp1(ixImin1:ixImax1,&
           ixImin2:ixImax2)),' at x=',x(lowindex(1),lowindex(2),1:ndim),&
           lowindex,' with limit=',small_e,' on time=',global_time
        crash=.true.
      end if
    else
       do ix2=ixImin2,ixImax2
       do ix1=ixImin1,ixImax1
         if(tmp1(ix1,ix2)<small_e) then
          tmp1(ix1,ix2)=small_e
         end if
       end do
       end do
    end if

    ! compute temperature before source addition
    Te(ixImin1:ixImax1,ixImin2:ixImax2)=tmp1(ixImin1:ixImax1,&
       ixImin2:ixImax2)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)*(tc_gamma-one)

    ! compute grad T and store grad T vector
    do idims=1,ndim
      ! idirth component of gradient of temperature at cell center
       select case(typegrad)
       case("central")
         ixmin1=ixImin1+1;ixmin2=ixImin2+1;ixmax1=ixImax1-1;ixmax2=ixImax2-1;
         call gradient(Te,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
            ixmax2,idims,tmp)
       case("limited")
         ixmin1=ixImin1+2;ixmin2=ixImin2+2;ixmax1=ixImax1-2;ixmax2=ixImax2-2;
         call gradientS(Te,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
            ixmax1,ixmax2,idims,tmp)
       end select
       qvec(ixmin1:ixmax1,ixmin2:ixmax2,idims)=tmp(ixmin1:ixmax1,&
          ixmin2:ixmax2)*tc_k_para*dsqrt(Te(ixmin1:ixmax1,ixmin2:ixmax2)**5)
    end do

    if(tc_saturate) then
      ! consider saturation with unsigned saturated TC flux = 5 phi rho c**3
      qsatflux(ixmin1:ixmax1,ixmin2:ixmax2)=5.d0*w(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)*dsqrt(Te(ixmin1:ixmax1,ixmin2:ixmax2)**3)
      qflux(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(qvec(ixmin1:ixmax1,&
         ixmin2:ixmax2,1)**2+qvec(ixmin1:ixmax1,ixmin2:ixmax2,2)**2)
      do ix2=ixmin2,ixmax2
      do ix1=ixmin1,ixmax1
        if(qflux(ix1,ix2)>qsatflux(ix1,ix2)) then
          qsatflux(ix1,ix2)=qsatflux(ix1,ix2)/qflux(ix1,ix2)
          do idims=1,ndim
            qvec(ix1,ix2,idims)=qsatflux(ix1,ix2)*qvec(ix1,ix2,idims)
          end do
        end if
      end do
      end do
    end if

    select case(typediv)
      case("central")
        call divvector(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,tmp)
      case("limited")
        call divvectorS(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,tmp)
    end select

  end subroutine hd_get_heatconduct

  subroutine hd_getdt_heatconduct(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    ! Check diffusion time limit dt < tc_dtpar * dx_i**2 / ((gamma-1)*tc_k_para_i/rho)
    use mod_global_parameters
    use mod_physics
    
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: dx1,dx2, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        dtnew
    
    double precision :: dxinv(1:ndim), tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
        Te(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: dtdiff_tcond,dtdiff_tsat
    integer          :: idim,ix1,ix2

    dxinv(1)=one/dx1;dxinv(2)=one/dx2;

    call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,tmp)

    Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tc_gamma-&
       one)*tc_k_para*dsqrt((Te(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))**5)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    
    do idim=1,ndim
       ! dt< tc_dtpar * dx_idim**2/((gamma-1)*tc_k_para_idim/rho)
       dtdiff_tcond=tc_dtpar/maxval(tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*dxinv(idim)**2)
       if(tc_saturate) then
         ! dt< tc_dtpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
         dtdiff_tsat=tc_dtpar/maxval((tc_gamma-1.d0)*dsqrt(Te(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2))*5.d0*dxinv(idim)**2)
         ! choose the slower flux (bigger time scale) between classic and saturated
         dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
       end if
       ! limit the time step
       dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  
  end subroutine hd_getdt_heatconduct

end module mod_thermal_conduction
