!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics
  use mod_global_parameters, only: name_len
  use mod_physics_hllc
  use mod_physics_roe
  use mod_physics_ppm

  implicit none
  public

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = ""

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0

  !> Whether the physics routines require diagonal ghost cells, for example for
  !> computing a curl.
  logical :: phys_req_diagonal = .true.

  !> Array per direction per variable, which can be used to specify that certain
  !> fluxes have to be treated differently
  integer, allocatable :: flux_type(:, :)

  !> Indicates a normal flux
  integer, parameter   :: flux_default        = 0
  !> Indicates the flux should be treated with tvdlf
  integer, parameter   :: flux_tvdlf          = 1
  !> Indicates dissipation should be omitted
  integer, parameter   :: flux_no_dissipation = 2


  procedure(sub_check_params), pointer    :: phys_check_params           => &
     null()
  procedure(sub_convert), pointer         :: phys_to_conserved           => &
     null()
  procedure(sub_convert), pointer         :: phys_to_primitive           => &
     null()
  procedure(sub_modify_wLR), pointer      :: phys_modify_wLR             => &
     null()
  procedure(sub_get_cmax), pointer        :: phys_get_cmax               => &
     null()
  procedure(sub_get_cbounds), pointer     :: phys_get_cbounds            => &
     null()
  procedure(sub_get_flux), pointer        :: phys_get_flux               => &
     null()
  procedure(sub_get_dt), pointer          :: phys_get_dt                 => &
     null()
  procedure(sub_add_source_geom), pointer :: phys_add_source_geom        => &
     null()
  procedure(sub_add_source), pointer      :: phys_add_source             => &
     null()
  procedure(sub_get_aux), pointer         :: phys_get_aux                => &
     null()
  procedure(sub_check_w), pointer         :: phys_check_w                => &
     null()
  procedure(sub_get_pthermal), pointer    :: phys_get_pthermal           => &
     null()
  procedure(sub_boundary_adjust), pointer :: phys_boundary_adjust        => &
     null()
  procedure(sub_write_info), pointer      :: phys_write_info             => &
     null()
  procedure(sub_angmomfix), pointer       :: phys_angmomfix              => &
     null()
  procedure(sub_small_values), pointer    :: phys_handle_small_values    => &
     null()

  abstract interface

     subroutine sub_check_params
     end subroutine sub_check_params

     subroutine sub_boundary_adjust
       use mod_global_parameters
     end subroutine sub_boundary_adjust

     subroutine sub_convert(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
           nw)
       double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:2)
     end subroutine sub_convert

     subroutine sub_modify_wLR(wLC, wRC, ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
       double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
          1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
     end subroutine sub_modify_wLR

     subroutine sub_get_cmax(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2, idim, cmax)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
       double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
           nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
       double precision, intent(inout) :: cmax(ixImin1:ixImax1,&
          ixImin2:ixImax2)
     end subroutine sub_get_cmax

     subroutine sub_get_cbounds(wLC, wRC, wLp, wRp, x, ixImin1,ixImin2,ixImax1,&
        ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
       double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
           nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2, nw)
       double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
           nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2, nw)
       double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:2)
       double precision, intent(inout) :: cmax(ixImin1:ixImax1,&
          ixImin2:ixImax2)
       double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
          ixImin2:ixImax2)
     end subroutine sub_get_cbounds

     subroutine sub_get_flux(wC, w, x, ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, f)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
       double precision, intent(in)    :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:nw)
       double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:nw)
       double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:2)
       double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
           nwflux)
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2, wCT, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,&
          ixImin2:ixImax2, 1:2)
       double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
     end subroutine sub_add_source_geom

     subroutine sub_add_source(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2, wCT, w, x, qsourcesplit, active)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(in)    :: qdt
       double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
       double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:nw)
       logical, intent(in)             :: qsourcesplit
       logical, intent(inout)          :: active !< Needs to be set to true when active
     end subroutine sub_add_source

     subroutine sub_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
          ixImin2:ixImax2, 1:2)
       double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
           1:nw)
       double precision, intent(inout) :: dtnew
     end subroutine sub_get_dt

     subroutine sub_get_aux(clipping,w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2,subname)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
          1:ndim)
       double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
          nw)
       logical, intent(in)             :: clipping
       character(len=*)                :: subname
     end subroutine sub_get_aux

     subroutine sub_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2, w, w_flag)
       use mod_global_parameters
       logical, intent(in)          :: primitive
       integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
       integer, intent(inout)       :: w_flag(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
     end subroutine sub_check_w

     subroutine sub_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2,pth)
       use mod_global_parameters
       integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
       double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
          1:ndim)
       double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2)
     end subroutine sub_get_pthermal

     subroutine sub_write_info(file_handle)
       integer, intent(in) :: file_handle
     end subroutine sub_write_info

     subroutine sub_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim)
       use mod_global_parameters
       integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
          1:ndim)
       double precision, intent(inout)    :: fC(ixImin1:ixImax1,&
          ixImin2:ixImax2,1:nwflux,1:ndim), wnew(ixImin1:ixImax1,&
          ixImin2:ixImax2,1:nw)
       integer, intent(in)                :: idim
     end subroutine sub_angmomfix

     subroutine sub_small_values(primitive, w, x, ixImin1,ixImin2,ixImax1,&
        ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, subname)
       use mod_global_parameters
       use mod_small_values
       logical, intent(in)             :: primitive
       integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2
       double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
          1:nw)
       double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
          1:ndim)
       character(len=*), intent(in)    :: subname
     end subroutine sub_small_values

  end interface

contains

  subroutine phys_check()
    use mod_global_parameters, only: nw, ndir

    use mod_physics_hllc, only: phys_hllc_check
    use mod_physics_roe, only: phys_roe_check
    use mod_physics_ppm, only: phys_ppm_check

    if (physics_type == "") call mpistop("Error: no physics module loaded")

    call phys_hllc_check()
    call phys_roe_check()
    call phys_ppm_check()

    ! Checks whether the required physics methods have been defined
    if (.not. associated(phys_check_params)) phys_check_params => &
       dummy_check_params

    if (.not. associated(phys_to_conserved)) call &
       mpistop("Error: phys_to_conserved not defined")

    if (.not. associated(phys_to_primitive)) call &
       mpistop("Error: phys_to_primitive not defined")

    if (.not. associated(phys_modify_wLR)) phys_modify_wLR => dummy_modify_wLR

    if (.not. associated(phys_get_cmax)) call &
       mpistop("Error: no phys_get_cmax not defined")

    if (.not. associated(phys_get_cbounds)) call &
       mpistop("Error: no phys_get_cbounds not defined")

    if (.not. associated(phys_get_flux)) call &
       mpistop("Error: no phys_get_flux not defined")

    if (.not. associated(phys_get_dt)) call &
       mpistop("Error: no phys_get_dt not defined")

    if (.not. associated(phys_add_source_geom)) phys_add_source_geom => &
       dummy_add_source_geom

    if (.not. associated(phys_add_source)) phys_add_source => dummy_add_source

    if (.not. associated(phys_get_aux)) phys_get_aux => dummy_get_aux

    if (.not. associated(phys_check_w)) phys_check_w => dummy_check_w

    if (.not. associated(phys_get_pthermal)) phys_get_pthermal => &
       dummy_get_pthermal

    if (.not. associated(phys_boundary_adjust)) phys_boundary_adjust => &
       dummy_boundary_adjust

    if (.not. associated(phys_write_info)) phys_write_info => dummy_write_info

    if (.not. associated(phys_angmomfix)) phys_angmomfix => dummy_angmomfix

    if (.not. associated(phys_handle_small_values)) phys_handle_small_values &
       => dummy_small_values

  end subroutine phys_check

  subroutine dummy_init_params
  end subroutine dummy_init_params

  subroutine dummy_check_params
  end subroutine dummy_check_params

  subroutine dummy_modify_wLR(wLC, wRC, ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(inout)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  end subroutine dummy_modify_wLR

  subroutine dummy_add_source_geom(qdt, ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, wCT, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:2)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, wCT, w, x, qsourcesplit, active)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active
    ! Don't have to set active, since it starts as .false.
  end subroutine dummy_add_source

  subroutine dummy_get_aux(clipping,w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,subname)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    logical, intent(in)             :: clipping
    character(len=*)                :: subname
  end subroutine dummy_get_aux

  subroutine dummy_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, w_flag)
    use mod_global_parameters
    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    integer, intent(inout)       :: w_flag(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

    w_flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0             ! All okay
  end subroutine dummy_check_w

  subroutine dummy_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, pth)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2)

    call mpistop("No get_pthermal method specified")
  end subroutine dummy_get_pthermal

  subroutine dummy_boundary_adjust
  end subroutine dummy_boundary_adjust

  subroutine dummy_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh !< File handle
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    ! Number of physics parameters
    integer, parameter                  :: n_par = 0

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine dummy_write_info

  subroutine dummy_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout)    :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim), wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    integer, intent(in)                :: idim
  end subroutine dummy_angmomfix

  subroutine dummy_small_values(primitive, w, x, ixImin1,ixImin2,ixImax1,&
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
  end subroutine dummy_small_values
  
end module mod_physics
