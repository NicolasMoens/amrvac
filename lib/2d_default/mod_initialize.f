!> This module handles the initialization of various components of amrvac
module mod_initialize

  implicit none
  private

  logical :: initialized_already = .false.

  ! Public methods
  public :: initialize_amrvac

contains

  !> Initialize amrvac: read par files and initialize variables
  subroutine initialize_amrvac()
    use mod_input_output
    use mod_physics, only: phys_check, phys_check_params
    use mod_usr_methods, only: usr_set_parameters

    if (initialized_already) return

    ! Check whether the user has loaded a physics module
    call phys_check()

    ! Read input files
    call read_par_files()
    call initialize_vars()
    call init_comm_types()

    if(associated(usr_set_parameters)) call usr_set_parameters()

    call phys_check_params()

    initialized_already = .true.
  end subroutine initialize_amrvac

  !> Initialize (and allocate) simulation and grid variables
  !> @todo Explain which ones are not initialized here
  subroutine initialize_vars
    use mod_forest
    use mod_global_parameters
    use mod_ghostcells_update

    integer :: igrid, level, ipe, ig1,ig2
    logical :: ok

    allocate(pw(max_blocks))
    allocate(pw_sub(max_blocks))
    allocate(neighbor(2,-1:1,-1:1,max_blocks),neighbor_child(2,0:3,0:3,&
       max_blocks))
    allocate(neighbor_type(-1:1,-1:1,max_blocks),neighbor_active(-1:1,-1:1,&
       max_blocks))
    if (phi_ > 0) allocate(neighbor_pole(-1:1,-1:1,max_blocks))
    allocate(igrids(max_blocks),igrids_active(max_blocks),&
       igrids_passive(max_blocks))
    allocate(rnode(rnodehi,max_blocks),rnode_sub(rnodehi,max_blocks),&
       dt_grid(max_blocks))
    allocate(node(nodehi,max_blocks),node_sub(nodehi,max_blocks),&
       phyboundblock(max_blocks))
    allocate(pflux(2,2,max_blocks))
    ! allocate mesh for particles
    if(use_particles) allocate(gridvars(max_blocks))

    it=it_init
    global_time=time_init

    dt=zero

    ! set all dt to zero
    dt_grid(1:max_blocks)=zero

    ! check resolution
    if (mod(ixGhi1,2)/=0.or.mod(ixGhi2,2)/=0) then
       call mpistop("mesh widths must give even number grid points")
    end if
    ixMlo1=ixGlo1+nghostcells;ixMlo2=ixGlo2+nghostcells
    ixMhi1=ixGhi1-nghostcells;ixMhi2=ixGhi2-nghostcells;

    if (nbufferx1>(ixMhi1-ixMlo1+1).or.nbufferx2>(ixMhi2-ixMlo2+1)) then
       write(unitterm,*) "nbufferx^D bigger than mesh size makes no sense."
       write(unitterm,*) "Decrease nbufferx or increase mesh size"
       call mpistop("")
    end if

    ! initialize dx arrays on finer (>1) levels
    do level=2,refine_max_level
       dx(1,level) = dx(1,level-1) * half
       dx(2,level) = dx(2,level-1) * half  ! refine ratio 2
    end do

    ! domain decomposition
    ! physical extent of a grid block at level 1, per dimension
    dg1(1)=dx(1,1)*dble(block_nx1)
    dg2(1)=dx(2,1)*dble(block_nx2)
    ! number of grid blocks at level 1 in simulation domain, per dimension
    ng1(1)=nint((xprobmax1-xprobmin1)/dg1(1))
    ng2(1)=nint((xprobmax2-xprobmin2)/dg2(1))
    ! total number of grid blocks at level 1
    nglev1=ng1(1)*ng2(1)

    do level=2,refine_max_level
       dg1(level)=half*dg1(level-1);dg2(level)=half*dg2(level-1);
       ng1(level)=ng1(level-1)*2;ng2(level)=ng2(level-1)*2;
    end do

    ! check that specified stepsize correctly divides domain
    ok=((abs(dble(ng1(1))*dg1(1)-(xprobmax1-xprobmin1))<=&
       smalldouble).and.(abs(dble(ng2(1))*dg2(1)-(xprobmax2-xprobmin2))<=&
       smalldouble))
    if (.not.ok) then
       write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
       call mpistop("domain cannot be divided by meshes of given gridsize")
    end if


    poleB=.false.
    if (.not.slab) call set_pole

    ! initialize connectivity data
    igridstail=0

    ! allocate memory for forest data structures
    allocate(level_head(refine_max_level),level_tail(refine_max_level))
    do level=1,refine_max_level
       nullify(level_head(level)%node,level_tail(level)%node)
    end do

    allocate(igrid_to_node(max_blocks,0:npe-1))
    do ipe=0,npe-1
       do igrid=1,max_blocks
          nullify(igrid_to_node(igrid,ipe)%node)
       end do
    end do

    allocate(sfc(1:3,max_blocks*npe))

    allocate(igrid_to_sfc(max_blocks))

    sfc=0
    allocate(Morton_start(0:npe-1),Morton_stop(0:npe-1))
    allocate(Morton_sub_start(0:npe-1),Morton_sub_stop(0:npe-1))

    allocate(nleafs_level(1:nlevelshi))

    allocate(coarsen(max_blocks,0:npe-1),refine(max_blocks,0:npe-1))
    coarsen=.false.
    refine=.false.
    if (nbufferx1/=0.or.nbufferx2/=0) then
       allocate(buffer(max_blocks,0:npe-1))
       buffer=.false.
    end if
    allocate(igrid_inuse(max_blocks,0:npe-1))
    igrid_inuse=.false.

    allocate(tree_root(1:ng1(1),1:ng2(1)))
    do ig2=1,ng2(1)
    do ig1=1,ng1(1)
    nullify(tree_root(ig1,ig2)%node)
    end do
    end do

    ! define index ranges and MPI send/receive derived datatype for ghost-cell swap
    call init_bc()
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    call create_bc_mpi_datatype(0,nwflux+nwaux)

  end subroutine initialize_vars


end module mod_initialize
