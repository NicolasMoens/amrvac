!> update ghost cells of all blocks including physical boundaries 
module mod_ghostcells_update

  implicit none
  ! Special buffer for pole copy
  type wbuffer
    double precision, dimension(:,:), allocatable :: w
  end type wbuffer

  ! A switch of update physical boundary or not
  logical, public :: bcphys=.true.
  integer :: ixMmin1,ixMmax1, ixCoGmin1,ixCoGmax1, ixCoMmin1,ixCoMmax1

  !> The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  ! index ranges to send (S) to sibling blocks, receive (R) from 
  ! sibling blocks, send restricted (r) ghost cells to coarser blocks 
  integer, dimension(-1:2,-1:1) :: ixS_srl_min1,ixS_srl_max1, ixR_srl_min1,&
     ixR_srl_max1, ixS_r_min1,ixS_r_max1

  ! index ranges to receive restriced ghost cells from finer blocks, 
  ! send prolongated (p) ghost cells to finer blocks, receive prolongated 
  ! ghost from coarser blocks
  integer, dimension(-1:1, 0:3) :: ixR_r_min1,ixR_r_max1, ixS_p_min1,&
     ixS_p_max1, ixR_p_min1,ixR_p_max1

  ! MPI derived datatype to send and receive subarrays of ghost cells to
  ! neighbor blocks in a different processor.
  !
  ! The first index goes from -1:2, where -1 is used when a block touches the
  ! lower boundary, 1 when a block touches an upper boundary, and 0 a situation
  ! away from boundary conditions, 2 when a block touched both lower and upper
  ! boundary
  !
  ! There are two variants, _f indicates that all flux variables are filled,
  ! whereas _p means that part of the variables is filled 
  ! Furthermore _r_ stands for restrict, _p_ for prolongation.
  integer, dimension(-1:2,-1:1), target :: type_send_srl_f, type_recv_srl_f
  integer, dimension(-1:1,-1:1), target :: type_send_r_f
  integer, dimension(-1:1, 0:3), target :: type_recv_r_f, type_send_p_f,&
      type_recv_p_f
  integer, dimension(-1:2,-1:1), target :: type_send_srl_p1, type_recv_srl_p1
  integer, dimension(-1:1,-1:1), target :: type_send_r_p1
  integer, dimension(-1:1, 0:3), target :: type_recv_r_p1, type_send_p_p1,&
      type_recv_p_p1
  integer, dimension(-1:2,-1:1), target :: type_send_srl_p2, type_recv_srl_p2
  integer, dimension(-1:1,-1:1), target :: type_send_r_p2
  integer, dimension(-1:1, 0:3), target :: type_recv_r_p2, type_send_p_p2,&
      type_recv_p_p2
  integer, dimension(:,:), pointer :: type_send_srl, type_recv_srl,&
      type_send_r
  integer, dimension(:,:), pointer :: type_recv_r, type_send_p, type_recv_p

contains

  subroutine init_bc()
    use mod_global_parameters 

    integer :: nghostcellsCo, interpolation_order
    integer :: nx1, nxCo1, ixGmin1,ixGmax1, i1, ic1, inc1, iib1

    ixGmin1=ixGlo1;ixGmax1=ixGhi1;
    ixMmin1=ixGmin1+nghostcells;ixMmax1=ixGmax1-nghostcells;
    ixCoGmin1=1;
    !ixCoGmax^D=ixGmax^D/2+nghostcells;
    ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells;

    ixCoMmin1=ixCoGmin1+nghostcells;ixCoMmax1=ixCoGmax1-nghostcells;

    nx1=ixMmax1-ixMmin1+1;
    nxCo1=nx1/2;
    
    select case (typeghostfill)
    case ("copy")
       interpolation_order=1
    case ("linear","unlimit")
       interpolation_order=2
    case default
       interpolation_order=2
       write (unitterm,*) "Undefined typeghostfill ",typeghostfill
       call mpistop("")
    end select
    nghostcellsCo=int((nghostcells+1)/2)
    
    if (nghostcellsCo+interpolation_order-1>nghostcells) then
       call mpistop("interpolation order for prolongation in getbc to high")
    end if
    
    ! (iib,i) index has following meanings: iib = 0 means it is not at any physical boundary
    ! iib=-1 means it is at the minimum side of a physical boundary  
    ! iib= 1 means it is at the maximum side of a physical boundary  
    ! i=-1 means subregion prepared for the neighbor at its minimum side 
    ! i= 1 means subregion prepared for the neighbor at its maximum side 
    
    ixS_srl_min1(:,-1)=ixMmin1
    ixS_srl_min1(:, 1)=ixMmax1+1-nghostcells
    ixS_srl_max1(:,-1)=ixMmin1-1+nghostcells
    ixS_srl_max1(:, 1)=ixMmax1
    
    ixS_srl_min1(-1,0)=1
    ixS_srl_min1( 0,0)=ixMmin1
    ixS_srl_min1( 1,0)=ixMmin1
    ixS_srl_min1( 2,0)=1
    ixS_srl_max1(-1,0)=ixMmax1
    ixS_srl_max1( 0,0)=ixMmax1
    ixS_srl_max1( 1,0)=ixGmax1
    ixS_srl_max1( 2,0)=ixGmax1
     
    ixR_srl_min1(:,-1)=1
    ixR_srl_min1(:, 1)=ixMmax1+1
    ixR_srl_max1(:,-1)=nghostcells
    ixR_srl_max1(:, 1)=ixGmax1
    
    ixR_srl_min1(-1,0)=1
    ixR_srl_min1( 0,0)=ixMmin1
    ixR_srl_min1( 1,0)=ixMmin1
    ixR_srl_min1( 2,0)=1
    ixR_srl_max1(-1,0)=ixMmax1
    ixR_srl_max1( 0,0)=ixMmax1
    ixR_srl_max1( 1,0)=ixGmax1
    ixR_srl_max1( 2,0)=ixGmax1
    
    ixS_r_min1(:,-1)=ixCoMmin1
    ixS_r_min1(:, 1)=ixCoMmax1+1-nghostcells
    ixS_r_max1(:,-1)=ixCoMmin1-1+nghostcells
    ixS_r_max1(:, 1)=ixCoMmax1
    
    ixS_r_min1(-1,0)=1
    ixS_r_min1( 0,0)=ixCoMmin1
    ixS_r_min1( 1,0)=ixCoMmin1
    ixS_r_max1(-1,0)=ixCoMmax1
    ixS_r_max1( 0,0)=ixCoMmax1
    ixS_r_max1( 1,0)=ixCoGmax1
    
    ixR_r_min1(:, 0)=1
    ixR_r_min1(:, 1)=ixMmin1
    ixR_r_min1(:, 2)=ixMmin1+nxCo1
    ixR_r_min1(:, 3)=ixMmax1+1
    ixR_r_max1(:, 0)=nghostcells
    ixR_r_max1(:, 1)=ixMmin1-1+nxCo1
    ixR_r_max1(:, 2)=ixMmax1
    ixR_r_max1(:, 3)=ixGmax1

    ixR_r_min1(-1,1)=1
    ixR_r_max1(-1,1)=ixMmin1-1+nxCo1
    ixR_r_min1( 1,2)=ixMmin1+nxCo1
    ixR_r_max1( 1,2)=ixGmax1

    ixS_p_min1(:, 0)=ixMmin1-(interpolation_order-1)
    ixS_p_min1(:, 1)=ixMmin1-(interpolation_order-1)
    ixS_p_min1(:, 2)=ixMmin1+nxCo1-nghostcellsCo-(interpolation_order-1)
    ixS_p_min1(:, 3)=ixMmax1+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max1(:, 0)=ixMmin1-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max1(:, 1)=ixMmin1-1+nxCo1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max1(:, 2)=ixMmax1+(interpolation_order-1)
    ixS_p_max1(:, 3)=ixMmax1+(interpolation_order-1)

    ixS_p_min1(-1,1)=1
    ixS_p_max1(-1,1)=ixMmin1-1+nxCo1+nghostcellsCo+(interpolation_order-1)
    ixS_p_min1( 1,2)=ixMmin1+nxCo1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max1( 1,2)=ixGmax1

    ixR_p_min1(:, 0)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 1)=ixCoMmin1-(interpolation_order-1)
    ixR_p_min1(:, 2)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 3)=ixCoMmax1+1-(interpolation_order-1)
    ixR_p_max1(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max1(:, 1)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)
    ixR_p_max1(:, 2)=ixCoMmax1+(interpolation_order-1)
    ixR_p_max1(:, 3)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)

    ixR_p_min1(-1,1)=1
    ixR_p_max1(-1,1)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)
    ixR_p_min1( 1,2)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_max1( 1,2)=ixCoGmax1
    

  end subroutine init_bc

  subroutine create_bc_mpi_datatype(nwstart,nwbc) 
    use mod_global_parameters 

    integer, intent(in) :: nwstart, nwbc
    integer :: i1, ic1, inc1, iib1

    do i1=-1,1
      if (i1==0) cycle
      do iib1=-1,2
         call get_bc_comm_type(type_send_srl(iib1,i1),ixS_srl_min1(iib1,i1),&
            ixS_srl_max1(iib1,i1),ixGlo1,ixGhi1,nwstart,nwbc)
         call get_bc_comm_type(type_recv_srl(iib1,i1),ixR_srl_min1(iib1,i1),&
            ixR_srl_max1(iib1,i1),ixGlo1,ixGhi1,nwstart,nwbc)
         if (iib1==2) cycle
         call get_bc_comm_type(type_send_r(iib1,i1),  ixS_r_min1(iib1,i1),&
            ixS_r_max1(iib1,i1),ixCoGmin1,ixCoGmax1,nwstart,nwbc)
         do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
            inc1=2*i1+ic1
            call get_bc_comm_type(type_recv_r(iib1,inc1),ixR_r_min1(iib1,inc1),&
               ixR_r_max1(iib1,inc1), ixGlo1,ixGhi1,nwstart,nwbc)
            call get_bc_comm_type(type_send_p(iib1,inc1),ixS_p_min1(iib1,inc1),&
               ixS_p_max1(iib1,inc1), ixGlo1,ixGhi1,nwstart,nwbc)
            call get_bc_comm_type(type_recv_p(iib1,inc1),ixR_p_min1(iib1,inc1),&
               ixR_p_max1(iib1,inc1),ixCoGmin1,ixCoGmax1,nwstart,nwbc)
         end do
      end do
    end do
  
  end subroutine create_bc_mpi_datatype

  subroutine get_bc_comm_type(comm_type,ixmin1,ixmax1,ixGmin1,ixGmax1,nwstart,&
     nwbc)
    use mod_global_parameters 
  
    integer, intent(inout) :: comm_type
    integer, intent(in) :: ixmin1,ixmax1, ixGmin1,ixGmax1, nwstart, nwbc
    
    integer, dimension(ndim+1) :: fullsize, subsize, start

    fullsize(1)=ixGmax1;
    fullsize(ndim+1)=nw
    subsize(1)=ixmax1-ixmin1+1;
    subsize(ndim+1)=nwbc
    start(1)=ixmin1-1;
    start(ndim+1)=nwstart
    
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,fullsize,subsize,start,&
       MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
    call MPI_TYPE_COMMIT(comm_type,ierrmpi)
    
  end subroutine get_bc_comm_type

  subroutine put_bc_comm_types()
    use mod_global_parameters 
 
    integer :: i1, ic1, inc1, iib1

    do i1=-1,1
       if (i1==0) cycle
       do iib1=-1,2
           call MPI_TYPE_FREE(type_send_srl(iib1,i1),ierrmpi)
           call MPI_TYPE_FREE(type_recv_srl(iib1,i1),ierrmpi)
           if (levmin==levmax) cycle
           if (iib1==2) cycle
           call MPI_TYPE_FREE(type_send_r(iib1,i1),ierrmpi)
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
              inc1=2*i1+ic1
              call MPI_TYPE_FREE(type_recv_r(iib1,inc1),ierrmpi)
              call MPI_TYPE_FREE(type_send_p(iib1,inc1),ierrmpi)
              call MPI_TYPE_FREE(type_recv_p(iib1,inc1),ierrmpi)
           end do
       end do
    end do
  
  end subroutine put_bc_comm_types

  subroutine getbc(time,qdt,nwstart,nwbc,req_diag)
    use mod_global_parameters
    use mod_physics

    double precision, intent(in)      :: time, qdt
    integer, intent(in)               :: nwstart ! Fill from nwstart+1
    integer, intent(in)               :: nwbc    ! Number of variables to fill
    logical, intent(in), optional     :: req_diag !If false, skip diagonal ghost cells

    integer :: my_neighbor_type, ipole, idims, iside, nwhead, nwtail
    integer :: iigrid, igrid, ineighbor, ipe_neighbor
    integer :: nrecvs, nsends, isizes
    integer :: ixGmin1,ixGmax1, ixRmin1,ixRmax1, ixSmin1,ixSmax1, ixBmin1,&
       ixBmax1, ixImin1,ixImax1, kmin1,kmax1
    integer :: i1, n_i1, ic1, inc1, n_inc1, iib1
    ! store physical boundary indicating index
    integer :: idphyb(ndim,max_blocks),bindex(ndim)
    integer :: isend_buf(npwbuf), ipwbuf, nghostcellsco,iB
    logical  :: req_diagonal
    type(wbuffer) :: pwbuf(npwbuf)

    double precision :: time_bcin
    ! Stretching grid parameters for coarsened block of the current block
    double precision :: logGl,qstl

    nwhead=nwstart+1
    nwtail=nwstart+nwbc

    req_diagonal = .true.
    if (present(req_diag)) req_diagonal = req_diag

    time_bcin=MPI_WTIME()
    ixGmin1=ixGlo1;ixGmax1=ixGhi1;
    
    if (internalboundary) then 
       call getintbc(time,ixGmin1,ixGmax1)
    end if
    ! fill ghost cells in physical boundaries
    if(bcphys) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         if(.not.phyboundblock(igrid)) cycle
         saveigrid=igrid
         block=>pw(igrid)
         dxlevel(1)=rnode(rpdx1_,igrid);
         do idims=1,ndim
            ! to avoid using as yet unknown corner info in more than 1D, we
            ! fill only interior mesh ranges of the ghost cell ranges at first,
            ! and progressively enlarge the ranges to include corners later
            
             kmin1=merge(0, 1, idims==1)
             kmax1=merge(0, 1, idims==1)
             ixBmin1=ixGmin1+kmin1*nghostcells
             ixBmax1=ixGmax1-kmax1*nghostcells
            
            
            
            do iside=1,2
               i1=kr(1,idims)*(2*iside-3);
               if (aperiodB(idims)) then 
                  if (neighbor_type(i1,igrid) /= neighbor_boundary .and. .not. &
                     pw(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
               else 
                  if (neighbor_type(i1,igrid) /= neighbor_boundary) cycle
               end if
               call bc_phys(iside,idims,time,qdt,pw(igrid)%wb,pw(igrid)%x,&
                  ixGmin1,ixGmax1,ixBmin1,ixBmax1)
            end do
         end do
      end do
    end if
    
    ! default : no singular axis
    ipole=0
    
    irecv=0
    isend=0
    isend_buf=0
    ipwbuf=1
    ! total number of times to call MPI_IRECV in each processor between sibling blocks or from finer neighbors
    nrecvs=nrecv_bc_srl+nrecv_bc_r
    ! total number of times to call MPI_ISEND in each processor between sibling blocks or to coarser neighors
    nsends=nsend_bc_srl+nsend_bc_r

    allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
    recvrequest=MPI_REQUEST_NULL

    allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
    sendrequest=MPI_REQUEST_NULL

    ! receiving ghost-cell values from sibling blocks and finer neighbors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       call identifyphysbound(igrid,iib1)   
       idphyb(1,igrid)=iib1;
       do i1=-1,1
          if (skip_direction([ i1 ])) cycle
          my_neighbor_type=neighbor_type(i1,igrid)
          select case (my_neighbor_type)
          case (neighbor_sibling)
             call bc_recv_srl
          case (neighbor_fine)
             call bc_recv_restrict
          end select
       end do
    end do
    
    ! sending ghost-cell values to sibling blocks and coarser neighbors
    nghostcellsco=ceiling(nghostcells*0.5d0)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>pw(igrid)

       ! Used stored data to identify physical boundaries
       iib1=idphyb(1,igrid);

       if (any(neighbor_type(:,igrid)==neighbor_coarse)) then
          dxlevel(1)=rnode(rpdx1_,igrid);
    
          call coarsen_grid(pw(igrid)%wb,pw(igrid)%x,ixGmin1,ixGmax1,ixMmin1,&
             ixMmax1,pw(igrid)%wcoarse,pw(igrid)%xcoarse,ixCoGmin1,ixCoGmax1,&
             ixCoMmin1,ixCoMmax1,igrid,igrid)
       end if
    
       do i1=-1,1
          if (skip_direction([ i1 ])) cycle
          if (phi_ > 0) ipole=neighbor_pole(i1,igrid)
          my_neighbor_type=neighbor_type(i1,igrid)
          select case (my_neighbor_type)
          case (neighbor_coarse)
             call bc_send_restrict
          case (neighbor_sibling)
             call bc_send_srl
          end select
       end do
    end do
    
    !if (irecv/=nrecvs) then
    !   call mpistop("number of recvs in phase1 in amr_ghostcells is incorrect")
    !end if
    !if (isend/=nsends) then
    !   call mpistop("number of sends in phase1 in amr_ghostcells is incorrect")
    !end if
    
    call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
    deallocate(recvstatus,recvrequest)

    call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do
    deallocate(sendstatus,sendrequest)

    irecv=0
    isend=0
    isend_buf=0
    ipwbuf=1
    ! total number of times to call MPI_IRECV in each processor from coarser neighbors
    nrecvs=nrecv_bc_p
    ! total number of times to call MPI_ISEND in each processor to finer neighbors
    nsends=nsend_bc_p

    allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
    recvrequest=MPI_REQUEST_NULL
    allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
    sendrequest=MPI_REQUEST_NULL

    ! receiving ghost-cell values from coarser neighbors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       iib1=idphyb(1,igrid);
       do i1=-1,1
          if (skip_direction([ i1 ])) cycle
          my_neighbor_type=neighbor_type(i1,igrid)
          if (my_neighbor_type==neighbor_coarse) call bc_recv_prolong
       end do
    end do
    ! sending ghost-cell values to finer neighbors 
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>pw(igrid)
       iib1=idphyb(1,igrid);
       dxlevel(1)=rnode(rpdx1_,igrid);
       if (any(neighbor_type(:,igrid)==neighbor_fine)) then
          do i1=-1,1
             if (skip_direction([ i1 ])) cycle
             if (phi_ > 0) ipole=neighbor_pole(i1,igrid)
             my_neighbor_type=neighbor_type(i1,igrid)
             if (my_neighbor_type==neighbor_fine) call bc_send_prolong
          end do
       end if
    end do
    
    !if (irecv/=nrecvs) then
    !   call mpistop("number of recvs in phase2 in amr_ghostcells is incorrect")
    !end if
    !if (isend/=nsends) then
    !   call mpistop("number of sends in phase2 in amr_ghostcells is incorrect")
    !end if
    
    call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
    deallocate(recvstatus,recvrequest)

    ! do prolongation on the ghost-cell values received from coarser neighbors 
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>pw(igrid)
       iib1=idphyb(1,igrid);
       dxlevel(1)=rnode(rpdx1_,igrid);
       if (any(neighbor_type(:,igrid)==neighbor_coarse)) then
          do i1=-1,1
             if (skip_direction([ i1 ])) cycle
             my_neighbor_type=neighbor_type(i1,igrid)
             if (my_neighbor_type==neighbor_coarse) call bc_prolong
          end do
       end if
    end do
    
    call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do
    deallocate(sendstatus,sendrequest)

     ! modify normal component of magnetic field to fix divB=0 
    if(bcphys .and. physics_type=='mhd' .and. ndim>1) call &
       phys_boundary_adjust()
    
    if (nwaux>0) call fix_auxiliary
    
    time_bc=time_bc+(MPI_WTIME()-time_bcin)
    
    contains

      logical function skip_direction(dir)
        integer, intent(in) :: dir(1)

        if (all(dir == 0)) then
           skip_direction = .true.
        else if (.not. req_diagonal .and. count(dir /= 0) > 1) then
           skip_direction = .true.
        else
           skip_direction = .false.
        end if
      end function skip_direction

      !> Send to sibling at same refinement level
      subroutine bc_send_srl

        ineighbor=neighbor(1,i1,igrid)
        ipe_neighbor=neighbor(2,i1,igrid)

        if (ipole==0) then
           n_i1=-i1;
           if (ipe_neighbor==mype) then
              ixSmin1=ixS_srl_min1(iib1,i1);ixSmax1=ixS_srl_max1(iib1,i1);
              ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmax1=ixR_srl_max1(iib1,n_i1);
              pw(ineighbor)%wb(ixRmin1:ixRmax1,&
                 nwhead:nwtail)=pw(igrid)%wb(ixSmin1:ixSmax1,nwhead:nwtail)
           else
              isend=isend+1
              itag=(3**1+4**1)*(ineighbor-1)+(n_i1+1)*3**(1-1)
              call MPI_ISEND(pw(igrid)%wb,1,type_send_srl(iib1,i1),&
                  ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
           end if
        else
           ixSmin1=ixS_srl_min1(iib1,i1);ixSmax1=ixS_srl_max1(iib1,i1);
           select case (ipole)
           case (1)
              n_i1=i1;
           end select
           if (ipe_neighbor==mype) then
              ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmax1=ixR_srl_max1(iib1,n_i1);
              call pole_copy(pw(ineighbor)%wb,ixGmin1,ixGmax1,ixRmin1,ixRmax1,&
                 pw(igrid)%wb,ixGmin1,ixGmax1,ixSmin1,ixSmax1)
           else
              if (isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(:,&
                    isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%w)
              end if
              allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,nwhead:nwtail))
              call pole_buf(pwbuf(ipwbuf)%w,ixSmin1,ixSmax1,ixSmin1,ixSmax1,&
                 pw(igrid)%wb,ixGmin1,ixGmax1,ixSmin1,ixSmax1)
              isend=isend+1
              isend_buf(ipwbuf)=isend
              itag=(3**1+4**1)*(ineighbor-1)+(n_i1+1)*3**(1-1)
              isizes=(ixSmax1-ixSmin1+1)*nwbc
              call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                  ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
              ipwbuf=1+modulo(ipwbuf,npwbuf)
           end if
        end if

      end subroutine bc_send_srl

      !> Send to coarser neighbor
      subroutine bc_send_restrict
        integer :: ii1

        ic1=1+modulo(node(pig1_,igrid)-1,2);
        if (.not.(i1==0.or.i1==2*ic1-3)) return
        if(phyboundblock(igrid)) then
          ! filling physical boundary ghost cells of a coarser representative block for
          ! sending swap region with width of nghostcells to its coarser neighbor
          do idims=1,ndim
             ! to avoid using as yet unknown corner info in more than 1D, we
             ! fill only interior mesh ranges of the ghost cell ranges at first,
             ! and progressively enlarge the ranges to include corners later
             kmin1=merge(0, 1, idims==1)
             kmax1=merge(0, 1, idims==1)
             ixBmin1=ixCoGmin1+kmin1*nghostcells
             ixBmax1=ixCoGmax1-kmax1*nghostcells
             
             
             if(i1==-1) then
               ixBmin1=ixCoGmin1+nghostcells
               ixBmax1=ixCoGmin1+2*nghostcells-1
             else if(i1==1) then
               ixBmin1=ixCoGmax1-2*nghostcells+1
               ixBmax1=ixCoGmax1-nghostcells
             end if
             do iside=1,2
                ii1=kr(1,idims)*(2*iside-3);
                if (abs(i1)==1.and.abs(ii1)==1) cycle
                if (neighbor_type(ii1,igrid)/=neighbor_boundary) cycle
                call bc_phys(iside,idims,time,0.d0,pw(igrid)%wcoarse,&
                   pw(igrid)%xcoarse,ixCoGmin1,ixCoGmax1,ixBmin1,ixBmax1)
             end do
          end do
        end if

        ineighbor=neighbor(1,i1,igrid)
        ipe_neighbor=neighbor(2,i1,igrid)

        if (ipole==0) then
           n_inc1=-2*i1+ic1;
           if (ipe_neighbor==mype) then
              ixSmin1=ixS_r_min1(iib1,i1);ixSmax1=ixS_r_max1(iib1,i1);
              ixRmin1=ixR_r_min1(iib1,n_inc1);ixRmax1=ixR_r_max1(iib1,n_inc1);
              pw(ineighbor)%wb(ixRmin1:ixRmax1,&
                 nwhead:nwtail)=pw(igrid)%wcoarse(ixSmin1:ixSmax1,&
                 nwhead:nwtail)
           else
              isend=isend+1
              itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
              call MPI_ISEND(pw(igrid)%wcoarse,1,type_send_r(iib1,i1),&
                  ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
           end if
        else
           ixSmin1=ixS_r_min1(iib1,i1);ixSmax1=ixS_r_max1(iib1,i1);
           select case (ipole)
           case (1)
              n_inc1=2*i1+(3-ic1);
           end select
           if (ipe_neighbor==mype) then
              ixRmin1=ixR_r_min1(iib1,n_inc1);ixRmax1=ixR_r_max1(iib1,n_inc1);
              call pole_copy(pw(ineighbor)%wb,ixGmin1,ixGmax1,ixRmin1,ixRmax1,&
                 pw(igrid)%wcoarse,ixCoGmin1,ixCoGmax1,ixSmin1,ixSmax1)
           else
              if (isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(:,&
                    isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%w)
              end if
              allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,nwhead:nwtail))
              call pole_buf(pwbuf(ipwbuf)%w,ixSmin1,ixSmax1,ixSmin1,ixSmax1,&
                 pw(igrid)%wcoarse,ixCoGmin1,ixCoGmax1,ixSmin1,ixSmax1)
              isend=isend+1
              isend_buf(ipwbuf)=isend
              itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
              isizes=(ixSmax1-ixSmin1+1)*nwbc
              call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                  ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
              ipwbuf=1+modulo(ipwbuf,npwbuf)
           end if
        end if

      end subroutine bc_send_restrict

      !> Send to finer neighbor
      subroutine bc_send_prolong
        integer :: ii1

        do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
           inc1=2*i1+ic1
           ixSmin1=ixS_p_min1(iib1,inc1);ixSmax1=ixS_p_max1(iib1,inc1);

           ineighbor=neighbor_child(1,inc1,igrid)
           ipe_neighbor=neighbor_child(2,inc1,igrid)
        
           if (ipole==0) then
              n_i1=-i1;
              n_inc1=ic1+n_i1;
              if (ipe_neighbor==mype) then
                 ixRmin1=ixR_p_min1(iib1,n_inc1)
                 ixRmax1=ixR_p_max1(iib1,n_inc1);
                 pw(ineighbor)%wcoarse(ixRmin1:ixRmax1,&
                    nwhead:nwtail) =pw(igrid)%wb(ixSmin1:ixSmax1,&
                    nwhead:nwtail)
              else
                 isend=isend+1
                 itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
                 call MPI_ISEND(pw(igrid)%wb,1,type_send_p(iib1,inc1),&
                     ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
              end if
           else
              select case (ipole)
              case (1)
                 n_inc1=inc1;
              end select
              if (ipe_neighbor==mype) then
                 ixRmin1=ixR_p_min1(iib1,n_inc1)
                 ixRmax1=ixR_p_max1(iib1,n_inc1);
                 call pole_copy(pw(ineighbor)%wcoarse,ixCoGmin1,ixCoGmax1,&
                    ixRmin1,ixRmax1,pw(igrid)%wb,ixGmin1,ixGmax1,ixSmin1,&
                    ixSmax1)
              else
                 if (isend_buf(ipwbuf)/=0) then
                    call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(:,&
                       isend_buf(ipwbuf)),ierrmpi)
                    deallocate(pwbuf(ipwbuf)%w)
                 end if
                 allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,nwhead:nwtail))
                 call pole_buf(pwbuf(ipwbuf)%w,ixSmin1,ixSmax1,ixSmin1,ixSmax1,&
                    pw(igrid)%wb,ixGmin1,ixGmax1,ixSmin1,ixSmax1)
                 isend=isend+1
                 isend_buf(ipwbuf)=isend
                 itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
                 isizes=(ixSmax1-ixSmin1+1)*nwbc
                 call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                     ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
                 ipwbuf=1+modulo(ipwbuf,npwbuf)
              end if
           end if
        end do

      end subroutine bc_send_prolong

      !> Receive from sibling at same refinement level
      subroutine bc_recv_srl

        ipe_neighbor=neighbor(2,i1,igrid)
        if (ipe_neighbor/=mype) then
           irecv=irecv+1
           itag=(3**1+4**1)*(igrid-1)+(i1+1)*3**(1-1)
           call MPI_IRECV(pw(igrid)%wb,1,type_recv_srl(iib1,i1), ipe_neighbor,&
              itag,icomm,recvrequest(irecv),ierrmpi)
        end if

      end subroutine bc_recv_srl

      !> Receive from fine neighbor
      subroutine bc_recv_restrict

        do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
           inc1=2*i1+ic1
           ipe_neighbor=neighbor_child(2,inc1,igrid)
           if (ipe_neighbor/=mype) then
              irecv=irecv+1
              itag=(3**1+4**1)*(igrid-1)+3**1+inc1*4**(1-1)
              call MPI_IRECV(pw(igrid)%wb,1,type_recv_r(iib1,inc1),&
                  ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)
           end if
        end do

      end subroutine bc_recv_restrict

      !> Receive from coarse neighbor
      subroutine bc_recv_prolong

        ic1=1+modulo(node(pig1_,igrid)-1,2);
        if (.not.(i1==0.or.i1==2*ic1-3)) return

        ipe_neighbor=neighbor(2,i1,igrid)
        if (ipe_neighbor/=mype) then
           irecv=irecv+1
           inc1=ic1+i1;
           itag=(3**1+4**1)*(igrid-1)+3**1+inc1*4**(1-1)
           call MPI_IRECV(pw(igrid)%wcoarse,1,type_recv_p(iib1,inc1),&
               ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)  
        end if

      end subroutine bc_recv_prolong

      subroutine bc_prolong
        use mod_physics, only: phys_to_primitive, phys_to_conserved

        integer :: ixFimin1,ixFimax1,ixComin1,ixComax1,ii1
        double precision :: dxFi1, dxCo1, xFimin1, xComin1, invdxCo1

        ixFimin1=ixR_srl_min1(iib1,i1);ixFimax1=ixR_srl_max1(iib1,i1);
        dxFi1=rnode(rpdx1_,igrid);
        dxCo1=two*dxFi1;
        invdxCo1=1.d0/dxCo1;

        ! compute the enlarged grid lower left corner coordinates
        ! these are true coordinates for an equidistant grid, 
        ! but we can temporarily also use them for getting indices 
        ! in stretched grids
        xFimin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxFi1;
        xComin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxCo1;

        if(prolongprimitive) then
           ! following line again assumes equidistant grid, but 
           ! just computes indices, so also ok for stretched case
           ! TO CHECK: previously had +1-1 and +1+1 here: needed? Changed to +1 always
           ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-&
              xComin1)*invdxCo1)+1;
           ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-&
              xComin1)*invdxCo1)+1;
           call phys_to_primitive(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,&
              pw(igrid)%wcoarse,pw(igrid)%xcoarse)
        endif

        select case (typeghostfill)
        case ("linear")
           call interpolation_linear(ixFimin1,ixFimax1,dxFi1,xFimin1,dxCo1,&
              invdxCo1,xComin1)
        case ("copy")
           call interpolation_copy(ixFimin1,ixFimax1,dxFi1,xFimin1,dxCo1,&
              invdxCo1,xComin1)
        case ("unlimit")
           call interpolation_unlimit(ixFimin1,ixFimax1,dxFi1,xFimin1,dxCo1,&
              invdxCo1,xComin1)
        case default
           write (unitterm,*) "Undefined typeghostfill ",typeghostfill
           call mpistop("")
        end select

        if(prolongprimitive) call phys_to_conserved(ixCoGmin1,ixCoGmax1,&
           ixComin1,ixComax1,pw(igrid)%wcoarse,pw(igrid)%xcoarse)

      end subroutine bc_prolong

      subroutine interpolation_linear(ixFimin1,ixFimax1,dxFi1,xFimin1, dxCo1,&
         invdxCo1,xComin1)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: ixFimin1,ixFimax1
        double precision, intent(in) :: dxFi1, xFimin1,dxCo1, invdxCo1,&
            xComin1

        integer :: ixCo1, jxCo1, hxCo1, ixFi1, ix1, iw, idims, nwmin,nwmax
        double precision :: xCo1, xFi1, eta1
        double precision :: slopeL, slopeR, slopeC, signC, signR
        double precision :: slope(1:nw,ndim)
        double precision :: local_invdxCo1

        if(prolongprimitive) then
          nwmin=1
          nwmax=nw
        else
          nwmin=nwhead
          nwmax=nwtail
        end if

        do ixFi1 = ixFimin1,ixFimax1
           ! cell-centered coordinates of fine grid point
           ! here we temporarily use an equidistant grid
           xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1
        
           ! indices of coarse cell which contains the fine cell
           ! since we computed lower left corner earlier 
           ! in equidistant fashion: also ok for stretched case
           ixCo1=int((xFi1-xComin1)*invdxCo1)+1
        
           ! actual cell-centered coordinates of fine grid point
           xFi1=block%x(ixFi1,1)
           ! actual cell-centered coordinates of coarse grid point
           xCo1=block%xcoarse(ixCo1,1)

           if(.not.slab) then
              local_invdxCo1=1.d0/block%dxcoarse(ixCo1,1)
           endif

           ! normalized distance between fine/coarse cell center
           ! in coarse cell: ranges from -0.5 to 0.5 in each direction
           ! (origin is coarse cell center)
           if(slab) then
             eta1=(xFi1-xCo1)*invdxCo1;
           else
             ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;
             eta1=(xFi1-xCo1)*local_invdxCo1 *two*(one-block%dvolume(ixFi1) &
                /sum(block%dvolume(ix1:ix1+1))) 
           end if
        
           do idims=1,ndim
              hxCo1=ixCo1-kr(1,idims)
              jxCo1=ixCo1+kr(1,idims)
        
              do iw=nwmin,nwmax
                 slopeL=pw(igrid)%wcoarse(ixCo1,iw)-pw(igrid)%wcoarse(hxCo1,&
                    iw)
                 slopeR=pw(igrid)%wcoarse(jxCo1,iw)-pw(igrid)%wcoarse(ixCo1,&
                    iw)
                 slopeC=half*(slopeR+slopeL)
        
                 ! get limited slope
                 signR=sign(one,slopeR)
                 signC=sign(one,slopeC)
                 select case(typeprolonglimit)
                 case('minmod')
                   slope(iw,idims)=signR*max(zero,min(dabs(slopeR),&
                       signR*slopeL))
                 case('woodward')
                   slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR),&
                       signR*slopeL,signR*half*slopeC))
                 case('koren')
                   slope(iw,idims)=signR*max(zero,min(two*signR*slopeL,&
                       (dabs(slopeR)+two*slopeL*signR)*third,&
                      two*dabs(slopeR)))
                 case default
                   slope(iw,idims)=signC*max(zero,min(dabs(slopeC),&
                       signC*slopeL,signC*slopeR))
                 end select
              end do
           end do
        
           ! Interpolate from coarse cell using limited slopes
           pw(igrid)%wb(ixFi1,nwmin:nwmax)=pw(igrid)%wcoarse(ixCo1,&
              nwmin:nwmax)+(slope(nwmin:nwmax,1)*eta1)
        
        end do
        
        if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGhi1,ixFimin1,&
           ixFimax1,pw(igrid)%wb,pw(igrid)%x)
      
      end subroutine interpolation_linear

      subroutine interpolation_copy(ixFimin1,ixFimax1,dxFi1,xFimin1, dxCo1,&
         invdxCo1,xComin1)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: ixFimin1,ixFimax1
        double precision, intent(in) :: dxFi1, xFimin1,dxCo1, invdxCo1,&
            xComin1

        integer :: ixCo1, ixFi1, nwmin,nwmax
        double precision :: xFi1

        if(prolongprimitive) then
          nwmin=1
          nwmax=nw
        else
          nwmin=nwhead
          nwmax=nwtail
        end if

        do ixFi1 = ixFimin1,ixFimax1
           ! cell-centered coordinates of fine grid point
           xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1
        
           ! indices of coarse cell which contains the fine cell
           ! note: this also works for stretched grids
           ixCo1=int((xFi1-xComin1)*invdxCo1)+1
        
           ! Copy from coarse cell
           pw(igrid)%wb(ixFi1,nwmin:nwmax)=pw(igrid)%wcoarse(ixCo1,&
              nwmin:nwmax)
        
        end do
        
        if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGhi1,ixFimin1,&
           ixFimax1,pw(igrid)%wb,pw(igrid)%x)
      
      end subroutine interpolation_copy

      subroutine interpolation_unlimit(ixFimin1,ixFimax1,dxFi1,xFimin1, dxCo1,&
         invdxCo1,xComin1)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: ixFimin1,ixFimax1
        double precision, intent(in) :: dxFi1, xFimin1, dxCo1,invdxCo1,&
            xComin1

        integer :: ixCo1, jxCo1, hxCo1, ixFi1, ix1, idims, nwmin, nwmax
        double precision :: xCo1, xFi1, eta1
        double precision :: slope(1:nw,ndim)
        double precision :: local_invdxCo1

        if(prolongprimitive) then
          nwmin=1
          nwmax=nw
        else
          nwmin=nwhead
          nwmax=nwtail
        end if

        do ixFi1 = ixFimin1,ixFimax1
           ! cell-centered coordinates of fine grid point
           xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1
        
           ! indices of coarse cell which contains the fine cell
           ixCo1=int((xFi1-xComin1)*invdxCo1)+1

           ! actual cell-centered coordinates of fine grid point
           xFi1=block%x(ixFi1,1)
           ! actual cell-centered coordinates of coarse grid point
           xCo1=block%xcoarse(ixCo1,1)

           if(.not.slab) then
              local_invdxCo1=1.d0/block%dxcoarse(ixCo1,1)
           endif
        
           ! normalized distance between fine/coarse cell center
           ! in coarse cell: ranges from -0.5 to 0.5 in each direction
           ! (origin is coarse cell center)
           if (slab) then
              eta1=(xFi1-xCo1)*invdxCo1;
           else
              ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;
              eta1=(xFi1-xCo1)*local_invdxCo1 *two*(one-block%dvolume(ixFi1) &
                 /sum(block%dvolume(ix1:ix1+1))) 
           end if
        
           do idims=1,ndim
              hxCo1=ixCo1-kr(1,idims)
              jxCo1=ixCo1+kr(1,idims)
        
              ! get centered slope
              slope(nwmin:nwmax,idims)=half*(pw(igrid)%wcoarse(jxCo1,&
                 1:nw)-pw(igrid)%wcoarse(hxCo1,1:nw))
           end do
        
           ! Interpolate from coarse cell using centered slopes
           pw(igrid)%wb(ixFi1,nwmin:nwmax)=pw(igrid)%wcoarse(ixCo1,&
              nwmin:nwmax)+(slope(nwmin:nwmax,1)*eta1)
        end do
        
        if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGhi1,ixFimin1,&
           ixFimax1,pw(igrid)%wb,pw(igrid)%x)
      
      end subroutine interpolation_unlimit

      subroutine pole_copy(wrecv,ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,wsend,&
         ixISmin1,ixISmax1,ixSmin1,ixSmax1)
      
        integer, intent(in) :: ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,ixISmin1,&
           ixISmax1,ixSmin1,ixSmax1
        double precision :: wrecv(ixIRmin1:ixIRmax1,1:nw),&
            wsend(ixISmin1:ixISmax1,1:nw)

        integer :: iw, iB

        select case (ipole)
        case (1)
           iside=int((i1+3)/2)
           iB=2*(1-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case ("symm")
               wrecv(ixRmin1:ixRmax1,iw) = wsend(ixSmax1:ixSmin1:-1,iw)
             case ("asymm")
               wrecv(ixRmin1:ixRmax1,iw) =-wsend(ixSmax1:ixSmin1:-1,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do 
        end select
      
      end subroutine pole_copy

      subroutine pole_buf(wrecv,ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,wsend,&
         ixISmin1,ixISmax1,ixSmin1,ixSmax1)
      
        integer, intent(in) :: ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,ixISmin1,&
           ixISmax1,ixSmin1,ixSmax1
        double precision :: wrecv(ixIRmin1:ixIRmax1,nwhead:nwtail),&
            wsend(ixISmin1:ixISmax1,1:nw)

        integer :: iw, iB

        select case (ipole)
        case (1)
           iside=int((i1+3)/2)
           iB=2*(1-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case ("symm")
               wrecv(ixRmin1:ixRmax1,iw) = wsend(ixSmax1:ixSmin1:-1,iw)
             case ("asymm")
               wrecv(ixRmin1:ixRmax1,iw) =-wsend(ixSmax1:ixSmin1:-1,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do 
        end select
      
      end subroutine pole_buf

      subroutine fix_auxiliary
        use mod_physics, only: phys_get_aux
      
        integer :: ixmin1,ixmax1

        do iigrid=1,igridstail; igrid=igrids(iigrid);
          saveigrid=igrid
          block=>pw(igrid)
          call identifyphysbound(igrid,iib1)   
             
          do i1=-1,1
             if (skip_direction([ i1 ])) cycle
        
             ixmin1=ixR_srl_min1(iib1,i1);ixmax1=ixR_srl_max1(iib1,i1);
             call phys_get_aux(.true.,pw(igrid)%wb,pw(igrid)%x,ixGmin1,ixGmax1,&
                ixmin1,ixmax1,"bc")
          end do
        end do
      
      end subroutine fix_auxiliary

  end subroutine getbc

  subroutine identifyphysbound(igrid,iib1)
    use mod_global_parameters

    integer, intent(in)  :: igrid
    integer, intent(out) :: iib1

    
    if(pw(igrid)%is_physical_boundary(2*1) .and. &
       pw(igrid)%is_physical_boundary(2*1-1)) then
      iib1=2
    else if(pw(igrid)%is_physical_boundary(2*1-1)) then
      iib1=-1
    else if(pw(igrid)%is_physical_boundary(2*1)) then
      iib1=1
    else
      iib1=0
    end if
    

  end subroutine identifyphysbound

end module mod_ghostcells_update
