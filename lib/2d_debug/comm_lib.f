!> \file
!> @todo Convert this file to a module

!> Initialize the MPI environment
!> @todo Check for errors in return code
subroutine comm_start
  use mod_global_parameters

  integer(kind=MPI_ADDRESS_KIND) :: lb
  integer(kind=MPI_ADDRESS_KIND) :: size

  ! Initialize MPI
  call MPI_INIT(ierrmpi)

  ! Each process stores its rank, which ranges from 0 to N-1, where N is the
  ! number of processes.
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)

  ! Store the number of processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)

  ! Use the default communicator, which contains all the processes
  icomm = MPI_COMM_WORLD

  ! Get size of double/integer
  call MPI_TYPE_GET_EXTENT(MPI_REAL,lb,size,ierrmpi)
  if (size /= size_real) call mpistop("Incompatible real size")
  call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size,ierrmpi)
  if (size /= size_double) call mpistop("Incompatible double size")
  call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,size,ierrmpi)
  if (size /= size_int) call mpistop("Incompatible integer size")
  call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,size,ierrmpi)
  if (size /= size_logical) call mpistop("Incompatible logical size")

end subroutine comm_start


!> Finalize (or shutdown) the MPI environment
subroutine comm_finalize

  use mod_global_parameters
  use mod_ghostcells_update

  call put_bc_comm_types
  call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
  call MPI_FINALIZE(ierrmpi)

end subroutine comm_finalize


!> Create and store the MPI types that will be used for parallel communication
subroutine init_comm_types

use mod_global_parameters

integer, dimension(ndim+1) :: sizes, subsizes, start
integer :: i1,i2, ic1,ic2, nx1,nx2, nxCo1,nxCo2, nxG1,nxG2
!-----------------------------------------------------------------------------
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
nxG1=ixGhi1-ixGlo1+1;nxG2=ixGhi2-ixGlo2+1;
nxCo1=nx1/2;nxCo2=nx2/2;

sizes(1)=ixGhi1;sizes(2)=ixGhi2;
sizes(ndim+1)=nw
subsizes(1)=nxG1;subsizes(2)=nxG2;
subsizes(ndim+1)=nw
start(1)=ixGlo1-1;start(2)=ixGlo2-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block,ierrmpi)
call MPI_TYPE_COMMIT(type_block,ierrmpi)
size_block=nxG1*nxG2*nw*size_double

sizes(1)=ixGhi1/2+nghostcells;sizes(2)=ixGhi2/2+nghostcells;
sizes(ndim+1)=nw
subsizes(1)=nxCo1;subsizes(2)=nxCo2;
subsizes(ndim+1)=nw
start(1)=ixMlo1-1;start(2)=ixMlo2-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_coarse_block,ierrmpi)
call MPI_TYPE_COMMIT(type_coarse_block,ierrmpi)

sizes(1)=ixGhi1;sizes(2)=ixGhi2;
sizes(ndim+1)=nw
do ic2=1,2
do ic1=1,2
   subsizes(1)=nxCo1;subsizes(2)=nxCo2;
   subsizes(ndim+1)=nw
   start(1)=ixMlo1-1+(ic1-1)*nxCo1;start(2)=ixMlo2-1+(ic2-1)*nxCo2;
   start(ndim+1)=0
   call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
       MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_sub_block(ic1,ic2),&
      ierrmpi)
   call MPI_TYPE_COMMIT(type_sub_block(ic1,ic2),ierrmpi)
end do
end do

sizes(1)=ixGhi1;sizes(2)=ixGhi2;
sizes(ndim+1)=nw
subsizes(1)=nx1;subsizes(2)=nx2;
subsizes(ndim+1)=nw
start(1)=ixMlo1-1;start(2)=ixMlo2-1;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_io,ierrmpi)
size_block_io=nx1*nx2*nw*size_double

sizes(1)=ixMhi1-ixMlo1+1;sizes(2)=ixMhi2-ixMlo2+1;
sizes(ndim+1)=2
subsizes(1)=sizes(1);subsizes(2)=sizes(2);
subsizes(ndim+1)=2
start(1)=0;start(2)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_xcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xcc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+2;sizes(2)=ixMhi2-ixMlo2+2;
sizes(ndim+1)=2
subsizes(1)=sizes(1);subsizes(2)=sizes(2);
subsizes(ndim+1)=2
start(1)=0;start(2)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_xc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_xc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+1;sizes(2)=ixMhi2-ixMlo2+1;
sizes(ndim+1)=nw+nwauxio
subsizes(1)=sizes(1);subsizes(2)=sizes(2);
subsizes(ndim+1)=nw+nwauxio
start(1)=0;start(2)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_wcc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wcc_io,ierrmpi)

sizes(1)=ixMhi1-ixMlo1+2;sizes(2)=ixMhi2-ixMlo2+2;
sizes(ndim+1)=nw+nwauxio
subsizes(1)=sizes(1);subsizes(2)=sizes(2);
subsizes(ndim+1)=nw+nwauxio
start(1)=0;start(2)=0;
start(ndim+1)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_block_wc_io,ierrmpi)
call MPI_TYPE_COMMIT(type_block_wc_io,ierrmpi)

end subroutine init_comm_types


!> Exit MPI-AMRVAC with an error message
subroutine mpistop(message)
  use mod_global_parameters

  character(len=*), intent(in) :: message !< The error message
  integer                      :: ierrcode

  write(*, *) "ERROR for processor", mype, ":"
  write(*, *) trim(message)

  call MPI_ABORT(icomm, ierrcode, ierrmpi)

end subroutine mpistop
