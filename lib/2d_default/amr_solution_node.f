!> Get first available igrid on processor ipe
integer function getnode(ipe)
use mod_forest, only: igrid_inuse
use mod_global_parameters

integer, intent(in) :: ipe
integer :: igrid, igrid_available

igrid_available=0

do igrid=1,max_blocks
   if (igrid_inuse(igrid,ipe)) cycle

   igrid_available=igrid
   exit
end do

if (igrid_available == 0) then
   getnode = -1
   print *, "Current maximum number of grid blocks:", max_blocks
   call mpistop("Insufficient grid blocks; increase max_blocks in meshlist")
else
   getnode=igrid_available
   igrid_inuse(igrid,ipe)=.true.
end if

if (ipe==mype) then
   ! initialize nodal block
   node(1:nodehi,getnode) = 0
   rnode(1:rnodehi,getnode) = zero
end if

end function getnode
!=============================================================================
subroutine putnode(igrid,ipe)
use mod_forest
implicit none

! putnode = return igrid node on processor ipe
 
integer, intent(in) :: igrid, ipe
!----------------------------------------------------------------------------
igrid_inuse(igrid,ipe)=.false.

end subroutine putnode
!=============================================================================
subroutine alloc_node(igrid)
use mod_forest
use mod_global_parameters

integer, intent(in) :: igrid

integer :: level, ig1,ig2, ign1,ign2, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
    ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmax1,ixCoCoGmax2, ix, i1,i2
integer :: imin, imax, index, ig1Co, ig2Co, ixshift
double precision :: rXmin1,rXmin2, dx1,dx2
logical, save:: first=.true.
!-----------------------------------------------------------------------------
ixCoGmin1=1;ixCoGmin2=1;
!ixCoGmax^D=ixGhi^D/2+nghostcells;
ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells;

if(.not. allocated(pw(igrid)%w)) then
  
  ! initialize solution space
  allocate(pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
      pw(igrid)%wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
      pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
      pw(igrid)%wcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nw))
  
  ! wio for visualization data
  allocate(pw(igrid)%wio(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwauxio))
  
  ! allocate temperary solution space
  select case (time_integrator)
  case("threestep","fourstep","jameson","twostep_trapezoidal")
    allocate(pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
  case("rk4","ssprk43")
    allocate(pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
    allocate(pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
  case("ssprk54")
    allocate(pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
    allocate(pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
    allocate(pw(igrid)%w4(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw))
  end select
  
  ! allocate coordinates
  allocate(pw(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),&
      pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:ndim))
end if

pw(igrid)%w=0.d0
pw(igrid)%wold=0.d0
pw(igrid)%w1=0.d0
pw(igrid)%wcoarse=0.d0
pw(igrid)%igrid=igrid
pw(igrid)%wio=0.d0

! wb is w by default
pw(igrid)%wb=>pw(igrid)%w


! block pointer to current block
block=>pw(igrid)

if(phys_energy .and. solve_internal_e) then
  block%e_is_internal=.true.
endif

! set level information
level=igrid_to_node(igrid,mype)%node%level
ig1=igrid_to_node(igrid,mype)%node%ig1;ig2=igrid_to_node(igrid,mype)%node%ig2;

node(plevel_,igrid)=level
node(pig1_,igrid)=ig1
node(pig2_,igrid)=ig2

! set dx information
rnode(rpdx1_,igrid)=dx(1,level)
rnode(rpdx2_,igrid)=dx(2,level)
dxlevel(:)=dx(:,level)

! uniform cartesian case as well as all unstretched coordinates
! determine the minimal and maximal corners
rnode(rpxmin1_,igrid)=xprobmin1+dble(ig1-1)*dg1(level)
rnode(rpxmin2_,igrid)=xprobmin2+dble(ig2-1)*dg2(level)
rnode(rpxmax1_,igrid)=xprobmax1-dble(ng1(level)-ig1)*dg1(level)
rnode(rpxmax2_,igrid)=xprobmax2-dble(ng2(level)-ig2)*dg2(level)

dx1=rnode(rpdx1_,igrid)
dx2=rnode(rpdx2_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-nghostcells*dx1
rXmin2=rnode(rpxmin2_,igrid)-nghostcells*dx2
do ix=ixGlo1,ixGhi1
   pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixGlo2,ixGhi2
   pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=rXmin2+(dble(ix)-half)*dx2
end do
dx1=2.0d0*rnode(rpdx1_,igrid)
dx2=2.0d0*rnode(rpdx2_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-nghostcells*dx1
rXmin2=rnode(rpxmin2_,igrid)-nghostcells*dx2
do ix=ixCoGmin1,ixCoGmax1
   pw(igrid)%xcoarse(ix,ixCoGmin2:ixCoGmax2,1)=rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixCoGmin2,ixCoGmax2
   pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ix,2)=rXmin2+(dble(ix)-half)*dx2
end do

if(stretched_grid) then
  if(slab_stretched)then
    ! cartesian and last coordinate stretching
    imin=(ig2-1)*block_nx2+1
    imax=ig2*block_nx2
    rnode(rpxmin2_,igrid)=xprobmin2+dxfirst_1mq(level) &
       *(1.0d0-qstretch(level)**(imin-1))
    rnode(rpxmax2_,igrid)=xprobmin2+dxfirst_1mq(level) &
       *(1.0d0-qstretch(level)**imax)
    ! fix possible out of bound due to precision
    if(rnode(rpxmax2_,igrid)>xprobmax2) then
       if(first) then
          write(*,*) 'Warning: edge beyond domain?',igrid,imax,rnode(rpxmax2_,&
             igrid)
          first=.false.
       endif
       rnode(rpxmax2_,igrid)=xprobmax2
    endif
    ixshift=(ig2-1)*block_nx2-nghostcells
    do ix=ixGlo2,ixGhi2
      index=ixshift+ix
      pw(igrid)%x(ixGlo1:ixGhi1,ix,2)=xprobmin2+&
         dxfirst_1mq(level)*(1.0d0-qstretch(level)**(index-1)) + &
         0.5d0*dxfirst(level)*qstretch(level)**(index-1)
    enddo
    ig2Co=(ig2-1)/2
    ixshift=ig2Co*block_nx2+(1-mod(ig2,2))*block_nx2/2-nghostcells
    do ix=ixCoGmin2,ixCoGmax2
      index=ixshift+ix
      pw(igrid)%xcoarse(ixCoGmin1:ixCoGmax1,ix,&
         2)=xprobmin2+dxfirst_1mq(level-1)*(1.0d0-qstretch(level-1)**(index-&
         1)) + 0.5d0*dxfirst(level)*qstretch(level-1)**(index-1)
    end do
  else
    ! cylindrical/spherical and radial stretching
    imin=(ig1-1)*block_nx1+1
    imax=ig1*block_nx1
    rnode(rpxmin1_,igrid)=xprobmin1+dxfirst_1mq(level) &
       *(1.0d0-qstretch(level)**(imin-1))
    rnode(rpxmax1_,igrid)=xprobmin1+dxfirst_1mq(level) &
       *(1.0d0-qstretch(level)**imax)
    ! fix possible out of bound due to precision
    if(rnode(rpxmax1_,igrid)>xprobmax1) then
       if(first) then
          write(*,*) 'Warning: edge beyond domain?',igrid,imax,rnode(rpxmax1_,&
             igrid)
          first=.false.
       endif
       rnode(rpxmax1_,igrid)=xprobmax1
    endif
    ixshift=(ig1-1)*block_nx1-nghostcells
    do ix=ixGlo1,ixGhi1
      index=ixshift+ix
      pw(igrid)%x(ix,ixGlo2:ixGhi2,1)=xprobmin1+&
         dxfirst_1mq(level)*(1.0d0-qstretch(level)**(index-1)) + &
         0.5d0*dxfirst(level)*qstretch(level)**(index-1)
    enddo
    ig1Co=(ig1-1)/2
    ixshift=ig1Co*block_nx1+(1-mod(ig1,2))*block_nx1/2-nghostcells
    do ix=ixCoGmin1,ixCoGmax1
      index=ixshift+ix
      pw(igrid)%xcoarse(ix,ixCoGmin2:ixCoGmax2,&
         1)=xprobmin1+dxfirst_1mq(level-1)*(1.0d0-qstretch(level-1)**(index-&
         1)) + 0.5d0*dxfirst(level)*qstretch(level-1)**(index-1)
    end do
  endif
endif

if (.not.slab) call getgridgeo(igrid)

if (B0field) then
   ! initialize background non-evolving solution
   call alloc_B0_grid(igrid)
   call set_B0_grid(igrid)
end if

! find the blocks on the boundaries
pw(igrid)%is_physical_boundary=.false.

do i1=-1,1
  if(i1==0) cycle
  ign1=ig1+i1
  ! blocks at periodic boundary have neighbors in the physical domain
  ! thus threated at internal blocks with no physical boundary
  if (periodB(1)) ign1=1+modulo(ign1-1,ng1(level))
  if (ign1 > ng1(level)) then
     if(phi_ > 0 .and. poleB(2,1)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*1)=.false.
     else
       pw(igrid)%is_physical_boundary(2*1)=.true.
     end if
  else if (ign1 < 1) then
     if(phi_ > 0 .and. poleB(1,1)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*1-1)=.false.
     else
       pw(igrid)%is_physical_boundary(2*1-1)=.true.
     end if
  end if
end do


do i2=-1,1
  if(i2==0) cycle
  ign2=ig2+i2
  ! blocks at periodic boundary have neighbors in the physical domain
  ! thus threated at internal blocks with no physical boundary
  if (periodB(2)) ign2=1+modulo(ign2-1,ng2(level))
  if (ign2 > ng2(level)) then
     if(phi_ > 0 .and. poleB(2,2)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*2)=.false.
     else
       pw(igrid)%is_physical_boundary(2*2)=.true.
     end if
  else if (ign2 < 1) then
     if(phi_ > 0 .and. poleB(1,2)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*2-1)=.false.
     else
       pw(igrid)%is_physical_boundary(2*2-1)=.true.
     end if
  end if
end do

if(any(pw(igrid)%is_physical_boundary)) then
  phyboundblock(igrid)=.true.
else
  phyboundblock(igrid)=.false.
end if

end subroutine alloc_node
!=============================================================================
subroutine dealloc_node(igrid)

use mod_global_parameters

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (igrid==0) then
   call mpistop("trying to delete a non-existing grid in dealloc_node")
end if

deallocate(pw(igrid)%w,pw(igrid)%w1,pw(igrid)%x)
deallocate(pw(igrid)%wcoarse,pw(igrid)%xcoarse)
deallocate(pw(igrid)%wold,pw(igrid)%wio)
! deallocate temperary solution space
select case (time_integrator)
case("threestep","fourstep","jameson","twostep_trapezoidal")
  deallocate(pw(igrid)%w2)
case("rk4","ssprk43")
  deallocate(pw(igrid)%w2)
  deallocate(pw(igrid)%w3)
case("ssprk54")
  deallocate(pw(igrid)%w2)
  deallocate(pw(igrid)%w3)
  deallocate(pw(igrid)%w4)
end select
if(allocated(pw(igrid)%w2)) deallocate(pw(igrid)%w2)
if(allocated(pw(igrid)%w3)) deallocate(pw(igrid)%w3)

if (.not.slab) call putgridgeo(igrid)

if (B0field) call dealloc_B0_grid(igrid)

end subroutine dealloc_node
!=============================================================================
