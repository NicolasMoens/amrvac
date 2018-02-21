!=============================================================================
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)

use mod_global_parameters

integer, intent(in) :: igrid, ipe
integer, dimension(2,2), intent(in) :: child_igrid, child_ipe
logical, intent(in) :: active

integer :: igridFi, ipeFi, ixComin1,ixComin2,ixComax1,ixComax2, ixCoGmin1,&
   ixCoGmin2,ixCoGmax1,ixCoGmax2, ixCoMmin1,ixCoMmin2,ixCoMmax1,ixCoMmax2, ic1,&
   ic2

if (ipe==mype) call alloc_node(igrid)

! New passive cell, coarsen from initial condition:
if (.not. active) then

   if (ipe == mype) then
      call initial_condition(igrid)
      
      do ic2=1,2
      do ic1=1,2
      igridFi=child_igrid(ic1,ic2)
      ipeFi=child_ipe(ic1,ic2)
      !if (ipeFi==mype) then
      !   ! remove solution space of child      
      !   call dealloc_node(igridFi)
      !end if
      end do
      end do
      
   end if

   return
end if



do ic2=1,2
do ic1=1,2
   igridFi=child_igrid(ic1,ic2)
   ipeFi=child_ipe(ic1,ic2)

   if (ipeFi==mype) then
      dxlevel(1)=rnode(rpdx1_,igridFi);dxlevel(2)=rnode(rpdx2_,igridFi);
      if (ipe==mype) then
         ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2
         ixComin2=ixMlo2+(ic2-1)*(ixMhi2-ixMlo2+1)/2;
         ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2
         ixComax2=ixMhi2+(ic2-2)*(ixMhi2-ixMlo2+1)/2;

         call coarsen_grid(pw(igridFi)%w,pw(igridFi)%x,ixGlo1,ixGlo2,ixGhi1,&
            ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,pw(igrid)%x,ixGlo1,&
            ixGlo2,ixGhi1,ixGhi2, ixComin1,ixComin2,ixComax1,ixComax2,igridFi,&
            igrid)

         ! remove solution space of child
         !call dealloc_node(igridFi)
      else
         ixCoGmin1=1;ixCoGmin2=1;
         ixCoGmax1=ixGhi1/2+nghostcells;ixCoGmax2=ixGhi2/2+nghostcells;
         ixCoMmin1=ixCoGmin1+nghostcells;ixCoMmin2=ixCoGmin2+nghostcells
         ixCoMmax1=ixCoGmax1-nghostcells;ixCoMmax2=ixCoGmax2-nghostcells;
         call coarsen_grid(pw(igridFi)%w,pw(igridFi)%x,ixGlo1,ixGlo2,ixGhi1,&
            ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igridFi)%wcoarse,&
            pw(igridFi)%xcoarse, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
            ixCoMmin1,ixCoMmin2,ixCoMmax1,ixCoMmax2,igridFi,igridFi)

         itag=ipeFi*max_blocks+igridFi
         isend=isend+1
         call MPI_ISEND(pw(igridFi)%wcoarse,1,type_coarse_block,ipe,itag,&
             icomm,sendrequest(isend),ierrmpi)
      end if
   else
      if (ipe==mype) then
         itag=ipeFi*max_blocks+igridFi
         irecv=irecv+1
         call MPI_IRECV(pw(igrid)%w,1,type_sub_block(ic1,ic2),ipeFi,itag,&
             icomm,recvrequest(irecv),ierrmpi)
      end if
   end if
end do
end do

end subroutine coarsen_grid_siblings
!=============================================================================
subroutine coarsen_grid(wFi,xFi,ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2,&
   ixFimin1,ixFimin2,ixFimax1,ixFimax2,wCo,xCo,ixCoGmin1,ixCoGmin2,ixCoGmax1,&
   ixCoGmax2,ixComin1,ixComin2,ixComax1,ixComax2,igridFi,igridCo)

use mod_global_parameters
use mod_physics, only: phys_to_primitive, phys_to_conserved

integer, intent(in) :: ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2, ixFimin1,&
   ixFimin2,ixFimax1,ixFimax2, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
    ixComin1,ixComin2,ixComax1,ixComax2, igridFi, igridCo
double precision, intent(inout) :: wFi(ixFiGmin1:ixFiGmax1,ixFiGmin2:ixFiGmax2,&
   1:nw), xFi(ixFiGmin1:ixFiGmax1,ixFiGmin2:ixFiGmax2,1:ndim)
double precision,intent(inout) :: wCo(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,&
   1:nw), xCo(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:ndim)

integer :: ixCo1,ixCo2, ixFi1,ixFi2, iw
double precision :: CoFiratio
!-----------------------------------------------------------------------------
! coarsen by 2 in every direction - conservatively

if(coarsenprimitive) call phys_to_primitive(ixFiGmin1,ixFiGmin2,ixFiGmax1,&
   ixFiGmax2,ixFimin1,ixFimin2,ixFimax1,ixFimax2,wFi,xFi)

if (slab) then
   CoFiratio=one/dble(2**ndim)
   do iw=1,nw
      do ixCo2 = ixComin2,ixComax2
         ixFi2=2*(ixCo2-ixComin2)+ixFimin2
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,ixCo2,iw)=sum(wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
            iw))*CoFiratio
      end do
      end do
   end do
else
   if(igridFi==igridCo) then
     do iw=1,nw
        do ixCo2 = ixComin2,ixComax2
           ixFi2=2*(ixCo2-ixComin2)+ixFimin2
        do ixCo1 = ixComin1,ixComax1
           ixFi1=2*(ixCo1-ixComin1)+ixFimin1
           wCo(ixCo1,ixCo2,iw)= sum(pw(igridFi)%dvolume(ixFi1:ixFi1+1,&
              ixFi2:ixFi2+1)*wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
              iw)) /pw(igridCo)%dvolumecoarse(ixCo1,ixCo2)
        end do
        end do
     end do
   else
     do iw=1,nw
        do ixCo2 = ixComin2,ixComax2
           ixFi2=2*(ixCo2-ixComin2)+ixFimin2
        do ixCo1 = ixComin1,ixComax1
           ixFi1=2*(ixCo1-ixComin1)+ixFimin1
           wCo(ixCo1,ixCo2,iw)= sum(pw(igridFi)%dvolume(ixFi1:ixFi1+1,&
              ixFi2:ixFi2+1)*wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
              iw)) /pw(igridCo)%dvolume(ixCo1,ixCo2)
        end do
        end do
     end do
   end if
end if

if(coarsenprimitive) then
  call phys_to_conserved(ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2,ixFimin1,&
     ixFimin2,ixFimax1,ixFimax2,wFi,xFi)
  call phys_to_conserved(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixComin1,&
     ixComin2,ixComax1,ixComax2,wCo,xCo)
end if

end subroutine coarsen_grid
!=============================================================================
