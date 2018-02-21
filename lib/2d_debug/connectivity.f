!=============================================================================
subroutine get_level_range
use mod_forest
use mod_global_parameters

integer :: level
!-----------------------------------------------------------------------------

! determine new finest level
do level=refine_max_level,1,-1
   if (associated(level_tail(level)%node)) then
      levmax=level
      exit
   end if
end do

! determine coarsest level
do level=1,levmax
   if (associated(level_tail(level)%node)) then
      levmin=level
      exit
   end if
end do

end subroutine get_level_range
!=============================================================================
subroutine getigrids
use mod_forest
use mod_global_parameters
implicit none

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
iigrid=0
do igrid=1,max_blocks
   if (igrid_inuse(igrid,mype)) then
      iigrid=iigrid+1
      igrids(iigrid)=igrid
   end if
end do

igridstail=iigrid

end subroutine getigrids
!=============================================================================
subroutine build_connectivity
use mod_forest
use mod_global_parameters

integer :: iigrid, igrid, i1,i2, my_neighbor_type
integer :: iside, idim, ic1,ic2, inc1,inc2, ih1,ih2, icdim
type(tree_node_ptr) :: tree, my_neighbor, child
logical, dimension(2) :: pole
logical :: nopole
!-----------------------------------------------------------------------------
nrecv_bc_srl=0; nsend_bc_srl=0
nrecv_bc_r=0; nsend_bc_r=0
nrecv_bc_p=0; nsend_bc_p=0
nrecv_fc=0; nsend_fc=0

do iigrid=1,igridstail; igrid=igrids(iigrid);
   tree%node => igrid_to_node(igrid,mype)%node

   do i2=-1,1
   do i1=-1,1
      ! skip the grid itself
      if (i1==0.and.i2==0) then
         neighbor_type(0,0,igrid)=0
         neighbor(1,0,0,igrid)=igrid
         neighbor(2,0,0,igrid)=mype
      else
         call find_neighbor(my_neighbor,my_neighbor_type,tree,i1,i2,pole)
         nopole=.not.any(pole)

         select case (my_neighbor_type)
         ! adjacent to physical boundary
         case (neighbor_boundary)
            neighbor(1,i1,i2,igrid)=0
            neighbor(2,i1,i2,igrid)=-1
         ! fine-coarse transition
         case (neighbor_coarse)
            neighbor(1,i1,i2,igrid)=my_neighbor%node%igrid
            neighbor(2,i1,i2,igrid)=my_neighbor%node%ipe
            if (my_neighbor%node%ipe/=mype) then
               ic1=1+modulo(tree%node%ig1-1,2)
               ic2=1+modulo(tree%node%ig2-1,2);
               if ((i1==0.or.i1==2*ic1-3).and.(i2==0.or.i2==2*ic2-3)) then
                  nrecv_bc_p=nrecv_bc_p+1
                  nsend_bc_r=nsend_bc_r+1
               end if
            end if
         ! same refinement level
         case (neighbor_sibling)
            neighbor(1,i1,i2,igrid)=my_neighbor%node%igrid
            neighbor(2,i1,i2,igrid)=my_neighbor%node%ipe
            if (my_neighbor%node%ipe/=mype) then
               nrecv_bc_srl=nrecv_bc_srl+1
               nsend_bc_srl=nsend_bc_srl+1
            end if
         ! coarse-fine transition
         case (neighbor_fine)
            neighbor(1,i1,i2,igrid)=0
            neighbor(2,i1,i2,igrid)=-1
            do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
               if (pole(2)) then
                  ih2=3-ic2
               else
                  ih2=ic2
               end if
            do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               if (pole(1)) then
                  ih1=3-ic1
               else
                  ih1=ic1
               end if
               child%node => my_neighbor%node%child(ih1,ih2)%node
               neighbor_child(1,inc1,inc2,igrid)=child%node%igrid
               neighbor_child(2,inc1,inc2,igrid)=child%node%ipe
               if (child%node%ipe/=mype) then
                  nrecv_bc_r=nrecv_bc_r+1
                  nsend_bc_p=nsend_bc_p+1
               end if
            end do
            end do
         end select

         ! flux fix for conservation only for pure directional shifts
         if (abs(i1)+abs(i2)==1) then
            if (i1/=0) then
               idim=1
               iside=int((i1+3)/2)
            end if
            if (i2/=0) then
               idim=2
               iside=int((i2+3)/2)
            end if
            select case (my_neighbor_type)
            ! only across fine-coarse or coarse-fine boundaries
            case (neighbor_coarse)
               if (my_neighbor%node%ipe/=mype) then
                  if (.not.pole(idim)) nsend_fc(idim)=nsend_fc(idim)+1
               end if
            case (neighbor_fine)
               if (pole(idim)) then
                  icdim=iside
               else
                  icdim=3-iside
               end if
               select case (idim)
               case (1)
                  do ic1=icdim,icdim
               do ic2=1,2
                     child%node => my_neighbor%node%child(ic1,ic2)%node
                     if (child%node%ipe/=mype) then
                        if (.not.pole(1)) nrecv_fc(1)=nrecv_fc(1)+1
                     end if
                  end do
               end do 
               case (2)
                  do ic1=1,2
               do ic2=icdim,icdim
                     child%node => my_neighbor%node%child(ic1,ic2)%node
                     if (child%node%ipe/=mype) then
                        if (.not.pole(2)) nrecv_fc(2)=nrecv_fc(2)+1
                     end if
                  end do
               end do 
               end select
            end select
         end if

         if (phi_ > 0) then
           neighbor_pole(i1,i2,igrid)=0
           if (my_neighbor_type>1) then
             do idim=1,2
               if (pole(idim)) then
                 neighbor_pole(i1,i2,igrid)=idim
                 exit ! there can only be one pole between two meshes
               end if
             end do
           end if
         end if
         neighbor_type(i1,i2,igrid)=my_neighbor_type

      end if
   end do
   end do
end do

end subroutine build_connectivity
!=============================================================================
