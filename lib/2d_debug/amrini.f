!=============================================================================
!> Generate and initialize all grids at the coarsest level (level one)
subroutine initlevelone


use mod_global_parameters

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
levmin=1
levmax=1

call init_forest_root

call getigrids
call build_connectivity

! fill solution space of all root grids
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   call alloc_node(igrid)
   ! in case gradient routine used in initial condition, ensure geometry known
   call initial_condition(igrid)
end do


end subroutine initlevelone
!=============================================================================
subroutine initial_condition(igrid)

! Need only to set the mesh values (can leave ghost cells untouched)
use mod_usr_methods, only: usr_init_one_grid
use mod_global_parameters

integer, intent(in) :: igrid

!----------------------------------------------------------------------------
pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=zero

saveigrid=igrid
! in case gradient routine used in initial condition, ensure geometry known
block=>pw(igrid)
dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
typelimiter=type_limiter(node(plevel_,igrid))
typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

if (.not. associated(usr_init_one_grid)) then
   call mpistop("usr_init_one_grid not defined")
else
   call usr_init_one_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
      ixMhi2,pw(igrid)%w,pw(igrid)%x)
end if

end subroutine initial_condition
!=============================================================================
subroutine modify_IC
use mod_usr_methods, only: usr_init_one_grid
use mod_global_parameters

integer :: iigrid, igrid

!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   block=>pw(igrid)
   dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
   typelimiter=type_limiter(node(plevel_,igrid))
   typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

   if (.not. associated(usr_init_one_grid)) then
      call mpistop("usr_init_one_grid not defined")
   else
      call usr_init_one_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
         ixMhi2,pw(igrid)%w,pw(igrid)%x)
   end if
end do

end subroutine modify_IC
!=============================================================================
