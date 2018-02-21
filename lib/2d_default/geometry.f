subroutine set_coordinate_system(geom)
  use mod_global_parameters

  character(len=*), intent(in) :: geom

  select case (geom)
  case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
    ndir = ndim
    typeaxial='slab'
  case ("Cartesian_1.5D")
    if (ndim /= 1) call mpistop("Geometry Cartesian_1.5D but ndim /= 1")
    typeaxial='slab'
    ndir = 2
  case ("Cartesian_1.75D")
    if (ndim /= 1) call mpistop("Geometry Cartesian_1.75D but ndim /= 1")
    typeaxial='slab'
    ndir = 3
  case ("Cartesian_2.5D")
    if (ndim /= 2) call mpistop("Geometry Cartesian_2.5D but ndim /= 2")
    typeaxial='slab'
    ndir = 3
  case ("cylindrical","cylindrical_2D","cylindrical_3D")
    ndir = ndim
    r_   = 1
    z_   = 2
    if(ndir==3) phi_ = 3
    typeaxial='cylindrical'
  case ("cylindrical_2.5D")
    if (ndim /= 2) call mpistop("Geometry cylindrical_2.5D but ndim /= 2")
    ndir = 3
    r_   = 1
    z_   = 2
    phi_ = 3
    typeaxial='cylindrical'
  case ("polar","polar_2D","polar_3D")
    ndir = ndim
    r_   = 1
    phi_ = 2
    if(ndir==3) z_ = 3
    typeaxial='cylindrical'
  case ("polar_2.5D")
    if (ndim /= 2) call mpistop("Geometry polar_2.5D but ndim /= 2")
    ndir = 3
    r_   = 1
    phi_ = 2
    z_   = 3
    typeaxial='cylindrical'
  case ("spherical","spherical_2D","spherical_3D")
    ndir = ndim
    r_   = 1
    if(ndir==3) phi_ = 3
    z_   = -1
    typeaxial='spherical'
  case ("spherical_2.5D")
    if (ndim /= 2) call mpistop("Geometry spherical_2.5D requires ndim == 2")
    ndir = 3
    r_   = 1
    phi_ = 3
    z_   = -1
    typeaxial='spherical'
  case default
    call mpistop("Unknown geometry specified")
  end select
end subroutine set_coordinate_system

subroutine set_pole

  use mod_global_parameters
  
  select case (typeaxial)
  case ("spherical") 
  case ("cylindrical")
    if (1 == phi_ .and. periodB(1)) then
      if(mod(ng1(1),2)/=0) then
        call mpistop("Number of meshes in phi-direction should be even!")
      end if

      if(abs(xprobmin1)<smalldouble) then
        if (mype==0) then
          write(unitterm,*) "Will apply pi-periodic conditions at r=0"
        end if
        poleB(1,1)=.true.
      else
        if (mype==0) then
          write(unitterm,*) "There is no cylindrical axis!"
        end if
      end if
    end if
    if (2 == phi_ .and. periodB(2)) then
      if(mod(ng2(1),2)/=0) then
        call mpistop("Number of meshes in phi-direction should be even!")
      end if

      if(abs(xprobmin1)<smalldouble) then
        if (mype==0) then
          write(unitterm,*) "Will apply pi-periodic conditions at r=0"
        end if
        poleB(1,1)=.true.
      else
        if (mype==0) then
          write(unitterm,*) "There is no cylindrical axis!"
        end if
      end if
    end if
  end select

end subroutine set_pole
!=============================================================================
subroutine getgridgeo(igrid)

use mod_global_parameters

integer, intent(in) :: igrid

integer :: ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
double precision :: xmin1,xmin2, dx1,dx2
!-----------------------------------------------------------------------------

ixCoGmin1=1;ixCoGmin2=1; 
!ixCoGmax^D=ixGhi^D/2+nghostcells;
ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells;

if(.not. allocated(pw(igrid)%surfaceC)) then
  ! allocate geometric info
  allocate(pw(igrid)%surfaceC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),&
      pw(igrid)%surface(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),&
      pw(igrid)%dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
      pw(igrid)%dx(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),&
     pw(igrid)%dxcoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:ndim),&
      pw(igrid)%dvolumecoarse(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2))
end if

dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);
xmin1=rnode(rpxmin1_,igrid);xmin2=rnode(rpxmin2_,igrid);

! first fill the grid itself
pw(igrid)%dvolumep=>pw(igrid)%dvolume
pw(igrid)%dxp=>pw(igrid)%dx
pw(igrid)%xCCp=>pw(igrid)%x
call fillgeo(igrid,ixGlo1,ixGlo2,ixGhi1,ixGhi2,xmin1,xmin2,dx1,dx2,.false.)

! then fill its coarse representation
dx1=2.d0*rnode(rpdx1_,igrid);dx2=2.d0*rnode(rpdx2_,igrid);
pw(igrid)%dvolumep=>pw(igrid)%dvolumecoarse
pw(igrid)%dxp=>pw(igrid)%dxcoarse
pw(igrid)%xCCp=>pw(igrid)%xcoarse
call fillgeo(igrid,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,xmin1,xmin2,dx1,dx2,&
   .true.)

end subroutine getgridgeo
!=============================================================================
subroutine putgridgeo(igrid)

  use mod_global_parameters

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
deallocate(pw(igrid)%surfaceC,pw(igrid)%surface,pw(igrid)%dvolume,pw(igrid)%dx,&
   pw(igrid)%dxcoarse,pw(igrid)%dvolumecoarse)

end subroutine putgridgeo
!=============================================================================
subroutine fillgeo(igrid,ixGmin1,ixGmin2,ixGmax1,ixGmax2,xmin1,xmin2,dx1,dx2,&
   need_only_volume)

use mod_global_parameters

integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmax1,ixGmax2
double precision, intent(in) :: xmin1,xmin2, dx1,dx2
logical, intent(in) :: need_only_volume

integer :: idims, ix, ixmin1,ixmin2,ixmax1,ixmax2, ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2
integer :: level, ig1,ig2, ig1Co, index, ixshift, ig2Co
double precision :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ndim),&
    drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
!-----------------------------------------------------------------------------
ixmin1=ixGmin1+1;ixmin2=ixGmin2+1;ixmax1=ixGmax1-1;ixmax2=ixGmax2-1;

select case (typeaxial)
case ("slabstretch")
   ! this case stretches the final ndim direction in a cartesian grid
  
   ! grid level and global index for filling the grid itself
   level=node(plevel_,igrid)
   ig1=node(pig1_,igrid)
   ig2=node(pig2_,igrid)
   ! when filling the coarse grid (through need_only_volume) then adjust
   if(need_only_volume)then
      level=level-1
      ig2Co=(ig2-1)/2
      ixshift=ig2Co*block_nx2+(1-mod(ig2,2))*block_nx2/2-nghostcells
   else
      ixshift=(ig2-1)*block_nx2-nghostcells
   endif
   do ix = ixGmin2,ixGmax2
      index=ixshift+ix
      drs(ixGmin1:ixGmax1,ix)=dxfirst(level)*qstretch(level)**(index-1)
   end do

   pw(igrid)%dvolumep(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=drs(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2)*dx1

   pw(igrid)%dxp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,2)=drs(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2)
   pw(igrid)%dxp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)=dx1;

   if (need_only_volume) return

   
   
   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2)
   pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1) =drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2)
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=dx1
   pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=dx1
  
   

case ("spherical")
   x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)=pw(igrid)%xCCp(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2,1)
   
   x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,2)=pw(igrid)%xCCp(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2,2)
   ! spherical grid stretches the radial (first) coordinate
   if(stretched_grid) then
     ! grid level and global index for filling the grid itself
     level=node(plevel_,igrid)
     ig1=node(pig1_,igrid)
     ig2=node(pig2_,igrid)
     ! when filling the coarse grid (through need_only_volume) then adjust
     if(need_only_volume)then
        level=level-1
        ig1Co=(ig1-1)/2
        ixshift=ig1Co*block_nx1+(1-mod(ig1,2))*block_nx1/2-nghostcells
     else
        ixshift=(ig1-1)*block_nx1-nghostcells
     endif
     do ix = ixGmin1,ixGmax1
        index=ixshift+ix
        drs(ix,ixGmin2:ixGmax2)=dxfirst(level)*qstretch(level)**(index-1)
     end do
   else
     drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=dx1
   end if

   pw(igrid)%dvolumep(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=(x(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2,1)**2+drs(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2)**2/12.0d0)*drs(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2) *two*dabs(dsin(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
      2)))*dsin(half*dx2)

   pw(igrid)%dxp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)=drs(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2)
    pw(igrid)%dxp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,2)=x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,1)*dx2
   

   if (need_only_volume) return

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;

   pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=(x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1)+half*drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2))**2  *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      2))*dsin(half*dx2)

   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2)*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+half*dx2)

   

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1)**2  *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      2))*dsin(half*dx2)
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2)*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2))

   

case ("cylindrical")
   x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)=pw(igrid)%xCCp(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2,1)
   ! cylindrical grid stretches the radial (first) coordinate
   if(stretched_grid) then
     ! this is the grid level and global index, when filling the grid itself
     level=node(plevel_,igrid)
     ig1=node(pig1_,igrid)
     ig2=node(pig2_,igrid)
     ! when filling the coarse grid (through need_only_volume) then adjust
     if(need_only_volume)then
        level=level-1
        ig1Co=(ig1-1)/2
        ixshift=ig1Co*block_nx1+(1-mod(ig1,2))*block_nx1/2-nghostcells
     else
        ixshift=(ig1-1)*block_nx1-nghostcells
     endif
     do ix = ixGmin1,ixGmax1
       index=ixshift+ix
       drs(ix,ixGmin2:ixGmax2)=dxfirst(level)*qstretch(level)**(index-1)
     end do
   else
     drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=dx1
   end if

   pw(igrid)%dvolumep(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=dabs(x(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2,1))*drs(ixGmin1:ixGmax1,ixGmin2:ixGmax2)*dx2
   ! following is also equivalent to the above
   !pw(igrid)%dvolumep(ixG^S)=dabs(half*&
   !       ((x(ixG^S,1)+half*drs(ixG^S))**2-&
   !        (x(ixG^S,1)-half*drs(ixG^S))**2)){^DE&*dx^DE }

   pw(igrid)%dxp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)=drs(ixGmin1:ixGmax1,&
      ixGmin2:ixGmax2)

   if (z_ > 0) then
     if (2==z_) pw(igrid)%dxp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,2)=dx2
   end if

   if (phi_ > 0) then
     if (2==phi_) pw(igrid)%dxp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        2)=x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1)*dx2
   end if

   if (need_only_volume) return

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      1)=dabs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)+half*drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2))*dx2
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   if (z_==2) pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      2)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2)
   if (phi_==2) pw(igrid)%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      2)=drs(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=dabs(x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1))*dx2
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   if (z_==2) pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      2)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)*drs(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2)
   if (phi_==2) pw(igrid)%surface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
      2)=drs(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   

case default
   call mpistop("Sorry, typeaxial unknown")
end select

end subroutine fillgeo
!=============================================================================
subroutine gradient(q,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

use mod_global_parameters

integer :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,ixmax1,ixmax2, idir
double precision :: q(ixImin1:ixImax1,ixImin2:ixImax2), gradq(ixImin1:ixImax1,&
   ixImin2:ixImax2)

double precision :: qC(ixImin1:ixImax1,ixImin2:ixImax2),invdx
integer :: jxmin1,jxmin2,jxmax1,jxmax2, hxmin1,hxmin2,hxmax1,hxmax2, ixCmin1,&
   ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2 

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
if (slab) then

   jxmin1=ixmin1+kr(idir,1);jxmin2=ixmin2+kr(idir,2);jxmax1=ixmax1+kr(idir,1)
   jxmax2=ixmax2+kr(idir,2);
   hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmax1=ixmax1-kr(idir,1)
   hxmax2=ixmax2-kr(idir,2);
   gradq(ixmin1:ixmax1,ixmin2:ixmax2) = half*(q(jxmin1:jxmax1,&
      jxmin2:jxmax2)-q(hxmin1:hxmax1,hxmin2:hxmax2))*invdx

else
   hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmax1=ixmax1-kr(idir,1)
   hxmax2=ixmax2-kr(idir,2);
   ixCmin1=hxmin1;ixCmin2=hxmin2;ixCmax1=ixmax1;ixCmax2=ixmax2;
   jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
   jxCmax1=ixCmax1+kr(idir,1);jxCmax2=ixCmax2+kr(idir,2);
   qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,idir)*half*(q(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2)+q(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
   gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(qC(ixmin1:ixmax1,&
      ixmin2:ixmax2)-qC(hxmin1:hxmax1,hxmin2:hxmax2))/block%dvolume(&
      ixmin1:ixmax1,ixmin2:ixmax2)
end if

end subroutine gradient
!=============================================================================
subroutine gradient_s(q,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idims,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

use mod_global_parameters

integer :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
    idims
double precision :: q(ixImin1:ixImax1,ixImin2:ixImax2), gradq(ixImin1:ixImax1,&
   ixImin2:ixImax2)

double precision :: qC(ixImin1:ixImax1,ixImin2:ixImax2),qpoint(ixImin1:ixImax1,&
   ixImin2:ixImax2),qface(ixImin1:ixImax1,ixImin2:ixImax2),invdx
integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2, ixBmin1,ixBmin2,ixBmax1,ixBmax2,&
    ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, ix1,ix2

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idims)
! ixC is cell-corner index
ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;
if (slab) then
  ! b unit vector at cell corner
  qpoint=0.d0
  do ix2=0,1
  do ix1=0,1
    ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;
    ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;
    qpoint(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qpoint(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2)+q(ixAmin1:ixAmax1,ixAmin2:ixAmax2)
  end do
  end do
  qpoint(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=qpoint(ixCmin1:ixCmax1,&
     ixCmin2:ixCmax2)*0.5d0**ndim
  ! values at cell face
  qface=0.d0
  ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
  ixBmax1=ixOmax1-kr(idims,1);ixBmax2=ixOmax2-kr(idims,2);
  ixAmax1=ixOmax1;ixAmax2=ixOmax2; ixAmin1=ixBmin1;ixAmin2=ixBmin2;
  ixBmin1=ixAmin1;ixBmin2=ixAmin2;ixBmax1=ixAmax1;ixBmax2=ixAmax2;
  do ix2=0,1 
  do ix1=0,1 
     if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims ) then
       ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2; 
       ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2; 
       qface(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=qface(ixAmin1:ixAmax1,&
          ixAmin2:ixAmax2)+qpoint(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
     end if
  end do
  end do
  qface(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=qface(ixAmin1:ixAmax1,&
     ixAmin2:ixAmax2)*0.5d0**(ndim-1)
  ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
  ixBmax1=ixOmax1-kr(idims,1);ixBmax2=ixOmax2-kr(idims,2);
  gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=invdx*(qface(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)-qface(ixBmin1:ixBmax1,ixBmin2:ixBmax2))
end if

end subroutine gradient_s
!=============================================================================
subroutine gradientS(q,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir
! first use limiter to go from cell center to edge

use mod_global_parameters
use mod_limiter
use mod_ppm

integer :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,ixmax1,ixmax2, idir
double precision :: q(ixImin1:ixImax1,ixImin2:ixImax2), gradq(ixImin1:ixImax1,&
   ixImin2:ixImax2)

double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2):: qC,qL,qR,dqC,ldq,&
   rdq
double precision :: invdx
integer :: hxmin1,hxmin2,hxmax1,hxmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,jxCmin1,&
   jxCmin2,jxCmax1,jxCmax2,gxCmin1,gxCmin2,gxCmax1,gxCmax2,hxCmin1,hxCmin2,&
   hxCmax1,hxCmax2

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmax1=ixmax1-kr(idir,1)
hxmax2=ixmax2-kr(idir,2);
ixCmin1=hxmin1;ixCmin2=hxmin2;ixCmax1=ixmax1;ixCmax2=ixmax2;
jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
jxCmax1=ixCmax1+kr(idir,1);jxCmax2=ixCmax2+kr(idir,2);
gxCmin1=ixCmin1-kr(idir,1);gxCmin2=ixCmin2-kr(idir,2);gxCmax1=jxCmax1
gxCmax2=jxCmax2;
hxCmin1=gxCmin1+kr(idir,1);hxCmin2=gxCmin2+kr(idir,2)
hxCmax1=gxCmax1+kr(idir,1);hxCmax2=gxCmax2+kr(idir,2);

! set the gradient limiter here
qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = q(hxCmin1:hxCmax1,hxCmin2:hxCmax2)
qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = q(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
if (typegradlimiter/=limiter_ppm) then
   dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2)= qR(gxCmin1:gxCmax1,&
      gxCmin2:gxCmax2)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
   call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,gxCmax1,&
      gxCmax2,idir,typegradlimiter,ldw=ldq,rdw=rdq)
   qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qL(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2) + half*ldq(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qR(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2) - half*rdq(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
else
   call PPMlimitervar(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,&
      idir,q,q,qL,qR)
endif

if (slab) then
   gradq(ixmin1:ixmax1,ixmin2:ixmax2)=half*(qR(ixmin1:ixmax1,&
      ixmin2:ixmax2)-qL(hxmin1:hxmax1,hxmin2:hxmax2))*invdx
else
   select case(idir)
   case(1)
    gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(qR(ixmin1:ixmax1,&
       ixmin2:ixmax2)-qL(hxmin1:hxmax1,hxmin2:hxmax2))/block%dx(ixmin1:ixmax1,&
       ixmin2:ixmax2,idir) 
   case(2)
    gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(qR(ixmin1:ixmax1,&
       ixmin2:ixmax2)-qL(hxmin1:hxmax1,hxmin2:hxmax2))/block%dx(ixmin1:ixmax1,&
       ixmin2:ixmax2,idir) 
   end select
end if

end subroutine gradientS
!=============================================================================
subroutine divvector(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,divq)

! Calculate divergence of a vector qvec within ixL

use mod_global_parameters

integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
    divq(ixImin1:ixImax1,ixImin2:ixImax2)

double precision :: qC(ixImin1:ixImax1,ixImin2:ixImax2), invdx(1:ndim)
integer :: jxOmin1,jxOmin2,jxOmax1,jxOmax2, hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
    ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, idims,&
    ixmin1,ixmin2,ixmax1,ixmax2 
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2) call &
   mpistop("Error in divvector: Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
if (slab) then
  do idims=1,ndim

     jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
     jxOmax1=ixOmax1+kr(idims,1);jxOmax2=ixOmax2+kr(idims,2);
     hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
     hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+half*(qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
        idims)-qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,idims))*invdx(idims)

  end do
else
  do idims=1,ndim
     hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
     hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
     ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
     jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
     jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
     qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,idims)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
        idims)+qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,idims))
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+qC(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-qC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)
  end do
  divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
end if


end subroutine divvector 
!=============================================================================
subroutine curlvector(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,curlvec,idirmin,idirmin0,ndir0)

! Calculate curl of a vector qvec within ixL

use mod_global_parameters

integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   idirmin,ixmin1,ixmin2,ixmax1,ixmax2,idir,jdir,kdir,hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2,jxOmin1,jxOmin2,jxOmax1,jxOmax2,ndir0,idirmin0
double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir0),&
   curlvec(ixImin1:ixImax1,ixImin2:ixImax2,idirmin0:3), invdx(1:ndim)
double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
   ixImin2:ixImax2),mydx(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2) call &
   mpistop("Error in curl: Non-conforming input limits")

! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
! Curl can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

invdx=1.d0/dxlevel
idirmin=4
curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idirmin0:3)=zero

do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixmin1:ixmax1,ixmin2:ixmax2)=qvec(ixmin1:ixmax1,ixmin2:ixmax2,kdir)
      hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
      hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
      jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
      jxOmax1=ixOmax1+kr(jdir,1);jxOmax2=ixOmax2+kr(jdir,2);
      select case(typeaxial)
        case('slab')


         tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
            jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*invdx(jdir)

        case('slabstretch')
         if(jdir==2) then
           call gradient(tmp,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
              ixOmax1,ixOmax2,jdir,tmp2)
         else
           tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
              jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
              hxOmin2:hxOmax2))*invdx(jdir)
         end if
        case('spherical')
         select case(jdir)
            case(1)
             tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                ixImin2:ixImax2)*block%x(ixImin1:ixImax1,ixImin2:ixImax2,1)
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
                jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2))/(block%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                1)*block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
    case(2)
             mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=block%dx(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,2)
             if(idir==1) then
                tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                   ixImin2:ixImax2)*dsin(block%x(ixImin1:ixImax1,&
                   ixImin2:ixImax2,2))
                mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsin(block%x(&
                   ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*mydx(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)
             endif
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
                jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2))/mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
 
         end select
        case('cylindrical')
         if(z_==3) then
           select case(jdir)
              case(1)
               mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=block%dx(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,1)
               if(idir==3) then
                  tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)*block%x(ixImin1:ixImax1,ixImin2:ixImax2,&
                     1)
                  mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=block%x(&
                     ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*mydx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)
               endif
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
                  jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      case(2)
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
                  jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  2)
   
           end select
         end if
         if(phi_==3) then
           select case(jdir)
              case(1)
               mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=block%dx(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,1)
               if(idir==2) then
                  tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)*block%x(ixImin1:ixImax1,ixImin2:ixImax2,&
                     1)
                  mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=block%x(&
                     ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*mydx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)
               endif
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-&
                  half*(tmp(jxOmin1:jxOmax1,&
                  jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      case(2)
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-&
                  half*(tmp(jxOmin1:jxOmax1,&
                  jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  2)
   
           end select
         end if
      end select
      if(lvc(idir,jdir,kdir)==1)then
         curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=curlvec(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,idir)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
         curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=curlvec(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,idir)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      endif
      if(idir<idirmin)idirmin=idir
   endif
enddo; enddo; enddo;

end subroutine curlvector 
!=============================================================================
subroutine divvectorS(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,divq)

! Calculate divergence of a vector qvec within ixL
! using limited extrapolation to cell edges

use mod_global_parameters
use mod_limiter
use mod_ppm

integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
double precision :: qvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir),&
    divq(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: qL,qR,dqC,ldq,rdq
double precision :: invdx(1:ndim)

integer :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
   jxCmin1,jxCmin2,jxCmax1,jxCmax2,idims,ixmin1,ixmin2,ixmax1,ixmax2,gxCmin1,&
   gxCmin2,gxCmax1,gxCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2
!-----------------------------------------------------------------------------
ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2) call &
   mpistop("Error in divvectorS: Non-conforming input limits")

invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
do idims=1,ndim
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
   ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
   gxCmin1=ixCmin1-kr(idims,1);gxCmin2=ixCmin2-kr(idims,2);gxCmax1=jxCmax1
   gxCmax2=jxCmax2;
   hxCmin1=gxCmin1+kr(idims,1);hxCmin2=gxCmin2+kr(idims,2)
   hxCmax1=gxCmax1+kr(idims,1);hxCmax2=gxCmax2+kr(idims,2);
   qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = qvec(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
      idims)
   qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = qvec(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
      idims)
   if(typegradlimiter/=limiter_ppm) then
      dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2)= qR(gxCmin1:gxCmax1,&
         gxCmin2:gxCmax2)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
      call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
         gxCmax1,gxCmax2,idims,typegradlimiter,ldw=ldq,rdw=rdq)
      qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2) + half*ldq(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2) - half*rdq(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
   else
      dqC(ixImin1:ixImax1,ixImin2:ixImax2)=qvec(ixImin1:ixImax1,&
         ixImin2:ixImax2,idims)
      call PPMlimitervar(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
         ixMhi2,idims,dqC,dqC,qL,qR)
   endif

   if (slab) then
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+half*(qR(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*invdx(idims)
   else
     qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,idims)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
     qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%surfaceC(ixCmin1:ixCmax1,&
        ixCmin2:ixCmax2,idims)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+qR(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2)
   end if
end do
if(.not.slab) divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

end subroutine divvectorS
!=============================================================================
!> cross product of two vectors
subroutine cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,a,b,axb)
  use mod_global_parameters

  integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2
  double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,3),&
      b(ixImin1:ixImax1,ixImin2:ixImax2,3)
  double precision, intent(out) :: axb(ixImin1:ixImax1,ixImin2:ixImax2,3)
!-------------------------------------------------------------------------

  axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)-a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)
  axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     3)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)-a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)
  axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)=a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     1)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)-a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     2)*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)

end subroutine cross_product
!=============================================================================
subroutine extremaq(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,q,nshift,qMax,qMin)

use mod_global_parameters

integer,intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in) :: q(ixImin1:ixImax1,ixImin2:ixImax2)
integer,intent(in)           :: nshift

double precision, intent(out) :: qMax(ixImin1:ixImax1,ixImin2:ixImax2),&
   qMin(ixImin1:ixImax1,ixImin2:ixImax2)

integer           :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,&
   ixsRmax1,ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,&
   ishift,i,j 
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
 ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
 ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
 if (ishift==1) then
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(q(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(q(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
 else
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMax(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(qMin(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
 end if
 
 idims=1
 jdims=idims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmax1=ixOmax1+i*ishift*kr(idims,1);ixsmax2=ixOmax2+i*ishift*kr(idims,2);
   ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
   ixsRmax1=ixsmax1+ishift*kr(jdims,1);ixsRmax2=ixsmax2+ishift*kr(jdims,2);
   ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
   ixsLmax1=ixsmax1-ishift*kr(jdims,1);ixsLmax2=ixsmax2-ishift*kr(jdims,2);
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMax(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(qMin(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
 end do

 
enddo

end subroutine  extremaq
!=============================================================================
subroutine extremaw(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,nshift,wMax,wMin)

use mod_global_parameters

integer,intent(in)            :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
integer,intent(in)            :: nshift

double precision, intent(out) :: wMax(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux),wMin(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux)

integer          :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,ixsRmax1,&
   ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
 ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
 ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
 if (ishift==1) then
    wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= max(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
    wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= min(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
 else
    wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= max(wMax(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
    wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= min(wMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
 end if
 
 idims=1
 jdims=idims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmax1=ixOmax1+i*ishift*kr(idims,1);ixsmax2=ixOmax2+i*ishift*kr(idims,2);
   ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
   ixsRmax1=ixsmax1+ishift*kr(jdims,1);ixsRmax2=ixsmax2+ishift*kr(jdims,2);
   ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
   ixsLmax1=ixsmax1-ishift*kr(jdims,1);ixsLmax2=ixsmax2-ishift*kr(jdims,2);
   wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= max(wMax(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
      1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
   wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= min(wMin(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
      1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
 end do

 
enddo

end subroutine  extremaw
!=============================================================================
subroutine extremaa(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,a,nshift,aMin)

use mod_global_parameters

integer,intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2)
integer,intent(in)           :: nshift

double precision, intent(out) :: aMin(ixImin1:ixImax1,ixImin2:ixImax2)

integer          :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,ixsRmax1,&
   ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
  idims=1
  ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
  ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
  ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
  ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
  aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(a(ixsRmin1:ixsRmax1,&
     ixsRmin2:ixsRmax2),a(ixOmin1:ixOmax1,ixOmin2:ixOmax2),a(ixsLmin1:ixsLmax1,&
     ixsLmin2:ixsLmax2))
  
  idims=1
  jdims=idims+1
  do i=-1,1
    ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
    ixsmax1=ixOmax1+i*ishift*kr(idims,1);ixsmax2=ixOmax2+i*ishift*kr(idims,2);
    ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
    ixsRmax1=ixsmax1+ishift*kr(jdims,1);ixsRmax2=ixsmax2+ishift*kr(jdims,2);
    ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
    ixsLmax1=ixsmax1-ishift*kr(jdims,1);ixsLmax2=ixsmax2-ishift*kr(jdims,2);
    aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(aMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),a(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
       a(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
  end do
 
  
end do

end subroutine extremaa
!=============================================================================
subroutine locate_in_table(xpoint,table,nxtp,ixtp,res)
!  Fast search table to find xpoint's location in the table
!  record the distance between xpoint and the its cloest element table(ixtp)
! INPUT:
!  xpoint : the point you want to locate in the table
!  table : 1D table you want to find location in
!  nxtp : number of elements in the table
! OUTPUT:
!  ixtp : closest element's index in table
!  res : offset(distance) from the closest element in table
use mod_global_parameters

integer, intent(in) :: nxtp
double precision,intent(in)   :: xpoint,table(nxtp)
double precision, intent(out) :: res
integer, intent(out) :: ixtp

integer :: jl,jc,jh
!-----------------------------------------------------------------------------
if(xpoint < table(1)) then
  ixtp=1
  res=(xpoint-table(1))/(table(2)-table(1))
else if (xpoint > table(nxtp)) then
  ixtp=nxtp
  res=(xpoint-table(nxtp))/(table(nxtp)-table(nxtp-1))
else
  jl=0
  jh=nxtp+1
  do
    if (jh-jl <= 1) exit
    jc=(jh+jl)/2
    if (xpoint >= table(jc)) then
        jl=jc
    else
        jh=jc
    end if
  end do
  res=(xpoint-table(jh-1))/(table(jh)-table(jh-1))
  if(res<=0.5d0) then
    ixtp=jh-1
  else
    ixtp=jh
    res=res-1.d0
  endif
end if
end subroutine locate_in_table
!=============================================================================
subroutine divvector_s(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,divq)

! Calculate divergence of a vector qvec within ixL

use mod_global_parameters

integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
    divq(ixImin1:ixImax1,ixImin2:ixImax2)

double precision :: invdx(1:ndim)
double precision :: qpoint(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
   qface(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2, ixBmin1,ixBmin2,ixBmax1,ixBmax2,&
    ixCmin1,ixCmin2,ixCmax1,ixCmax2, idims, ixmin1,ixmin2,ixmax1,ixmax2,ix1,&
   ix2

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2<ixmax2) call &
   mpistop("Error in divvector: Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
! ixC is cell-corner index
ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;
if (slab) then
  ! vector at cell corner
  qpoint=0.d0
  do ix2=0,1
  do ix1=0,1
    ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;
    ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;
    qpoint(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:ndim)=qpoint(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,1:ndim)+qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndim)
  end do
  end do
  qpoint(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:ndim)=qpoint(ixCmin1:ixCmax1,&
     ixCmin2:ixCmax2,1:ndim)*0.5d0**ndim
  ! values at cell face
  qface=0.d0
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
         qface(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=qface(ixAmin1:ixAmax1,&
            ixAmin2:ixAmax2,idims)+qpoint(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
            idims)
       end if
    end do
    end do
    qface(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims)=qface(ixAmin1:ixAmax1,&
       ixAmin2:ixAmax2,idims)*0.5d0**(ndim-1)
  end do
  divq=0.d0
  do idims=1,ndim
    ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
    ixBmax1=ixOmax1-kr(idims,1);ixBmax2=ixOmax2-kr(idims,2);
    divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+invdx(idims)*(qface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idims)-qface(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idims))
  end do
else
  divq=0.d0
  do idims=1,ndim
    ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
    ixBmax1=ixOmax1-kr(idims,1);ixBmax2=ixOmax2-kr(idims,2);
    ixCmin1=ixBmin1;ixCmin2=ixBmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
    ixAmin1=ixCmin1+kr(idims,1);ixAmin2=ixCmin2+kr(idims,2)
    ixAmax1=ixCmax1+kr(idims,1);ixAmax2=ixCmax2+kr(idims,2);
    qface(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
       idims)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
       idims)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
       idims)+qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idims))
    divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+qface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idims)-qface(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idims)
  end do
  divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
end if

end subroutine divvector_s
