subroutine set_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  call set_B0_cell(pw(igrid)%B0(:,:,:,0),pw(igrid)%x,ixGlo1,ixGlo2,ixGhi1,&
     ixGhi2,ixGlo1,ixGlo2,ixGhi1,ixGhi2)
  call set_J0_cell(igrid,pw(igrid)%J0,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1-1,&
     ixMlo2-1,ixMhi1+1,ixMhi2+1)
  call set_B0_face(igrid,pw(igrid)%x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,&
     ixMhi1,ixMhi2)

end subroutine set_B0_grid

subroutine set_B0_cell(wB0,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2)
  use mod_usr_methods, only: usr_set_B0
  use mod_global_parameters
  
  integer, intent(in):: ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
     ixmax2
  double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndir)
  double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

  wB0(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero

  ! approximate cell-averaged B0 as cell-centered B0
  select case (typeaxial)
  case ("spherical")
     
     if (dabs(Bdip)>smalldouble) then
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=2.0d0*Bdip*dcos(x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**3
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=Bdip*dsin(x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**3
     end if
  
     if (abs(Bquad)>smalldouble) then
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           1) +Bquad*0.5d0*(1.0d0+3.0d0*dcos(2.0d0*x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2)))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**4
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           2)+Bquad*dsin(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
           2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**4
     end if
     if (abs(Boct)>smalldouble) then
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           1) +Boct*(10.0d0*dcos(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
           2))-2.0d0) *dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,2))/x(ixmin1:ixmax1,&
           ixmin2:ixmax2,1)**5
        wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,&
           2) +Boct*1.5d0*(3.0d0+5.0d0*dcos(2.0d0*x(ixmin1:ixmax1,&
           ixmin2:ixmax2,2))) *dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,&
           2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**5
     end if
    
  end select
  
  if (associated(usr_set_B0)) call usr_set_B0(ixImin1,ixImin2,ixImax1,ixImax2,&
     ixmin1,ixmin2,ixmax1,ixmax2,x,wB0)

end subroutine set_B0_cell

subroutine set_J0_cell(igrid,wJ0,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2)
  use mod_usr_methods, only: usr_set_J0
  use mod_global_parameters

  integer, intent(in):: igrid,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
     ixmax1,ixmax2
  double precision, intent(inout) :: wJ0(ixImin1:ixImax1,ixImin2:ixImax2,&
     7-2*ndir:3)
  integer :: idirmin0, idirmin

  if(associated(usr_set_J0)) then
    call usr_set_J0(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
       ixmax2,pw(igrid)%x,wJ0)
  else
    idirmin0 = 7-2*ndir
    call curlvector(pw(igrid)%B0(:,:,:,0),ixImin1,ixImin2,ixImax1,ixImax2,&
       ixmin1,ixmin2,ixmax1,ixmax2,wJ0,idirmin,idirmin0,ndir)
  end if

end subroutine set_J0_cell

subroutine set_B0_face(igrid,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2)
  use mod_global_parameters

  integer, intent(in) :: igrid, ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,&
     ixmax1,ixmax2
  double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

  double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),dxu(ndim),&
     xshift(ndim)
  integer :: idims, ixCmin1,ixCmin2,ixCmax1,ixCmax2, idims2

  dxu(1)=rnode(rpdx1_,igrid);dxu(2)=rnode(rpdx2_,igrid);

  do idims=1,ndim
     ixCmin1=ixmin1-kr(1,idims);ixCmin2=ixmin2-kr(2,idims); ixCmax1=ixmax1
     ixCmax2=ixmax2;
     xshift(1:ndim)=half*(one-kr(1:ndim,idims))
     do idims2=1,ndim
       xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idims2)=x(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,idims2)+(half-xshift(idims2))*dxu(idims2)
     end do
     if(stretched_grid) then
       if(slab_stretched.and.idims==2) then
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=x(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,2)+0.5d0*pw(igrid)%dx(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,2)
       else if(.not.slab_stretched.and.idims==1) then
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=x(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,1)+0.5d0*pw(igrid)%dx(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,1)
       end if
     end if
     call set_B0_cell(pw(igrid)%B0(:,:,:,idims),xC,ixImin1,ixImin2,ixImax1,&
        ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2)
  end do

end subroutine set_B0_face

subroutine alloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  if(.not. allocated(pw(igrid)%B0)) then
    allocate(pw(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir,0:ndim))
    allocate(pw(igrid)%J0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2*ndir:3))
  end if

end subroutine alloc_B0_grid

subroutine dealloc_B0_grid(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  deallocate(pw(igrid)%B0)
  deallocate(pw(igrid)%J0)

end subroutine dealloc_B0_grid
