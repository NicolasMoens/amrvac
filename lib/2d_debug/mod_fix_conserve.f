!> Module for flux conservation near refinement boundaries
module mod_fix_conserve
  implicit none
  private

  integer, save                        :: nrecv, nsend
  double precision, allocatable, save  :: recvbuffer(:)
  integer, save                        :: isize(2)
  integer, dimension(:), allocatable   :: fc_recvreq, fc_sendreq
  integer, dimension(:,:), allocatable :: fc_recvstat, fc_sendstat

  public :: init_comm_fix_conserve
  public :: allocateBflux
  public :: deallocateBflux
  public :: sendflux
  public :: recvflux
  public :: storeflux
  public :: fix_conserve

 contains

   subroutine init_comm_fix_conserve(idimmin,idimmax,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax,nwfluxin

     integer :: iigrid, igrid, idims, iside, i1,i2, nxCo1,nxCo2
     integer :: ic1,ic2, inc1,inc2, ipe_neighbor
     integer :: ibuf, recvsize

     nsend    = 0
     nrecv    = 0
     recvsize = 0

     do idims= idimmin,idimmax
       select case (idims)
         case (1)
         nrecv=nrecv+nrecv_fc(1)
         nsend=nsend+nsend_fc(1)
         nxCo1=1;nxCo2=ixGhi2/2-nghostcells;
         isize(1)=nxCo1*nxCo2*(nwfluxin)
         recvsize=recvsize+nrecv_fc(1)*isize(1)
         case (2)
         nrecv=nrecv+nrecv_fc(2)
         nsend=nsend+nsend_fc(2)
         nxCo1=ixGhi1/2-nghostcells;nxCo2=1;
         isize(2)=nxCo1*nxCo2*(nwfluxin)
         recvsize=recvsize+nrecv_fc(2)*isize(2)
       end select
     end do

     ! Reallocate buffers when size differs
     if (allocated(recvbuffer)) then
       if (recvsize /= size(recvbuffer)) then
         deallocate(recvbuffer)
         allocate(recvbuffer(recvsize))
       end if
     else
       allocate(recvbuffer(recvsize))
     end if

     if (allocated(fc_recvreq)) then
       if (nrecv /= size(fc_recvreq)) then
         deallocate(fc_recvreq, fc_recvstat)
         allocate(fc_recvstat(MPI_STATUS_SIZE,nrecv), fc_recvreq(nrecv))
       end if
     else
       allocate(fc_recvstat(MPI_STATUS_SIZE,nrecv), fc_recvreq(nrecv))
     end if

     if (allocated(fc_sendreq)) then
       if (nsend /= size(fc_sendreq)) then
         deallocate(fc_sendreq, fc_sendstat)
         allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
       end if
     else
       allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
     end if

   end subroutine init_comm_fix_conserve

   subroutine recvflux(idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, i1,i2, nxCo1,nxCo2
     integer :: ic1,ic2, inc1,inc2, ipe_neighbor
     integer :: ibuf, recvsize

     if (nrecv>0) then
       fc_recvreq=MPI_REQUEST_NULL
       ibuf=1
       irecv=0

       do iigrid=1,igridstail; igrid=igrids(iigrid);
         do idims= idimmin,idimmax
           do iside=1,2
             i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);

             if (phi_ > 0) then
               if (neighbor_pole(i1,i2,igrid)/=0) cycle
             end if

             if (neighbor_type(i1,i2,igrid)/=4) cycle
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
             inc2=2*i2+ic2
             do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
             if (ipe_neighbor/=mype) then
               irecv=irecv+1
               itag=4**2*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)
               call MPI_IRECV(recvbuffer(ibuf),isize(idims),&
                   MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                  fc_recvreq(irecv),ierrmpi)
               ibuf=ibuf+isize(idims)
             end if
             end do
             end do
           end do
         end do
       end do
     end if

   end subroutine recvflux

   subroutine allocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside, i1,i2, nx1,nx2, nxCo1,nxCo2
     !-----------------------------------------------------------------------------
     nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
     nxCo1=nx1/2;nxCo2=nx2/2;

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do iside=1,2
         i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

         if (phi_ > 0) then
           if (neighbor_pole(i1,i2,igrid)/=0) cycle
         end if

         select case (neighbor_type(i1,i2,igrid))
         case (4)
           allocate(pflux(iside,1,igrid)%flux(1,1:nx2,1:nwflux))
         case (2)
           allocate(pflux(iside,1,igrid)%flux(1,1:nxCo2,1:nwflux))
         end select
       end do
       do iside=1,2
         i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

         if (phi_ > 0) then
           if (neighbor_pole(i1,i2,igrid)/=0) cycle
         end if

         select case (neighbor_type(i1,i2,igrid))
         case (4)
           allocate(pflux(iside,2,igrid)%flux(1:nx1,1,1:nwflux))
         case (2)
           allocate(pflux(iside,2,igrid)%flux(1:nxCo1,1,1:nwflux))
         end select
       end do
     end do

   end subroutine allocateBflux

   subroutine deallocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside
     !-----------------------------------------------------------------------------
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do iside=1,2
         if (associated(pflux(iside,1,igrid)%flux)) then
           deallocate(pflux(iside,1,igrid)%flux)
           nullify(pflux(iside,1,igrid)%flux)
         end if
       end do
       do iside=1,2
         if (associated(pflux(iside,2,igrid)%flux)) then
           deallocate(pflux(iside,2,igrid)%flux)
           nullify(pflux(iside,2,igrid)%flux)
         end if
       end do
     end do

   end subroutine deallocateBflux

   subroutine fix_conserve(idimmin,idimmax,nw0,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax, nw0, nwfluxin

     integer :: iigrid, igrid, idims, iside, iotherside, i1,i2, ic1,ic2, inc1,&
        inc2, ixmin1,ixmin2,ixmax1,ixmax2
     integer :: nxCo1,nxCo2, iw, ix, ipe_neighbor, ineighbor, ibuf, ibufnext,&
         nw1
     double precision :: CoFiratio
     !-----------------------------------------------------------------------------

     nw1=nw0-1+nwfluxin
     if (slab) then
       ! The flux is divided by volume of fine cell. We need, however,
       ! to divide by volume of coarse cell => muliply by volume ratio
       CoFiratio=one/dble(2**ndim)
     end if

     if (nrecv>0) then
       call MPI_WAITALL(nrecv,fc_recvreq,fc_recvstat,ierrmpi)
       ibuf=1
     end if

     nxCo1=(ixMhi1-ixMlo1+1)/2;nxCo2=(ixMhi2-ixMlo2+1)/2;

     !if (.false.) then 

     ! for all grids: perform flux update at Coarse-Fine interfaces
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims= idimmin,idimmax
         select case (idims)
           case (1)
           do iside=1,2
             i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

             if (phi_ > 0) then
               if (neighbor_pole(i1,i2,igrid)/=0) cycle
             end if

             if (neighbor_type(i1,i2,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is 
             ! not computed, keep the buffer counter up to date:
             if (.not.neighbor_active(i1,i2,igrid).or..not.neighbor_active(0,0,&
                igrid) ) then
               do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               if (ipe_neighbor/=mype) then
                 ibufnext=ibuf+isize(1)
                 ibuf=ibufnext
               end if
               end do
           end do
               cycle
             end if
             !

             select case (iside)
             case (1)
               ix=ixMlo1
             case (2)
               ix=ixMhi1
             end select

             ! remove coarse flux
             if (slab) then
               pw(igrid)%wb(ix,ixMlo2:ixMhi2,nw0:nw1) = pw(igrid)%wb(ix,&
                  ixMlo2:ixMhi2,nw0:nw1) -pflux(iside,1,igrid)%flux(1,:,&
                  1:nwfluxin)
             else
               do iw=nw0,nw1
                 pw(igrid)%wb(ix,ixMlo2:ixMhi2,iw)=pw(igrid)%wb(ix,&
                    ixMlo2:ixMhi2,iw)-pflux(iside,1,igrid)%flux(1,:,&
                    iw-nw0+1) /pw(igrid)%dvolume(ix,ixMlo2:ixMhi2)
               end do
             end if


             ! add fine flux
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
             inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ineighbor=neighbor_child(1,inc1,inc2,igrid)
             ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
             ixmin1=ix;ixmin2=ixMlo2+(ic2-1)*nxCo2;
             ixmax1=ix;ixmax2=ixmin2-1+nxCo2;
             if (ipe_neighbor==mype) then
               iotherside=3-iside
               if (slab) then
                 pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1) = pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1) + pflux(iotherside,1,ineighbor)%flux(:,:,&
                    1:nwfluxin)* CoFiratio
               else
                 do iw=nw0,nw1
                   pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw)=pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw) +pflux(iotherside,1,ineighbor)%flux(:,:,&
                      iw-nw0+1) /pw(igrid)%dvolume(ixmin1:ixmax1,&
                      ixmin2:ixmax2)
                 end do
               end if
             else
               if (slab) then
                 ibufnext=ibuf+isize(1)
                 pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1) = pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1)+CoFiratio *reshape(source=recvbuffer(&
                    ibuf:ibufnext-1), shape=shape(pw(igrid)%wb(ixmin1:ixmax1,&
                    ixmin2:ixmax2,nw0:nw1)))
                 ibuf=ibufnext
               else
                 do iw=nw0,nw1
                   ibufnext=ibuf+isize(1)/(nwfluxin)
                   pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw)=pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw) +reshape(source=recvbuffer(ibuf:ibufnext-1),&
                       shape=shape(pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw))) /pw(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2)
                   ibuf=ibufnext
                 end do
               end if
             end if
             end do
           end do
           end do
           case (2)
           do iside=1,2
             i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

             if (phi_ > 0) then
               if (neighbor_pole(i1,i2,igrid)/=0) cycle
             end if

             if (neighbor_type(i1,i2,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is 
             ! not computed, keep the buffer counter up to date:
             if (.not.neighbor_active(i1,i2,igrid).or..not.neighbor_active(0,0,&
                igrid) ) then
               do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               if (ipe_neighbor/=mype) then
                 ibufnext=ibuf+isize(2)
                 ibuf=ibufnext
               end if
               end do
           end do
               cycle
             end if
             !

             select case (iside)
             case (1)
               ix=ixMlo2
             case (2)
               ix=ixMhi2
             end select

             ! remove coarse flux
             if (slab) then
               pw(igrid)%wb(ixMlo1:ixMhi1,ix,&
                  nw0:nw1) = pw(igrid)%wb(ixMlo1:ixMhi1,ix,&
                  nw0:nw1) -pflux(iside,2,igrid)%flux(:,1,1:nwfluxin)
             else
               do iw=nw0,nw1
                 pw(igrid)%wb(ixMlo1:ixMhi1,ix,iw)=pw(igrid)%wb(ixMlo1:ixMhi1,&
                    ix,iw)-pflux(iside,2,igrid)%flux(:,1,&
                    iw-nw0+1) /pw(igrid)%dvolume(ixMlo1:ixMhi1,ix)
               end do
             end if


             ! add fine flux
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
             inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ineighbor=neighbor_child(1,inc1,inc2,igrid)
             ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
             ixmin1=ixMlo1+(ic1-1)*nxCo1;ixmin2=ix;
             ixmax1=ixmin1-1+nxCo1;ixmax2=ix;
             if (ipe_neighbor==mype) then
               iotherside=3-iside
               if (slab) then
                 pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1) = pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1) + pflux(iotherside,2,ineighbor)%flux(:,:,&
                    1:nwfluxin)* CoFiratio
               else
                 do iw=nw0,nw1
                   pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw)=pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw) +pflux(iotherside,2,ineighbor)%flux(:,:,&
                      iw-nw0+1) /pw(igrid)%dvolume(ixmin1:ixmax1,&
                      ixmin2:ixmax2)
                 end do
               end if
             else
               if (slab) then
                 ibufnext=ibuf+isize(2)
                 pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1) = pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                    nw0:nw1)+CoFiratio *reshape(source=recvbuffer(&
                    ibuf:ibufnext-1), shape=shape(pw(igrid)%wb(ixmin1:ixmax1,&
                    ixmin2:ixmax2,nw0:nw1)))
                 ibuf=ibufnext
               else
                 do iw=nw0,nw1
                   ibufnext=ibuf+isize(2)/(nwfluxin)
                   pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw)=pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw) +reshape(source=recvbuffer(ibuf:ibufnext-1),&
                       shape=shape(pw(igrid)%wb(ixmin1:ixmax1,ixmin2:ixmax2,&
                      iw))) /pw(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2)
                   ibuf=ibufnext
                 end do
               end if
             end if
             end do
           end do
           end do
         end select
       end do
     end do

     if (nsend>0) then
       call MPI_WAITALL(nsend,fc_sendreq,fc_sendstat,ierrmpi)
     end if

   end subroutine fix_conserve

   subroutine storeflux(igrid,fC,idimmin,idimmax,nwfluxin)
     use mod_global_parameters

     integer, intent(in)          :: igrid, idimmin,idimmax, nwfluxin
     double precision, intent(in) :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwfluxin,&
        1:ndim)

     integer :: idims, iside, i1,i2, ic1,ic2, inc1,inc2, ix1,ix2, ixCo1,ixCo2,&
         nxCo1,nxCo2, iw
     !integer :: ineighbor, ipe_neighbor
     !-----------------------------------------------------------------------------
     do idims = idimmin,idimmax
       select case (idims)
         case (1)
         do iside=1,2
           i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

           if (phi_ > 0) then
             if (neighbor_pole(i1,i2,igrid)/=0) cycle
           end if

           select case (neighbor_type(i1,i2,igrid))
           case (4)
             select case (iside)
             case (1)
               pflux(iside,1,igrid)%flux(1,:,1:nwfluxin) = -fC(nghostcells,&
                  ixMlo2:ixMhi2,1:nwfluxin,1)
             case (2)
               pflux(iside,1,igrid)%flux(1,:,1:nwfluxin) = fC(ixMhi1,&
                  ixMlo2:ixMhi2,1:nwfluxin,1)
             end select
           case (2)
             nxCo1=1;nxCo2=ixGhi2/2-nghostcells;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                 do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                 ix1=nghostcells;ix2=ixMlo2+2*(ixCo2-1);
                 pflux(iside,1,igrid)%flux(ixCo1,ixCo2,iw) = sum(fC(ix1,&
                    ix2:ix2+1,iw,1))
                 end do
         end do
               end do
             case (2)
               do iw=1,nwfluxin
                 do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                 ix1=ixMhi1;ix2=ixMlo2+2*(ixCo2-1);
                 pflux(iside,1,igrid)%flux(ixCo1,ixCo2,iw) =-sum(fC(ix1,&
                    ix2:ix2+1,iw,1))
                 end do
         end do
               end do
             end select
           end select
         end do
         case (2)
         do iside=1,2
           i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

           if (phi_ > 0) then
             if (neighbor_pole(i1,i2,igrid)/=0) cycle
           end if

           select case (neighbor_type(i1,i2,igrid))
           case (4)
             select case (iside)
             case (1)
               pflux(iside,2,igrid)%flux(:,1,1:nwfluxin) = -fC(ixMlo1:ixMhi1,&
                  nghostcells,1:nwfluxin,2)
             case (2)
               pflux(iside,2,igrid)%flux(:,1,1:nwfluxin) = fC(ixMlo1:ixMhi1,&
                  ixMhi2,1:nwfluxin,2)
             end select
           case (2)
             nxCo1=ixGhi1/2-nghostcells;nxCo2=1;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                 do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                 ix1=ixMlo1+2*(ixCo1-1);ix2=nghostcells;
                 pflux(iside,2,igrid)%flux(ixCo1,ixCo2,iw) = sum(fC(ix1:ix1+1,&
                    ix2,iw,2))
                 end do
         end do
               end do
             case (2)
               do iw=1,nwfluxin
                 do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                 ix1=ixMlo1+2*(ixCo1-1);ix2=ixMhi2;
                 pflux(iside,2,igrid)%flux(ixCo1,ixCo2,iw) =-sum(fC(ix1:ix1+1,&
                    ix2,iw,2))
                 end do
         end do
               end do
             end select
           end select
         end do
       end select
     end do

   end subroutine storeflux

   subroutine sendflux(idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax

     integer :: idims, iside, i1,i2, ic1,ic2, inc1,inc2, ix1,ix2, ixCo1,ixCo2,&
         nxCo1,nxCo2, iw
     integer :: ineighbor, ipe_neighbor, igrid, iigrid

     fc_sendreq = MPI_REQUEST_NULL
     isend      = 0

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims = idimmin,idimmax
         select case (idims)
           case (1)
           do iside=1,2
             i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

             if (phi_ > 0) then
               if (neighbor_pole(i1,i2,igrid)/=0) cycle
             end if

             if (neighbor_type(i1,i2,igrid)==2) then

               ineighbor=neighbor(1,i1,i2,igrid)
               ipe_neighbor=neighbor(2,i1,i2,igrid)
               if (ipe_neighbor/=mype) then
                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2);
                 inc1=-2*i1+ic1;inc2=ic2;
                 itag=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                 isend=isend+1

                 call MPI_ISEND(pflux(iside,1,igrid)%flux,isize(1),&
                     MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                    fc_sendreq(isend),ierrmpi)
               end if
             end if
           end do
           case (2)
           do iside=1,2
             i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

             if (phi_ > 0) then
               if (neighbor_pole(i1,i2,igrid)/=0) cycle
             end if

             if (neighbor_type(i1,i2,igrid)==2) then

               ineighbor=neighbor(1,i1,i2,igrid)
               ipe_neighbor=neighbor(2,i1,i2,igrid)
               if (ipe_neighbor/=mype) then
                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2);
                 inc1=ic1;inc2=-2*i2+ic2;
                 itag=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                 isend=isend+1

                 call MPI_ISEND(pflux(iside,2,igrid)%flux,isize(2),&
                     MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                    fc_sendreq(isend),ierrmpi)
               end if
             end if
           end do
         end select
       end do
     end do
   end subroutine sendflux

end module mod_fix_conserve
