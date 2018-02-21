!> Module with finite difference methods for fluxes
module mod_finite_difference

  implicit none
  private

  public :: fd
  public :: centdiff
  public :: centdiff4

contains

  subroutine fd(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,x)
    use mod_physics
    use mod_source, only: addsource2
    use mod_global_parameters

    character(len=*), intent(in)                                     :: method
    double precision, intent(in)                                     :: qdt,&
        qtC, qt, dx1,dx2
    integer, intent(in)                                              :: &
       ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
        idimmin,idimmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        intent(in)            :: x

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        intent(inout)           :: wCT, wnew, wold
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
       1:ndim), intent(out)  :: fC

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux)                      :: fCT
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)                          :: fm, fp, fmR, fpL, wprim
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)                               :: v
    double precision                                                 :: &
       dxinv(1:ndim)
    logical                                                          :: &
       transport
    integer                                                          :: idims,&
        iw, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixmin1,ixmin2,ixmax1,ixmax2,&
        hxOmin1,hxOmin2,hxOmax1,hxOmax2, ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
    do idims= idimmin,idimmax

       block%iw0=idims

       ! Get fluxes for the whole grid (mesh+nghostcells)
        ixCmin1 = ixOmin1 - nghostcells * kr(idims,1)
         ixCmin2 = ixOmin2 - nghostcells * kr(idims,2)
        ixCmax1 = ixOmax1 + nghostcells * kr(idims,1)
         ixCmax2 = ixOmax2 + nghostcells * kr(idims,2)

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
       ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixmax1=ixOmax1;ixmax2=ixOmax2; ixmin1=hxOmin1;ixmin2=hxOmin2;

        ixCRmin1 = ixCmin1 + kr(idims,1)*phys_wider_stencil
         ixCRmin2 = ixCmin2 + kr(idims,2)*phys_wider_stencil
        ixCRmax1 = ixCmax1 - kr(idims,1)*phys_wider_stencil
         ixCRmax2 = ixCmax2 - kr(idims,2)*phys_wider_stencil

       wprim=wCT
       call phys_to_primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCRmin1,ixCRmin2,&
          ixCRmax1,ixCRmax2,wprim,x)
       call phys_get_flux(wCT,wprim,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCRmin1,&
          ixCRmin2,ixCRmax1,ixCRmax2,idims,fCT)

       do iw=1,nwflux
          ! Lax-Friedrich splitting:
          fp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,&
             iw) = half * (fCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,&
             iw) + tvdlfeps * cmax_global * wCT(ixCRmin1:ixCRmax1,&
             ixCRmin2:ixCRmax2,iw))
          fm(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,&
             iw) = half * (fCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,&
             iw) - tvdlfeps * cmax_global * wCT(ixCRmin1:ixCRmax1,&
             ixCRmin2:ixCRmax2,iw))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructL(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
          ixmax2,idims,fp,fpL)
       call reconstructR(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
          ixmax2,idims,fm,fmR)

       do iw=1,nwflux
          if (slab) then
             fC(ixmin1:ixmax1,ixmin2:ixmax2,iw,&
                idims) = dxinv(idims) * (fpL(ixmin1:ixmax1,ixmin2:ixmax2,&
                iw) + fmR(ixmin1:ixmax1,ixmin2:ixmax2,iw))
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
          else
             fC(ixmin1:ixmax1,ixmin2:ixmax2,iw,&
                idims)=-qdt*block%surfaceC(ixmin1:ixmax1,ixmin2:ixmax2,&
                idims) * (fpL(ixmin1:ixmax1,ixmin2:ixmax2,&
                iw) + fmR(ixmin1:ixmax1,ixmin2:ixmax2,iw))
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
                idims))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          end if
       end do ! iw loop

    end do !idims loop

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'fd')

  end subroutine fd

  subroutine reconstructL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,idims,w,wLC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision, intent(out)   :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw) 

    double precision                :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
        dwC(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: jxRmin1,jxRmin2,jxRmax1,jxRmax2,&
        ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, iw

    select case (typelimiter)
    case (limiter_mp5)
       call MP5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wLC)
    case default 

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
       kxCmax2=ixImax2-kr(idims,2);

       wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux) = w(kxCmin1:kxCmax1,&
          kxCmin2:kxCmax2,1:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
       jxRmax1=iLmax1+kr(idims,1);jxRmax2=iLmax2+kr(idims,2);

       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=iLmin1-kr(idims,1)
       ixCmin2=iLmin2-kr(idims,2);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

       do iw=1,nwflux
          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)

          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,typelimiter,ldw=ldw)

          wLC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wLC(iLmin1:iLmax1,iLmin2:iLmax2,&
             iw)+half*ldw(iLmin1:iLmax1,iLmin2:iLmax2)
       end do
    end select

  end subroutine reconstructL

  subroutine reconstructR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,idims,w,wRC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision, intent(out)   :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw) 

    double precision                :: rdw(ixImin1:ixImax1,ixImin2:ixImax2),&
        dwC(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: jxRmin1,jxRmin2,jxRmax1,jxRmax2,&
        ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,kxRmax2, iw

    select case (typelimiter)
    case (limiter_mp5)
       call MP5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wRC)
    case default 

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
       kxCmax2=ixImax2-kr(idims,2);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

       wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=w(kxRmin1:kxRmax1,&
          kxRmin2:kxRmax2,1:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
       jxRmax1=iLmax1+kr(idims,1);jxRmax2=iLmax2+kr(idims,2);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=iLmin1-kr(idims,1)
       ixCmin2=iLmin2-kr(idims,2);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

       do iw=1,nwflux
          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,typelimiter,rdw=rdw)

          wRC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wRC(iLmin1:iLmax1,iLmin2:iLmax2,&
             iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)
       end do
    end select

  end subroutine reconstructR

  subroutine centdiff(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,fC,dx1,dx2,x)

    ! Advance the iws flow variables from global_time to global_time+qdt within ixO^L by centered 
    ! differencing in space the dw/dt+dF_i(w)/dx_i=S type equation. 
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_physics
    use mod_limiter
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idimmin,idimmax
    double precision, intent(in) :: qdt, qtC, qt, dx1,dx2
    double precision :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw), wprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)

    double precision :: f(ixImin1:ixImax1,ixImin2:ixImax2, nwflux)
    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,&
       jxCmax2
    logical :: transport

    ! An extra layer is needed in each direction for which fluxes are added.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimmin,idimmax
       ixmin1=ixmin1-kr(idims,1);ixmin2=ixmin2-kr(idims,2)
       ixmax1=ixmax1+kr(idims,1);ixmax2=ixmax2+kr(idims,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       then
       call mpistop("Error in CentDiff: Non-conforming input limits")
    end if

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;

    ! Add fluxes to w
    do idims= idimmin,idimmax
       ixmin1=ixOmin1-kr(idims,1);ixmin2=ixOmin2-kr(idims,2)
       ixmax1=ixOmax1+kr(idims,1);ixmax2=ixOmax2+kr(idims,2); ixCmin1=ixmin1
       ixCmin2=ixmin2; ixCmax1=ixOmax1;ixCmax2=ixOmax2;
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2); 
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       call phys_get_flux(wCT,wprim,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,&
          ixmin2,ixmax1,ixmax2,idims,f)

       do iw=1,nwflux
          ! Center flux to interface
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=half*(f(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2, iw)+f(jxCmin1:jxCmax1,jxCmin2:jxCmax2, iw))
          if (slab) then
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
                idims)=dxinv(idims)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
                idims)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
          else
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
                idims)=-qdt*block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                idims)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
                idims))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          end if
       end do    !next iw
    end do       !next idims

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,w,x,&
       .false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,'centdiff')

  end subroutine centdiff

  subroutine centdiff4(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,fC,dx1,dx2,x)

    ! Advance the flow variables from global_time to global_time+qdt within ixO^L by
    ! fourth order centered differencing in space 
    ! for the dw/dt+dF_i(w)/dx_i=S type equation.
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_physics
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idimmin,idimmax
    double precision, intent(in) :: qdt, qtC, qt, dx1,dx2
    double precision :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)

    double precision :: v(ixImin1:ixImax1,ixImin2:ixImax2,ndim),&
        f(ixImin1:ixImax1,ixImin2:ixImax2, nwflux)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim,&
        wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: vLC,&
        phi, cmaxLC, cmaxRC

    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,&
       jxCmax2, hxCmin1,hxCmin2,hxCmax1,hxCmax2, kxCmin1,kxCmin2,kxCmax1,&
       kxCmax2, kkxCmin1,kkxCmin2,kkxCmax1,kkxCmax2, kkxRmin1,kkxRmin2,&
       kkxRmax1,kkxRmax2
    logical :: transport, new_cmax, patchw(ixImin1:ixImax1,ixImin2:ixImax2)

    ! two extra layers are needed in each direction for which fluxes are added.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimmin,idimmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
    end do

    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       then
       call mpistop("Error in CentDiff4: Non-conforming input limits")
    end if

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    ! Add fluxes to w
    do idims= idimmin,idimmax
       block%iw0=idims

       ixmin1=ixOmin1-2*kr(idims,1);ixmin2=ixOmin2-2*kr(idims,2)
       ixmax1=ixOmax1+2*kr(idims,1);ixmax2=ixOmax2+2*kr(idims,2); 
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       ixCmin1=hxOmin1;ixCmin2=hxOmin2; ixCmax1=ixOmax1;ixCmax2=ixOmax2;
       hxCmin1=ixCmin1-kr(idims,1);hxCmin2=ixCmin2-kr(idims,2)
       hxCmax1=ixCmax1-kr(idims,1);hxCmax2=ixCmax2-kr(idims,2); 
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2); 
       kxCmin1=ixCmin1+2*kr(idims,1);kxCmin2=ixCmin2+2*kr(idims,2)
       kxCmax1=ixCmax1+2*kr(idims,1);kxCmax2=ixCmax2+2*kr(idims,2); 

       kkxCmin1=ixImin1;kkxCmin2=ixImin2; kkxCmax1=ixImax1-kr(idims,1)
       kkxCmax2=ixImax2-kr(idims,2);
       kkxRmin1=kkxCmin1+kr(idims,1);kkxRmin2=kkxCmin2+kr(idims,2)
       kkxRmax1=kkxCmax1+kr(idims,1);kkxRmax2=kkxCmax2+kr(idims,2);
       wRC(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,&
          1:nwflux)=wprim(kkxRmin1:kkxRmax1,kkxRmin2:kkxRmax2,1:nwflux)
       wLC(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,&
          1:nwflux)=wprim(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,1:nwflux)

       call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wprim,wprim,&
          wLC,wRC,wLp,wRp,x,.false.)

       ! Calculate velocities from upwinded values
       call phys_get_cmax(wLC,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,idims,cmaxLC)
       call phys_get_cmax(wRC,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,idims,cmaxRC)
       ! now take the maximum of left and right states
       vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=max(cmaxRC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2),cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

       call phys_get_flux(wCT,wprim,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,&
          ixmin2,ixmax1,ixmax2,idims,f)

       do iw=1,nwflux
          ! Center flux to interface
          ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=(-f(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2, iw)+7.0d0*(f(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
              iw) + f(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))-f(hxCmin1:hxCmax1,&
             hxCmin2:hxCmax2, iw))/12.0d0
          ! add rempel dissipative flux, only second order version for now
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw,idims)-half*vLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2) *(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))

          if (slab) then
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
                idims)=dxinv(idims)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
                idims)
             ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)+f_(i-1))-f_(i-2)]/12
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
          else
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
                idims)=-qdt*block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                idims)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
                idims))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
          end if
       end do    !next iw
    end do       !next idims

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,w,x,&
       .false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,w,x,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,'centdiff4')

  end subroutine centdiff4

end module mod_finite_difference
