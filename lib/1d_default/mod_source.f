!> Module for handling split source terms (split from the fluxes)
module mod_source

  implicit none
  private

  public :: add_split_source
  public :: addsource2

contains

  subroutine add_split_source(prior)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_thermal_conduction, only: phys_thermal_conduction
    use mod_physics, only: phys_req_diagonal

    logical, intent(in) :: prior

    double precision :: qdt, qt
    integer :: iigrid, igrid, i1
    logical :: src_active

    ! add thermal conduction
    if(associated(phys_thermal_conduction)) call phys_thermal_conduction()

    src_active = .false.

    if ((.not.prior).and.(typesourcesplit=='sf' .or. typesourcesplit=='ssf')) &
       return

    if (prior) then
       qt=global_time
    else
       qt=global_time+dt
    end if
    !$OMP PARALLEL DO PRIVATE(igrid,qdt,i1)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       qdt=dt_grid(igrid)
       block=>pw(igrid)
       do i1=-1,1
       if (i1==0) cycle
       if (neighbor_type(i1,igrid)==2 .or. neighbor_type(i1,igrid)==4) then
          leveljump(i1)=.true.
       else
          leveljump(i1)=.false.
       end if
       end do
       !print*, BEFORE CALLING addsource1_grid", pw(igrid)%w(5,5,iw_mom(1))
       call addsource1_grid(igrid,qdt,qt,pw(igrid)%w,src_active)
       !print*, AFTER CALLING addsource1_grid", pw(igrid)%w(5,5,iw_mom(1))
    end do
    !$OMP END PARALLEL DO

    if (src_active) then
       call getbc(qt,0.d0,0,nwflux+nwaux, phys_req_diagonal)
    end if

  end subroutine add_split_source

  subroutine addsource1_grid(igrid,qdt,qt,w,src_active)

    use mod_global_parameters

    integer, intent(in) :: igrid
    double precision, intent(in) :: qdt, qt
    double precision, intent(inout) :: w(ixGlo1:ixGhi1,nw)
    logical, intent(inout) :: src_active

    double precision :: w1(ixGlo1:ixGhi1,nw)

    saveigrid=igrid
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

    dxlevel(1)=rnode(rpdx1_,igrid);

    w1(ixGlo1:ixGhi1,1:nw)=w(ixGlo1:ixGhi1,1:nw)

    select case (typesourcesplit)
    case ('sf')
       call addsource2(qdt  ,ixGlo1,ixGhi1,ixMlo1,ixMhi1,1,nw,qt,w1,qt,w,&
          pw(igrid)%x,.true.,src_active)
    case ('sfs')
       !print*, ===============================================================", it
       !print*, "CALLING addsource2 from addsource1_grid", w(5,5,iw_mom(1)), "active", src_active
       call addsource2(qdt/2,ixGlo1,ixGhi1,ixMlo1,ixMhi1,1,nw,qt,w1,qt,w,&
          pw(igrid)%x,.true.,src_active)
       !print*, "AFTER CALLING addsource2 from addsource1_grid", w(5,5,iw_mom(1)), "active", src_active
    case ('ssf')
       call addsource2(qdt/2,ixGlo1,ixGhi1,ixGlo1,ixGhi1,1,nw,qt,w,qt,w1,&
          pw(igrid)%x,.true.,src_active)
       call addsource2(qdt  ,ixGlo1,ixGhi1,ixMlo1,ixMhi1,1,nw,qt,w1,qt,w,&
          pw(igrid)%x,.true.,src_active)
    case ('ssfss')
       call addsource2(qdt/4,ixGlo1,ixGhi1,ixGlo1,ixGhi1,1,nw,qt,w,qt,w1,&
          pw(igrid)%x,.true.,src_active)
       call addsource2(qdt/2,ixGlo1,ixGhi1,ixMlo1,ixMhi1,1,nw,qt,w1,qt,w,&
          pw(igrid)%x,.true.,src_active)
    case default
       write(unitterm,*)'No such typesourcesplit=',typesourcesplit
       call mpistop("Error: Unknown typesourcesplit!")
    end select

  end subroutine addsource1_grid

  subroutine addsource2(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
     wCT,qt,w,x,qsourcesplit,src_active)

    ! Add source within ixO for iws: w=w+qdt*S[wCT]

    use mod_global_parameters
    use mod_physics, only: phys_add_source
    use mod_usr_methods, only: usr_source
    ! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

    integer, intent(in)              :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
        iwmin,iwmax
    double precision, intent(in)     :: qdt, qtC, qt
    double precision, intent(in)     :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout)  :: w(ixImin1:ixImax1,1:nw)
    logical, intent(in)              :: qsourcesplit
    logical, intent(inout), optional :: src_active
    logical                          :: tmp_active

    tmp_active = .false.
    ! user defined sources, typically explicitly added
    if ((qsourcesplit .eqv. source_split_usr) .and. associated(usr_source)) &
       then
       tmp_active = .true.
       call usr_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,&
          qt,w,x)
    end if

    ! physics defined sources, typically explicitly added,
    ! along with geometrical source additions
    !print*, "CALLING PHYS_ADD_SOURCE", w(5,5,iw_mom(1)), qsourcesplit, tmp_active
    call phys_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
       qsourcesplit,tmp_active)
    !print*, "AFTER CALLING PHYS_ADD_SOURCE", w(5,5,iw_mom(1)), qsourcesplit, tmp_active
    if (present(src_active)) src_active = src_active .or. tmp_active
  end subroutine addsource2

end module mod_source
