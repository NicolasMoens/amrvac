!> Module containing all the time stepping schemes
module mod_advance

  implicit none
  private

  logical :: firstsweep, lastsweep

  !> Whether to conserve fluxes at the current partial step
  logical :: fix_conserve_at_step = .true.

  public :: advance
  public :: process
  public :: process_advanced

contains

  !> Advance all the grids over one time step, including all sources
  subroutine advance(iit)
    use mod_global_parameters
    use mod_particles, only: handle_particles
    use mod_source, only: add_split_source

    integer, intent(in) :: iit

    integer :: iigrid, igrid, idimsplit

    ! old solution values at t_n-1 no longer needed: make copy of w(t_n)
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pw(igrid)%wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
          1:nwflux+nwaux)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
          1:nwflux+nwaux)
    end do
    !$OMP END PARALLEL DO

    

    ! split source addition
    !print*, "BEFORE add_split_source(prior=.true.)", pw(igrid)%w(5,5,iw_mom(1))
    call add_split_source(prior=.true.)
    !print*, "BEFORE add_split_source(prior=.true.)", pw(igrid)%w(5,5,iw_mom(1))

    firstsweep=.true.
    if (dimsplit) then
       if ((iit/2)*2==iit .or. typedimsplit=='xy') then
          ! do the sweeps in order of increasing idim,
          do idimsplit=1,ndim
             lastsweep= idimsplit==ndim
             call advect(idimsplit,idimsplit)
          end do
       else
          ! If the parity of "iit" is odd and typedimsplit=xyyx,
          ! do sweeps backwards
          do idimsplit=ndim,1,-1
             lastsweep= idimsplit==1
             call advect(idimsplit,idimsplit)
          end do
       end if
    else
       ! Add fluxes from all directions at once
       lastsweep= .true.
       !print*, BEFORE CALLING advect", pw(igrid)%w(5,5,iw_mom(1))
       call advect(1,ndim)
       !print*, AFTER CALLING advect", pw(igrid)%w(5,5,iw_mom(1)), pw(igrid)%w1(5,5,iw_mom(1))
    end if

    ! split source addition
    !print*, "BEFORE add_split_source(prior=.false.)"
    call add_split_source(prior=.false.)
    !print*, "AFTER add_split_source(prior=.false.)"

    if(use_particles) call handle_particles

  end subroutine advance

  !> Advance all grids over one time step, but without taking dimensional
  !> splitting or split source terms into account
  subroutine advect(idimmin,idimmax)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: idimmin,idimmax
    integer             :: iigrid, igrid

    call init_comm_fix_conserve(idimmin,idimmax,nwflux)
    fix_conserve_at_step = time_advance .and. levmax>levmin

    ! copy w instead of wold because of potential use of dimsplit or sourcesplit
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pw(igrid)%w1=pw(igrid)%w
    end do
    !$OMP END PARALLEL DO

    istep = 0

    select case (time_integrator)
    case ("onestep")
       !print*, "BEFORE CALLING advect1", pw(igrid)%w(5,5,iw_mom(1)), pw(igrid)%w1(5,5,iw_mom(1))
       call advect1(flux_scheme,one,idimmin,idimmax,global_time,1,global_time,&
          0)
       !print*, "AFTER CALLING advect1", pw(igrid)%w(5,5,iw_mom(1)), pw(igrid)%w1(5,5,iw_mom(1))

    case ("twostep")
      ! predictor step
       fix_conserve_at_step = .false.
       call advect1(typepred1,half, idimmin,idimmax,global_time,0,global_time,&
          1)

       ! corrector step
       fix_conserve_at_step = time_advance .and. levmax>levmin
       call advect1(flux_scheme,one,idimmin,idimmax,global_time+half*dt,1,&
          global_time,0)

    case ("twostep_trapezoidal")
       ! Explicit trapezoidal rule / Euler's method / Heun's method
       ! In the future, this could be implemented using just w and w1
       call advect1(flux_scheme, 1.0d0, idimmin,idimmax, global_time, 0,&
           global_time, 1)

       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w2 = pw(igrid)%w
          pw(igrid)%w = pw(igrid)%w1
       end do

       call advect1(flux_scheme, 1.0d0, idimmin,idimmax, global_time+dt, 1,&
           global_time+dt, 0)

       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w = 0.5d0 * (pw(igrid)%w + pw(igrid)%w2)
       end do
    case ("threestep")
       ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
       call advect1(flux_scheme,one, idimmin,idimmax,global_time,0,global_time,&
          1)

       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=0.75d0*pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)+0.25d0*pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)
          if (nw>nwflux) pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw)
       end do

       call advect1(flux_scheme,0.25d0, idimmin,idimmax,global_time+dt,1,&
          global_time+dt*0.25d0,2)

       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=1.0d0/3.0d0*pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)+2.0d0/3.0d0*pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,2.0d0/3.0d0, idimmin,idimmax,&
          global_time+dt/2.0d0,2,global_time+dt/3.0d0,0)

    case ("ssprk43")
       ! Strong stability preserving 4 stage RK 3rd order method by Ruuth and Spiteri
       !
       ! Ruuth & Spiteri
       ! J. S C, 17 (2002) p. 211 - 220
       !
       ! supposed to be stable up to CFL=2.
       ! don't use time-dependent sources since i did not bother to set those intermediate times.
       ! oliver.

       ! === First step ===
       call advect1(flux_scheme,0.5d0, idimmin,idimmax,global_time,0,&
          global_time,1)

       ! === Second step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
          if (nw>nwflux) pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw)
       end do
       call advect1(flux_scheme,0.5d0, idimmin,idimmax,global_time,1,&
          global_time+dt,2)

       ! === Third step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=2.0d0/3.0d0 * pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux) + 1.0d0/3.0d0 * pw(igrid)%w2(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,1:nwflux)
          if (nw>nwflux) pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw)
       end do
       call advect1(flux_scheme,1.0d0/6.0d0, idimmin,idimmax,global_time,2,&
          global_time+dt,3)

       ! === Fourth step ===
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.5d0, idimmin,idimmax, global_time,3,&
           global_time+dt,0)

    case ("ssprk54")
       ! Strong stability preserving 5 stage RK 4th order method by Ruuth and Spiteri
       !
       ! SIAM J. NUMER. ANAL.
       ! Vol. 40, No. 2, pp. 469–491
       ! c 2002 Society for Industrial and Applied Mathematics
       ! A NEW CLASS OF OPTIMAL
       ! HIGH-ORDER STRONG-STABILITY-PRESERVING
       ! TIME DISCRETIZATION METHODS
       !
       ! E.g. Table A.2
       !
       ! I have the actual coefficients however from the overview article by Gottlieb, JoSC 25 (2005)
       ! ON HIGH ORDER STRONG STABILITY PRESERVING RUNGE-KUTTA AND MULTI STEP TIME DISCRETIZATIONS
       !
       ! there are slight differences in the coefficients (~8th digit after the .)
       ! This is SSP till CFL number 1.508 which makes the effective CFL per step ceff=0.377,
       ! in contrast to the ceff=0.3333 of the classical RK3.
       !
       ! coded by oliver on 11/05/2013.  Enjoy!

       ! === First step ===
       call advect1(flux_scheme,0.391752226571890d0, idimmin,idimmax,&
           global_time,0,global_time,1)

       ! === Second step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=0.444370493651235d0 * pw(igrid)%w(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,1:nwflux) + 0.555629506348765d0 * &
             pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
          if (nw>nwflux) pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw)
       end do
       call advect1(flux_scheme,0.368410593050371d0, idimmin,idimmax,&
           global_time+0.391752226571890d0*dt,1,&
           global_time+0.2176690962611688d0*dt,2)

       ! === Third step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=0.620101851488403d0 * pw(igrid)%w(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,1:nwflux) + 0.379898148511597d0 * &
             pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
          if (nw>nwflux) pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw)
       end do
       call advect1(flux_scheme,0.251891774271694d0, idimmin,idimmax,&
           global_time+0.5860796893115398d0*dt,2,&
           global_time+0.222650588849706d0*dt,3)

       ! === Fourth step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w4(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=0.178079954393132d0 * pw(igrid)%w(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,1:nwflux) + 0.821920045606868d0 * &
             pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
          if (nw>nwflux) pw(igrid)%w4(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw) = pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             nwflux+1:nw)
       end do
       call advect1(flux_scheme,0.544974750228521d0, idimmin,idimmax,&
           global_time+0.4745423631214d0*dt,3,&
           global_time+0.390035880739132d0*dt,4)
       ! Now recover back the dt*L(u3), store in w1:
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux) = ( pw(igrid)%w4(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux) - (0.178079954393132d0 * pw(igrid)%w(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,1:nwflux) + 0.821920045606868d0 * &
             pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)) ) / 0.544974750228521d0
       end do
       !$OMP END PARALLEL DO

       ! === Fifth step ===
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)= 0.517231671970585d0 * pw(igrid)%w2(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,1:nwflux) + 0.096059710526147d0 * &
             pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux) + 0.063692468666290d0 * pw(igrid)%w1(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,1:nwflux) + 0.386708617503269d0 * &
             pw(igrid)%w4(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.226007483236906d0, idimmin,idimmax,&
           global_time+0.935010630967653d0*dt,4,&
           global_time+0.710300048096804d0*dt,0)


    case ("rk4")
       ! classical RK4 four step scheme
       call advect1(flux_scheme,0.5d0, idimmin,idimmax,global_time,0,&
          global_time,1)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       call advect1(flux_scheme,0.5d0, idimmin,idimmax,global_time+dt/2d0,1,&
          global_time,2)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       call advect1(flux_scheme,one,   idimmin,idimmax,global_time+dt/2d0,2,&
          global_time,3)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=(pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux) +two*pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux) +pw(igrid)%w3(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux) -4.0d0*pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux))/3.0d0
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,1.0d0/6.0d0, idimmin,idimmax,global_time+dt,3,&
          global_time,0)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)+pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       !$OMP END PARALLEL DO
    case ("fourstep")
       ! four step scheme, variant Hans De Sterck
       call advect1(flux_scheme,0.12d0, idimmin,idimmax,global_time,0,&
          global_time,1)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       call advect1(flux_scheme,0.25d0, idimmin,idimmax,global_time+dt*0.12d0,&
          1,global_time,2)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.5d0,  idimmin,idimmax,global_time+dt/4d0   ,&
          2,global_time,1)
       call advect1(flux_scheme,one,    idimmin,idimmax,global_time+dt/2d0   ,&
          1,global_time,0)
    case ("jameson")
       ! four step scheme, variant jameson
       call advect1(flux_scheme,0.25d0, idimmin,idimmax,global_time,0,&
          global_time,1)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w2(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       call advect1(flux_scheme,(1.0d0/3.0d0), idimmin,idimmax,&
          global_time+dt*0.25d0,1,global_time,2)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          pw(igrid)%w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             1:nwflux)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.5d0,  idimmin,idimmax,global_time+dt/3d0,2,&
          global_time,1)
       call advect1(flux_scheme,one,    idimmin,idimmax,global_time+dt/2d0,1,&
          global_time,0)
    case default
       write(unitterm,*) "time_integrator=",time_integrator
       write(unitterm,*) "Error in advect: Unknown time integration method"
       call mpistop("Correct time_integrator")
    end select

    firstsweep=.false.
  end subroutine advect

  !> Integrate all grids by one partial step
  subroutine advect1(method,dtfactor,idimmin,idimmax,qtC,a,qt,b)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics, only: phys_req_diagonal

    integer, intent(in) :: idimmin,idimmax
    integer, intent(in) :: a !< Compute fluxes based on this w
    integer, intent(in) :: b !< Update solution on this w
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: qt
    character(len=*), intent(in) :: method(nlevelshi)

    double precision :: qdt
    integer :: iigrid, igrid, level

    logical :: setigrid

    istep = istep+1

    ! opedit: Just advance the active grids:
    !$OMP PARALLEL DO PRIVATE(igrid,level,qdt)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       level=node(plevel_,igrid)
       qdt=dtfactor*dt_grid(igrid)
       select case(a)
       case(0)
         pw(igrid)%wa=>pw(igrid)%w
       case(1)
         pw(igrid)%wa=>pw(igrid)%w1
       case(2)
         pw(igrid)%wa=>pw(igrid)%w2
       case(3)
         pw(igrid)%wa=>pw(igrid)%w3
       case(4)
         pw(igrid)%wa=>pw(igrid)%w4
       end select

       select case(b)
       case(0)
         pw(igrid)%wb=>pw(igrid)%w
       case(1)
         pw(igrid)%wb=>pw(igrid)%w1
       case(2)
         pw(igrid)%wb=>pw(igrid)%w2
       case(3)
         pw(igrid)%wb=>pw(igrid)%w3
       case(4)
         pw(igrid)%wb=>pw(igrid)%w4
       end select

       !print*, "BEFORE CALLING process1_grid from advect1", pw(igrid)%wb(5,5,iw_mom(1)), pw(igrid)%wa(5,5,iw_mom(1))
       call process1_grid(method(level),igrid,qdt,ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
          idimmin,idimmax,qtC,pw(igrid)%wa,qt,pw(igrid)%wb,pw(igrid)%wold)
       !print*, "AFTER CALLING process1_grid from advect1", pw(igrid)%wb(5,5,iw_mom(1)), pw(igrid)%wa(5,5,iw_mom(1))

       !print*, "------------------a---------------", a, pw(igrid)%wa(5,5,iw_mom(1))
       !print*, "------------------b---------------", b, pw(igrid)%wb(5,5,iw_mom(1))

    end do
    !$OMP END PARALLEL DO

    ! opedit: Send flux for all grids, expects sends for all
    ! nsend_fc(^D), set in connectivity.t.

    if (fix_conserve_at_step) then
      call recvflux(idimmin,idimmax)
      call sendflux(idimmin,idimmax)
      call fix_conserve(idimmin,idimmax,1,nwflux)
    end if

    ! For all grids: fill ghost cells
    qdt = dtfactor*dt
    call getbc(qt+qdt,qdt,0,nwflux+nwaux, phys_req_diagonal)

  end subroutine advect1

  !> Prepare to advance a single grid over one partial time step
  subroutine process1_grid(method,igrid,qdt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
     idimmin,idimmax,qtC,wCT,qt,w,wold)
    use mod_global_parameters
    use mod_fix_conserve

    character(len=*), intent(in) :: method
    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmax1,ixGmax2, idimmin,&
       idimmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision :: wCT(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw),&
        w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw), wold(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,1:nw)

    double precision :: dx1,dx2
    double precision :: fC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nwflux,1:ndim)

    dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);
    dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
    saveigrid=igrid

    block=>pw(igrid)
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

    fC=0.d0

    !print*, "BEFORE CALLING advect1_grid from process1_grid", w(5,5,iw_mom(1))
    call advect1_grid(method,qdt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,idimmin,&
       idimmax,qtC,wCT,qt,w,wold,fC,dx1,dx2, pw(igrid)%x)
    !print*, "AFTER CALLING advect1_grid from process1_grid", w(5,5,iw_mom(1))

    ! opedit: Obviously, flux is stored only for active grids.
    ! but we know in fix_conserve wether there is a passive neighbor
    ! but we know in conserve_fix wether there is a passive neighbor
    ! via neighbor_active(i^D,igrid) thus we skip the correction for those.
    ! This violates strict conservation when the active/passive interface
    ! coincides with a coarse/fine interface.
    if (fix_conserve_at_step) then
      call storeflux(igrid,fC,idimmin,idimmax,nwflux)
    end if

  end subroutine process1_grid

  !> Advance a single grid over one partial time step
  subroutine advect1_grid(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,idimmin,&
     idimmax,qtC,wCT,qt,w,wold,fC,dx1,dx2,x)

    !  integrate one grid by one partial step
    use mod_finite_volume
    use mod_finite_difference
    use mod_tvd
    use mod_source, only: addsource2
    use mod_global_parameters

    character(len=*), intent(in) :: method
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, idimmin,idimmax
    double precision, intent(in) :: qdt, qtC, qt, dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    double precision :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw), wold(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)

    integer :: ixOmin1,ixOmin2,ixOmax1,ixOmax2

    ixOmin1=ixImin1+nghostcells;ixOmin2=ixImin2+nghostcells
    ixOmax1=ixImax1-nghostcells;ixOmax2=ixImax2-nghostcells;

    select case (method)
    case ('cd')
       call centdiff(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,fC,dx1,dx2,x)
    case ('cd4')
       call centdiff4(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,fC,dx1,dx2,x)
    case ('hancock')
       call hancock(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,dx1,dx2,x)
    case ('fd')
       call fd(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,wold,fC,dx1,dx2,x)
    case ('tvdmu','tvdlf','hll','hllc','hllcd','hlld')
       !print*, "BEFORE CALLING finite_volume from advect1_grid", w(5,5,iw_mom(1))
       call finite_volume(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,wold,fC,dx1,dx2,&
          x)
       !print*, "AFTER CALLING finite_volume from advect1_grid", w(5,5,iw_mom(1))
    case ('tvd')
       call centdiff(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,fC,dx1,dx2,x)
       call tvdlimit(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,idimmin,idimmax,wCT,qt+qdt,w,fC,dx1,dx2,x)
    case ('source')
       !print*, "CALLING addsource2 from mod_advance"
       call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImin2,&
          ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,w,x,&
          .false.)
    case ('nul')
       ! There is nothing to do
    case default
       write(unitterm,*)'Error in advect1_grid:',method,' is unknown!'
       call mpistop("")
    end select

  end subroutine advect1_grid

  !> process is a user entry in time loop, before output and advance
  !>         allows to modify solution, add extra variables, etc.
  !> Warning: CFL dt already determined (and is not recomputed)!
  subroutine process(iit,qt)
    use mod_usr_methods, only: usr_process_grid, usr_process_global
    use mod_global_parameters
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt

    integer:: iigrid, igrid,level

    if (associated(usr_process_global)) then
       call usr_process_global(iit,qt)
    end if

    !$OMP PARALLEL DO PRIVATE(igrid,level)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       level=node(plevel_,igrid)
       ! next few lines ensure correct usage of routines like divvector etc
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
       block=>pw(igrid)
       typelimiter=type_limiter(node(plevel_,igrid))
       typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

       if (associated(usr_process_grid)) then
          call usr_process_grid(igrid,level,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,&
             ixMlo2,ixMhi1,ixMhi2, qt,pw(igrid)%w,pw(igrid)%x)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine process

  !> process_advanced is user entry in time loop, just after advance
  !>           allows to modify solution, add extra variables, etc.
  !>           added for handling two-way coupled PIC-MHD
  !> Warning: w is now at t^(n+1), global time and iteration at t^n, it^n
  subroutine process_advanced(iit,qt)
    use mod_usr_methods, only: usr_process_adv_grid, usr_process_adv_global
    use mod_global_parameters
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt

    integer:: iigrid, igrid,level

    if (associated(usr_process_adv_global)) then
       call usr_process_adv_global(iit,qt)
    end if

    !$OMP PARALLEL DO PRIVATE(igrid,level)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       level=node(plevel_,igrid)
       ! next few lines ensure correct usage of routines like divvector etc
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
       block=>pw(igrid)
       typelimiter=type_limiter(node(plevel_,igrid))
       typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

       if (associated(usr_process_adv_grid)) then
          call usr_process_adv_grid(igrid,level,ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
             ixMlo1,ixMlo2,ixMhi1,ixMhi2, qt,pw(igrid)%w,pw(igrid)%x)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine process_advanced

end module mod_advance
