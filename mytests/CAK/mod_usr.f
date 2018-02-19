!> This is a template for a new user problem
!$AMRVAC_DIR/setup.pl -d=1 -p=hd

module mod_usr

  ! Include a physics module
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  ! All in cgs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  double precision, parameter :: solar_mass = 1.98d33
  double precision, parameter :: solar_radius = 6.96d10
  double precision, parameter :: solar_lum = 3.9d33
  
  double precision, parameter :: mass_proton = 1.67d-24
  double precision, parameter :: boltzman_cst = 1.38d-16
  
  double precision, parameter :: alpha = 0.67
  double precision, parameter :: qbar = 2000.0
  double precision, parameter :: mass = 50.0*solar_mass
  double precision, parameter :: radius = 1.39d12 !20.0*solar_radius
  double precision, parameter :: tempw = 40000.0
  double precision, parameter :: lum = 8.d5*solar_lum
  double precision, parameter :: kappa_e = 0.34
  double precision, parameter :: beta = 0.7
  
  double precision, parameter :: pi_dp = 3.141592
  double precision, parameter :: c_dp = 2.99d10
  double precision, parameter :: G_dp = 6.67d-8 
  double precision, parameter :: c_sound = 2.3d6

  double precision, parameter :: Gamma_e = &
     kappa_e*lum/(4.*pi_dp*G_dp*mass*c_dp)
  double precision, parameter :: typical_speed = &
     (2.*boltzman_cst*tempw/mass_proton)**0.5  
  double precision, parameter :: escape_speed = (2.*G_dp*mass/radius)**0.5  
  double precision, parameter :: M_dot = lum/c_dp**2. * alpha/(1-alpha) &
     *(qbar*Gamma_e/(1.-Gamma_e))**((1.-alpha)/alpha)

  !compute dinflo from massloss
  !double precision, parameter :: dinflo = 4. * M_dot/(4*pi_dp*radius**2*c_sound)
   double precision, parameter :: dinflo = 2.21d-12
contains

!==========================================================================================

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_usr_methods
    use mod_global_parameters
 
    call set_coordinate_system("spherical")
    
    print*,"soundspeed", c_sound
    print*,"M_dot", M_dot*365*24*60*60/solar_mass
    print*,"dinflo", dinflo
    
    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions
    usr_special_bc => special_bound
    usr_source => special_source
    usr_get_dt => special_dt
    !usr_internal_bc => internal_bound
        
    ! Active the physics module
    call hd_activate()
  end subroutine usr_init

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixGmin1,ixGmax1, ixmin1,ixmax1, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in)    :: x(ixGmin1:ixGmax1, ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1, nw)
    integer :: i
    double precision :: v_inf, velocity_field(ixGmin1:ixGmax1)

    ! Set initial values for w
 
    v_inf = (1-Gamma_e)**0.5*escape_speed*(alpha/(1-alpha))**0.5
   

    do i = ixGmin1,ixGmax1
      if (x(i,1) .gt. radius) then
        !Set initial velocity acc to beta law
        velocity_field(i) = v_inf*(1-(radius/x(i,1)))**beta 
 
        !Set initial density
        w(i, rho_) = M_dot/(4*pi_dp*velocity_field(i)*x(i,1)**2)

        !set momentum
        w(i, mom(1)) = w(i,rho_)*velocity_field(i)

      else 

        !Set initial velocity acc to beta law
        w(i,mom(1)) = 0

        !Set initial density
        w(i,rho_) = dinflo        
        
      end if
    end do
 
    open(3,file='test_init')
    do i = ixmin1,ixmax1
      write(3,222),i,x(i,1)/(radius),w(i, mom(1))/w(i, rho_)/c_sound,w(i,&
          rho_)
    end do
    close(3)
    222 format(i8,5e15.5E3)

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================
 
  subroutine special_bound(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,iB,w,x)  

    use mod_global_parameters
    
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iB
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: v(ixImin1:ixImax1)
    !integer :: ixI^L, ix2
    integer :: i
    
    !---------------------------------------------

    w(1,rho_) = dinflo
    w(2,rho_) = dinflo
    
    !>Equal
    !v(2) = v(3)
    !v(1) = v(3)

    !>interpolation
    v(2) = 2*(w(3,mom(1))/w(3,rho_)) - (w(4,mom(1))/w(4,rho_))  
    v(1) = 3*(w(3,mom(1))/w(3,rho_)) - 2*(w(4,mom(1))/w(4,rho_)) 
    
    !>Conservation of momentum
    !v(2) = (x(3,1)**2*w(3,mom(1)))/(x(2,1)**2*dinflo)
    !v(1) = (x(3,1)**2*w(3,mom(1)))/(x(1,1)**2*dinflo)
    
    !>Combination
    !v(2) = (x(3,1)**2*w(3,mom(1)))/(x(2,1)**2*dinflo)
    !v(1) = 2*(w(2,mom(1))/w(2,rho_)) - (w(3,mom(1))/w(3,rho_))

    !>V should be < c_sound
    if (v(2) > c_sound) then 
    v(2) = c_sound
    v(1) = c_sound
    end if
    if (v(2) < -c_sound) then
    v(2) = -c_sound
    v(1) = -c_sound
    end if

    !Velocity to momentum
    w(2,mom(1)) = v(2)*dinflo
    w(1,mom(1)) = v(1)*dinflo

  end subroutine special_bound

!==========================================================================================

  subroutine special_source(qdt,ixImin1,ixImax1,ixmin1,ixmax1,iwmin,iwmax,qtC,&
     wCT,qt,w,x)
    use mod_global_parameters
    
    integer, intent(in)             :: ixImin1,ixImax1, ixmin1,ixmax1, iwmin,&
       iwmax
    integer :: i
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: grad_fwd(ixImin1:ixImax1),grad_bwd(ixImin1:ixImax1),&
       grad_centCT(ixImin1:ixImax1)
    double precision :: g_CAK(ixImin1:ixImax1), g_grav_sc(ixImin1:ixImax1),&
        ppp(300)
    double precision :: vel(ixImin1:ixImax1)
    double precision :: fin_disc(ixImin1:ixImax1), mu_star(ixImin1:ixImax1),&
        sigma_star(ixImin1:ixImax1)
    double precision :: beta_fd(ixImin1:ixImax1), F_fd(ixImin1:ixImax1)
           
    !--------------------------------------------------------------------------------------
 
    vel(ixImin1:ixImax1) = wCT(ixImin1:ixImax1,mom(1))/wCT(ixImin1:ixImax1,&
       rho_)
    !print*, vel

    !--------------------------------------------------------------------------------------
    ! calculate grad_v
    
    do i = ixmin1+1,ixmax1-1 
      ! forward difference
      grad_fwd(i) = (vel(i+1)-vel(i))/(x(i+1,1)-x(i,1))   
      
      ! backward difference
      grad_bwd(i) = (vel(i)-vel(i-1))/(x(i,1)-x(i-1,1)) 
      
      ! central difference
      grad_centCT(i) = 0.5*abs(grad_fwd(i)+grad_bwd(i))

    end do
      grad_centCT(ixmin1) = abs( (wCT(ixmin1+1,mom(1))/wCT(ixmin1+1,&
          rho_)-wCT(ixmin1,mom(1))/wCT(ixmin1, rho_)))/(x(ixmin1+1,1)-x(ixmin1,&
         1))
      !grad_centCT(ixmin1) = abs( (vel(ixmin1+1)-vel(ixmin1))/(x(ixmin1+1,1)-x(ixmin1,1)) )
      grad_centCT(ixmax1) = abs( (wCT(ixmax1,mom(1))/wCT(ixmax1,&
          rho_)-wCT(ixmax1-1,mom(1))/wCT(ixmax1-1, rho_)))/(x(ixmax1,&
         1)-x(ixmax1-1,1))
      !grad_centCT(ixmin1) = abs( (vel(ixmax1)-vel(ixmax1-1))/(x(ixmax1,1)-x(ixmax1-1,1)) )   
   
    !--------------------------------------------------------------------------------------
    
    do i = ixImin1, ixImax1
      beta_fd(i) = (1 - (vel(i)/x(i,1))/grad_centCT(i))*(radius/x(i,1))**2
      if (beta_fd(i) .ge. 1.) then
        F_fd(i) = 1./(1+alpha)        
      else if (beta_fd(i).lt.-1.d10) then
        F_fd(i) = abs(beta)**alpha/(1+alpha)
      else if (abs(beta_fd(i)).gt.1.d-3) then
        F_fd(i) = (1.-(1.-beta_fd(i))**(1+alpha))/(beta_fd(i)*(1+alpha))
      else
        F_fd(i) = 1.-0.5*alpha*beta_fd(i)*(1.+0.333333*(1.-alpha)*beta_fd(i))
      end if
    end do
    F_fd(1) = F_fd(3) 
    F_fd(2) = F_fd(3)

    !--------------------------------------------------------------------------------------

    ! calculate g_CAK
    ! this is WHITOUT GEOMETRICAL CORRECTION!!!!!
    g_CAK(ixImin1:ixImax1) =0.2* 1./(1.-alpha) * &
       kappa_e*lum*qbar/(4.*pi_dp*x(ixImin1:ixImax1,&
       1)**2*c_dp) *(grad_centCT(ixImin1:ixImax1)/wCT(ixImin1:ixImax1,&
       rho_))**alpha *1./(c_dp*qbar*kappa_e)**alpha 
  
    g_grav_sc(ixImin1:ixImax1) = G_dp*mass*(1.-Gamma_e)/((x(ixImin1:ixImax1,&
       1))**2)
    g_CAK(1) = 0 
    g_CAK(2) = 0  

    !FINITE DISC CORRECTION
    do i = ixImin1,ixImax1
    g_CAK(i) = g_CAK(i) * F_fd(i)
    end do

    !--------------------------------

    do i = ixImin1,ixImax1
      if (g_CAK(i) /= g_CAK(i)) then
        print*, 'iteration', it
        print*, i,w(i, rho_),grad_centCT(i),wCT(i,rho_),g_CAK(i)
        print*, i,1/(1-alpha) * kappa_e*lum*qbar/(4*pi_dp*(x(i,1))**2*c_dp) *(grad_centCT(i)/(wCT(i,rho_)*c_dp*qbar*kappa_e))**alpha
        print*, i,(grad_centCT(i)/(wCT(i,rho_)))**alpha
        print*, fin_disc(i)
        print*, i,grad_centCT(i),wCT(i,rho_),(grad_centCT(i)/(wCT(i,rho_)))
        stop
      end if
    end do

    !--------------------------------
    !if (it == it_max-1) then
    if (global_time >= time_max-dt) then 
      print*, 'Calculation done, writing to file', it
      open(2,file='test_out')
      do i = ixImin1,ixImax1
        !print*,i,x(i,1)/(radius),w(i, mom(1))/w(i, rho_)/c_sound,w(i, rho_),g_CAK(i),g_grav_sc(i), g_grav_sc(i)-g_CAK(i)
        write(2,222),i,x(i,1)/(radius),w(i, mom(1))/w(i, rho_)/c_sound,w(i,&
            rho_),g_CAK(i),g_grav_sc(i),grad_centCT(i), F_fd(i), beta_fd(i)
      end do
      close(2) 
      222 format(i8,6e15.5E3,2e15.5)  
    print*, it
    end if  
    
    !w = w + qdt*gsource 
    w(ixmin1:ixmax1,mom(1)) = w(ixmin1:ixmax1,&
       mom(1)) + qdt * (g_CAK(ixmin1:ixmax1) - g_grav_sc(ixmin1:ixmax1)) * &
       wCT(ixmin1:ixmax1,rho_) 
   
  end subroutine special_source

!==========================================================================================

  subroutine special_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)

    use mod_global_parameters
    
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: dtnew    

    integer :: i
    double precision :: grad_fwd(ixImin1:ixImax1),grad_bwd(ixImin1:ixImax1),&
       grad_centCT(ixImin1:ixImax1)
    double precision :: g_CAK(ixImin1:ixImax1), g_grav_sc(ixImin1:ixImax1),&
        g_sc(ixImin1:ixImax1), g_grav(ixImin1:ixImax1), g_tot(ixImin1:ixImax1)

    !--------------------------------------------------------------------------------------
    !RECALCULATE G_CAK FOR DT DETIRMINATION, VERY INEFFICIENT
    !CALCULATED AT T, NOT CT. THIS IS A GOOD PROXY
      !--------------------------------------------------------------------------------------
      ! calculate grad_v

      do i = iximin1+1,iximax1-1

        ! forward difference
        grad_fwd(i) = (w(i+1,mom(1))/w(i+1, rho_)-w(i,mom(1))/w(i,&
            rho_))/(x(i+1,1)-x(i,1))
        ! backward difference
        grad_bwd(i) = (w(i,mom(1))/w(i, rho_)-w(i-1,mom(1))/w(i-1, rho_))/(x(i,&
           1)-x(i-1,1))
        ! central difference
        grad_centCT(i) = 0.5*(grad_fwd(i)+grad_bwd(i))
      end do
        grad_centCT(iximin1) = (w(iximin1+1,mom(1))/w(iximin1+1,&
            rho_)-w(iximin1,mom(1))/w(iximin1, rho_))/(x(iximin1+1,&
           1)-x(iximin1,1))
        grad_centCT(iximax1) = (w(iximax1,mom(1))/w(iximax1, rho_)-w(iximax1-1,&
           mom(1))/w(iximax1-1, rho_))/(x(iximax1,1)-x(iximax1-1,1))


      !---------------------------------------------------------------------------------------
      !calculate g_tot

      g_CAK(iximin1:iximax1) = 1/(1-alpha) * &
         kappa_e*lum*qbar/(4*pi_dp*(x(iximin1:iximax1,&
         1))**2*c_dp)* (grad_centCT(iximin1:iximax1)/(w(iximin1:iximax1,&
         rho_)*c_dp*qbar*kappa_e))**alpha
      g_CAK(1) = 0
      g_CAK(2) = 0
      g_CAK(3) = 0

      !g_grav_sc(ixi^S) = G_dp*mass*(1-Gamma_e)/((x(ixi^S,1))**2)
      g_grav(iximin1:iximax1) = G_dp*mass*(1)/((x(iximin1:iximax1,1))**2)
      g_sc(iximin1:iximax1) = G_dp*mass*(Gamma_e)/((x(iximin1:iximax1,1))**2)

      g_tot(iximin1:iximax1) = -g_CAK(iximin1:iximax1) - g_sc(iximin1:iximax1) &
         + g_grav(iximin1:iximax1)

      !------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------
 
    !set dt for g_cak
    
    !dt = sqrt(dx/g)
    dtnew = 0.3* minval( abs(dx(:,1)/g_tot(iximin1:iximax1))**0.5 )   

  if (global_time >= time_max-dt) then 
  !if (it == it_max) then 
    print*, 'g_grav'
    print*, minval(g_grav),',',maxval(g_grav), 'm/s^2'
    print*, 'g_sc'
    print*, minval(g_sc),',',maxval(g_sc), 'm/s^2'
    print*, 'g_CAK'
    print*, minval(g_CAK),',',maxval(g_CAK) , 'm/s^2'
    print*, 'old dt'
    print*, dt, 's'
    print*,'new dt'
    print*, dtnew, 's'
    print*, 'global time'
    print*, global_time, 's'
    print*,'dx/c'
    print*, dx(:,1)/(c_sound + maxval(w(ixImin1:ixImax1,mom(1))/w(ixImin1:ixImax1,rho_))), 's'
    print*, '10 dynamical wind timescales'
    print*, 10.d0*radius/maxval(w(ixImin1:ixImax1,mom(1))/w(ixImin1:ixImax1,rho_)), 's'
    print*,'scale height 1-Gamma'
    print*, c_sound**2/(g_grav(3)-g_sc(3)), 'cm'
    print*,'scale height in stellar radii'
    print*, c_sound**2/(g_grav(3)-g_sc(3))/radius
    print*, 'stellar radius in scaleheights'
    print*, radius/c_sound**2*(g_grav(3)-g_sc(3))
    print*, 'distance in grid'
    print*, (x(iximax1,1)-x(iximin1,1))
    print*, 'cellwidth'
    print*, dx(1,1)
    print*, 'nr of cells'
    print*, (x(iximax1,1)-x(iximin1,1))/dx(:,1)
    print*, 'advised nr of cells'
    print*, (x(iximax1,1)-x(iximin1,1))*(g_grav(3)-g_sc(3))/c_sound**2
 !stop
 end if

  end subroutine special_dt

!===========================================================================================

end module mod_usr
