{#IFDEF PARTICLES
!=============================================================================
subroutine integrate_particles()

!-----------------------------------------------------------------------------
{#IFDEF PARTICLES_LORENTZ
call integrate_particles_lorentz
}
{#IFDEF PARTICLES_ADVECT
call integrate_particles_advect
}
{#IFDEF PARTICLES_GCA
call integrate_particles_gca
}
{#IFDEF PARTICLES_USER
call integrate_particles_user
}

end subroutine integrate_particles
!=============================================================================
subroutine set_particles_dt()

!-----------------------------------------------------------------------------
{#IFDEF PARTICLES_LORENTZ
call set_particles_dt_lorentz
}
{#IFDEF PARTICLES_ADVECT
call set_particles_dt_advect
}
{#IFDEF PARTICLES_GCA
call set_particles_dt_gca
}
{#IFDEF PARTICLES_USER
call set_particles_dt_user
}

end subroutine set_particles_dt
!=============================================================================





{#IFDEF PARTICLES_GCA
!=============================================================================
subroutine integrate_particles_gca()

use mod_odeint
use mod_particles
use constants
use mod_gridvars
include 'amrvacdef.f'

logical                             :: int_choice
integer                             :: ipart, iipart, seed, ic^D,igrid_particle, ipe_particle, ipe_working
logical                                   :: BC_applied
double precision                    :: lfac, absS
double precision                    :: B_new(1:ndir), absB_new, u(1:ndir)
double precision                    :: r^C, s^C, prob, theta
double precision                    :: dt_p, tloc, dydt(1:^NC+2),ytmp(1:^NC+2), euler_cfl, int_factor
double precision, dimension(1:ndir) :: x, ue, e, b, bhat, x_new
double precision, dimension(1:ndir) :: drift1, drift2
double precision, dimension(1:ndir) :: drift3, drift4, drift5, drift6, drift7
double precision, dimension(1:ndir) :: bdotgradb, uedotgradb, gradkappaB
double precision, dimension(1:ndir) :: bdotgradue, uedotgradue
double precision, dimension(1:ndir) :: gradBdrift, reldrift, bdotgradbdrift
double precision, dimension(1:ndir) :: uedotgradbdrift, bdotgraduedrift
double precision, dimension(1:ndir) :: uedotgraduedrift
double precision                    :: kappa, Mr, upar, m, absb, gamma, q, mompar, vpar, ueabs
double precision                    :: gradBdrift_abs, reldrift_abs, epar
double precision                    :: bdotgradbdrift_abs, uedotgradbdrift_abs
double precision                    :: bdotgraduedrift_abs, uedotgraduedrift_abs
double precision                    :: momentumpar1, momentumpar2, momentumpar3, momentumpar4
! for odeint:
integer                             :: nok, nbad, ic1^D, ic2^D
double precision                    :: h1, hmin, h_old

!!!!! Precision of time-integration: !!!!!!!!!!!!!
double precision,parameter          :: eps=1.0d-6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter                  :: nvar = ^NC+2
double precision,dimension(1:nvar)  :: y
external derivs_gca, derivs_gca_rk
!-----------------------------------------------------------------------------
do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);
int_choice=.false.   
   dt_p = particle(ipart)%self%dt
   igrid_working = particle(ipart)%igrid
   ipart_working = particle(ipart)%self%index
   tloc = particle(ipart)%self%t
   x(1:ndir) = particle(ipart)%self%x(1:ndir)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> Adaptive stepwidth RK4:

   !> initial solution vector:
   y(1:ndir) = x(1:ndir) !< position of guiding center
   y(ndir+1) = particle(ipart)%self%u(1) !< parallel momentum component (gamma v||)
   y(ndir+2) = particle(ipart)%self%u(2) !< conserved magnetic moment Mr
!   y(ndir+3) = particle(ipart)%self%u(3) !< Lorentz factor of particle

!>we temporarily save the solution vector, to replace the one from the euler
!>timestep after euler integration
   ytmp=y

call derivs_gca(particle(ipart)%self%t,y,dydt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> make an Euler step with the proposed timestep:

   !> factor to ensure we capture all particles near the internal ghost cells. Can be adjusted during a run, after an interpolation error
   euler_cfl=2.5d0

   !> new solution vector:
   y(1:ndir+2) = y(1:ndir+2) + euler_cfl * dt_p * dydt(1:ndir+2)
   particle(ipart)%self%x(1:ndir) = y(1:ndir) ! position of guiding center
   particle(ipart)%self%u(1)      = y(ndir+1) ! parallel momentum component(gamma v||)
   particle(ipart)%self%u(2)      = y(ndir+2) ! conserved magnetic moment

!> check if the particle is in the internal ghost cells
int_factor =1.0d0

if (.not. particle_in_igrid(ipart_working,igrid_working)) then

!>if particle is not in the grid an euler timestep is taken instead of a RK4
!>timestep. Then based on that we do an interpolation and check how much further
!>the timestep for the RK4 has to be restricted.

     !>factor to make integration more accurate for particles near the internal
     !>ghost cells. This factor can be changed during integration after an
     !>interpolation error. But one should be careful with timesteps for i/o

!> flat interpolation:
{ic^D = int((y(^D)-rnode(rpxmin^D_,igrid_working))/rnode(rpdx^D_,igrid_working)) + 1 + dixB\}

!> linear interpolation:
{
if (px(igrid_working)%x({ic^DD},^D) .lt. y(^D)) then
   ic1^D = ic^D
else
   ic1^D = ic^D -1
end if
ic2^D = ic1^D + 1
\}

int_factor =0.5d0

{^D&
if (ic1^D .le. ixGlo^D-2 .or. ic2^D .ge. ixGhi^D+2) then
  int_factor = 0.05d0
end if
\}

{^D&
if (ic1^D .eq. ixGlo^D-1 .or. ic2^D .eq. ixGhi^D+1) then
  int_factor = 0.1d0
end if
\}

     dt_p=int_factor*dt_p
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> replace the solution vector with the original as it was before the Euler timestep
   y(1:ndir+2) = ytmp(1:ndir+2)

   particle(ipart)%self%x(1:ndir) = ytmp(1:ndir) !< position of guiding center
   particle(ipart)%self%u(1)      = ytmp(ndir+1) !< parallel momentum component (gamma v||)
   particle(ipart)%self%u(2)      = ytmp(ndir+2) !< conserved magnetic moment


!< specify a minimum step hmin. If the timestep reaches this minimum, multiply by
!< a factor 100 to make sure the RK integration doesn't crash
   h1 = dt_p/2.0d0; hmin=1.0d-9; h_old=dt_p/2.0d0

   if(h1 .lt. hmin)then
           h1=hmin
           dt_p=2.0d0*h1
   endif

!< RK4 integration with adaptive stepwidth
   call odeint(y,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_gca_rk,rkqs)

!< original RK integration without interpolation in ghost cells
! call odeint(y,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_gca,rkqs)

!< final solution vector after rk integration
   particle(ipart)%self%x(1:ndir) = y(1:ndir)
   particle(ipart)%self%u(1)      = y(ndir+1)
   particle(ipart)%self%u(2)      = y(ndir+2)
   !particle(ipart)%self%u(3)      = y(ndir+3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> now calculate other quantities, mean Lorentz factor, drifts, perpendicular velocity:
   call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,b,bp1_,bp^NC_)
   call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,e,ep1_,ep^NC_)

   absb         = sqrt( {^C& b(^C)**2 |+} )
   bhat(1:ndir) = b(1:ndir) / absb

   epar         = {^C& e(^C)*bhat(^C) |+}

   call cross(e,bhat,ue)
   ue(1:ndir)   = ue(1:ndir)*CONST_c / absb
   ueabs = sqrt({^C& ue(^C)**2|+})

   kappa = sqrt(1.0d0 - ({^C& ue(^C)**2|+})/CONST_c**2)

   Mr = y(ndir+2); upar = y(ndir+1); m=particle(ipart)%self%m; q=particle(ipart)%self%q

   gamma = sqrt(1.0d0+upar**2/CONST_c**2+2.0d0*Mr*absb/m/CONST_c**2)/kappa

   particle(ipart)%self%u(3)      = gamma

   vpar = particle(ipart)%self%u(1)/particle(ipart)%self%u(3)
   mompar = particle(ipart)%self%u(1)
   
   call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,bdotgradb,b_dot_grad_b1_,b_dot_grad_b^NC_)
   call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,uedotgradb,ue_dot_grad_b1_,ue_dot_grad_b^NC_)
   call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,gradkappaB,grad_kappa_B1_,grad_kappa_B^NC_)
   call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,bdotgradue,b_dot_grad_ue1_,b_dot_grad_ue^NC_)
   call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,uedotgradue,ue_dot_grad_ue1_,ue_dot_grad_ue^NC_)

   drift1(1:ndir) = bhat(1:ndir)/(absb*kappa**2)

   drift2(1:ndir) = Mr*CONST_C/(gamma*q)*gradkappaB(1:ndir)
   call cross(drift1,drift2,gradBdrift)
   gradBdrift_abs = sqrt({^C& gradBdrift(^C)**2|+})

   drift3(1:ndir) = upar*epar/(gamma*CONST_C)*ue(1:ndir)
   call cross(drift1,drift3,reldrift)
   reldrift_abs = sqrt({^C& reldrift(^C)**2|+})

   drift4(1:ndir) = m*CONST_C/q* ( upar**2/gamma*bdotgradb(1:ndir))
   call cross(drift1,drift4,bdotgradbdrift)
   bdotgradbdrift_abs = sqrt({^C& bdotgradbdrift(^C)**2|+})

   drift5(1:ndir) = m*CONST_C/q* ( upar*uedotgradb(1:ndir))
   call cross(drift1,drift5,uedotgradbdrift)
   uedotgradbdrift_abs = sqrt({^C& uedotgradbdrift(^C)**2|+})

   drift6(1:ndir) = m*CONST_C/q* ( upar*bdotgradue(1:ndir))
   call cross(drift1,drift6,bdotgraduedrift)
   bdotgraduedrift_abs = sqrt({^C& bdotgraduedrift(^C)**2|+})

   drift7(1:ndir) = m*CONST_C/q* (gamma*uedotgradue(1:ndir))
   call cross(drift1,drift7,uedotgraduedrift)
   uedotgraduedrift_abs = sqrt({^C& uedotgraduedrift(^C)**2|+})

   momentumpar1 = m*gamma*vpar**2*({^C& ue(^C)*bdotgradb(^C)|+})
   momentumpar2 = m*gamma*vpar*({^C& ue(^C)*uedotgradb(^C)|+})
   momentumpar3 = q*epar
   momentumpar4 = -(Mr/gamma)*({^C& bhat(^C)*gradkappaB(^C)|+})



! **************************************************
!> Payload update
! **************************************************

      !> current gyroradius
      particle(ipart)%self%payload(1) = sqrt(2.0d0*m*Mr*absb)/abs(q*absb)*CONST_c

      !> pitch angle
      particle(ipart)%self%payload(2) = datan(sqrt((2.0d0*Mr*absb)/(m*gamma**2))/vpar)

      !> particle v_perp
      particle(ipart)%self%payload(3) = sqrt((2.0d0*Mr*absb)/(m*gamma**2))

      !> particle parallel momentum term 1
      particle(ipart)%self%payload(4) = momentumpar1

      !> particle parallel momentum term 2
      particle(ipart)%self%payload(5) = momentumpar2

      !> particle parallel momentum term 3
      particle(ipart)%self%payload(6) = momentumpar3

      !> particle parallel momentum term 4
      particle(ipart)%self%payload(7) = momentumpar4
     
      !> particle ExB drift
      particle(ipart)%self%payload(8) = ueabs 

      !> relativistic drift
      particle(ipart)%self%payload(9) = reldrift_abs

      !> gradB drift
      particle(ipart)%self%payload(10) = gradBdrift_abs

      !> bdotgradb drift
      particle(ipart)%self%payload(11) = bdotgradbdrift_abs

      !> uedotgradb drift
      particle(ipart)%self%payload(12) = uedotgradbdrift_abs
      
      !> bdotgradue drift
      particle(ipart)%self%payload(13) = bdotgraduedrift_abs

      !> uedotgradue drift
      particle(ipart)%self%payload(14) = uedotgraduedrift_abs


   ! **************************************************
   !> Time update
   ! **************************************************
   particle(ipart)%self%t = particle(ipart)%self%t + dt_p

end do

end subroutine integrate_particles_gca
!=============================================================================
subroutine derivs_gca_rk(t_s,y,dydt)

use mod_odeint
use mod_particles
use mod_gridvars
use mod_particles, only: particle
use constants
include 'amrvacdef.f'

integer                         :: ic^D
double precision                :: t_s, y(*)
double precision                :: dydt(*)
! .. local ..
double precision,dimension(^NC) :: ue, b, e, x, bhat, bdotgradb, uedotgradb, gradkappaB
double precision,dimension(^NC) :: bdotgradue, uedotgradue, u, utmp1, utmp2, utmp3
double precision                :: upar, Mr, gamma, absb, q, m, epar, kappa
!-----------------------------------------------------------------------------
!> Here the terms in the guiding centre equations of motion are interpolated for
!> the RK integration. The interpolation is also done in the ghost cells such
!> that the RK integration does not give an error

q = particle(ipart_working)%self%q
m = particle(ipart_working)%self%m

x(1:ndir) = y(1:ndir)
upar      = y(ndir+1) !< gamma v||
Mr        = y(ndir+2)
!gamma     = y(ndir+3)

call get_vec_rk(igrid_working,x,t_s,b,bp1_,bp^NC_)
call get_vec_rk(igrid_working,x,t_s,e,ep1_,ep^NC_)

call get_vec_rk(igrid_working,x,t_s,bdotgradb,b_dot_grad_b1_,b_dot_grad_b^NC_)
call get_vec_rk(igrid_working,x,t_s,uedotgradb,ue_dot_grad_b1_,ue_dot_grad_b^NC_)
call get_vec_rk(igrid_working,x,t_s,gradkappaB,grad_kappa_B1_,grad_kappa_B^NC_)
call get_vec_rk(igrid_working,x,t_s,bdotgradue,b_dot_grad_ue1_,b_dot_grad_ue^NC_)
call get_vec_rk(igrid_working,x,t_s,uedotgradue,ue_dot_grad_ue1_,ue_dot_grad_ue^NC_)

absb         = sqrt( {^C& b(^C)**2 |+} )
bhat(1:ndir) = b(1:ndir) / absb

epar         = {^C& e(^C)*bhat(^C) |+}

call cross(e,bhat,ue)
ue(1:ndir)   = ue(1:ndir)*CONST_c / absb

kappa = sqrt(1.0d0 - ({^C& ue(^C)**2|+})/CONST_c**2)

gamma = sqrt(1.0d0+upar**2/CONST_c**2+2.0d0*Mr*absb/m/CONST_c**2)/kappa

utmp1(1:ndir) = bhat(1:ndir)/(absb*kappa**2)
utmp2(1:ndir) = Mr*CONST_C/(gamma*q)*gradkappaB(1:ndir) &
 + upar*epar/(gamma*CONST_C)*ue(1:ndir) &
     + m*CONST_C/q* ( upar**2/gamma*bdotgradb(1:ndir) + upar*uedotgradb(1:ndir) &
     + upar*bdotgradue(1:ndir) + gamma*uedotgradue(1:ndir))

call cross(utmp1,utmp2,utmp3)

u(1:ndir) = ue(1:ndir) + utmp3(1:ndir)

!> done assembling the terms, now write rhs:

dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )/ UNIT_LENGTH

dydt(ndir+1) = ({^C& ue(^C) * ( upar*bdotgradb(^C) + gamma*uedotgradb(^C) ) |+}) &
     + q/m*epar - Mr/(m*gamma) * ({^C& bhat(^C)*gradkappaB(^C)|+})


dydt(ndir+2) = 0.0d0 !< magnetic moment is conserved

!dydt(ndir+3) = q/(m*CONST_c**2) * ({^C& dydt(^C)*e(^C)|+}) * UNIT_LENGTH

end subroutine derivs_gca_rk
!=============================================================================
subroutine derivs_gca(t_s,y,dydt)

use mod_odeint
use mod_particles
use mod_gridvars
use mod_particles, only: particle
use constants

include 'amrvacdef.f'

integer                         :: ic^D, igrid_particle, ipe_particle, ipe
double precision                :: t_s, y(*)
double precision                :: dydt(*)
! .. local ..
double precision,dimension(^NC) :: ue, b, e, x, bhat, bdotgradb, uedotgradb, gradkappaB
double precision,dimension(^NC) :: bdotgradue, uedotgradue, u, utmp1, utmp2, utmp3
double precision                :: upar, Mr, gamma, absb, q, m, epar, kappa
!-----------------------------------------------------------------------------
!> Here the normal interpolation is done for the terms in the GCA equations of
!> motion

q = particle(ipart_working)%self%q
m = particle(ipart_working)%self%m

x(1:ndir) = y(1:ndir)
upar      = y(ndir+1) !< gamma v||
Mr        = y(ndir+2)
!gamma     = y(ndir+3)

call get_vec(igrid_working,x,t_s,b,bp1_,bp^NC_)
call get_vec(igrid_working,x,t_s,e,ep1_,ep^NC_)

call get_vec(igrid_working,x,t_s,bdotgradb,b_dot_grad_b1_,b_dot_grad_b^NC_)
call get_vec(igrid_working,x,t_s,uedotgradb,ue_dot_grad_b1_,ue_dot_grad_b^NC_)
call get_vec(igrid_working,x,t_s,gradkappaB,grad_kappa_B1_,grad_kappa_B^NC_)
call get_vec(igrid_working,x,t_s,bdotgradue,b_dot_grad_ue1_,b_dot_grad_ue^NC_)
call get_vec(igrid_working,x,t_s,uedotgradue,ue_dot_grad_ue1_,ue_dot_grad_ue^NC_)

absb         = sqrt( {^C& b(^C)**2 |+} )
bhat(1:ndir) = b(1:ndir) / absb

epar         = {^C& e(^C)*bhat(^C) |+}

call cross(e,bhat,ue)
ue(1:ndir)   = ue(1:ndir)*CONST_c / absb

kappa = sqrt(1.0d0 - ({^C& ue(^C)**2|+})/CONST_c**2)

gamma = sqrt(1.0d0+upar**2/CONST_c**2+2.0d0*Mr*absb/m/CONST_c**2)/kappa

utmp1(1:ndir) = bhat(1:ndir)/(absb*kappa**2)
utmp2(1:ndir) = Mr*CONST_C/(gamma*q)*gradkappaB(1:ndir) &
 + upar*epar/(gamma*CONST_C)*ue(1:ndir) &
     + m*CONST_C/q* ( upar**2/gamma*bdotgradb(1:ndir) + upar*uedotgradb(1:ndir) &
     + upar*bdotgradue(1:ndir) + gamma*uedotgradue(1:ndir)) 

call cross(utmp1,utmp2,utmp3)

u(1:ndir) = ue(1:ndir) + utmp3(1:ndir)

!> done assembling the terms, now write rhs:

dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )/ UNIT_LENGTH

dydt(ndir+1) = ({^C& ue(^C) * ( upar*bdotgradb(^C) + gamma*uedotgradb(^C) ) |+}) &
     + q/m*epar - Mr/(m*gamma) * ({^C& bhat(^C)*gradkappaB(^C)|+})

dydt(ndir+2) = 0.0d0 !< magnetic moment is conserved

!dydt(ndir+3) = q/(m*CONST_c**2) * ({^C& dydt(^C)*e(^C)|+}) * UNIT_LENGTH

end subroutine derivs_gca
!=============================================================================
subroutine set_particles_dt_gca()

use mod_particles
use constants
use mod_gridvars, only: igrid_working, ipart_working
include 'amrvacdef.f'

logical                                   :: BC_applied
integer                         :: ipart, iipart, nout, ic^D, igrid_particle, ipe_particle, ipe
double precision                :: t_min_mype, tout, dt_particles_mype, dt_cfl0, dt_cfl1, dt_a
double precision                :: dxmin, vp, a, gammap
double precision                :: v(1:ndir), y(1:^NC+2),ytmp(1:^NC+2), dydt(1:^NC+2), v0(1:ndir), v1(1:ndir), dydt1(1:^NC+2)
double precision                :: ap0, ap1, dt_cfl_ap0, dt_cfl_ap1, dt_cfl_ap
double precision                :: dt_max_output, dt_max_time, dt_euler, dt_tmp
double precision, parameter     :: cfl=0.8d0, uparcfl=0.8d0 !< make these particle cfl conditions more restrictive if you are interpolating out of the grid
double precision, parameter     :: uparmin=1.0d-6*CONST_C
!-----------------------------------------------------------------------------
!> Here the timestep for the guiding centre integration is chosen

dt_particles      = bigdouble
dt_particles_mype = bigdouble
t_min_mype        = bigdouble
dt_max_output     = bigdouble
dt_max_time       = bigdouble

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

   igrid_working = particle(ipart)%igrid
   ipart_working = particle(ipart)%self%index
  
   dt_tmp = (tmax_particles - particle(ipart)%self%t)
   if (dt_tmp .le. 0.0d0) &
        dt_tmp = smalldouble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> make sure we step only one cell at a time, first check CFL at current location
!> then we make an Euler step to the new location and check the new CFL 
!> we simply take the minimum of the two timesteps. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> added safety factor cfl:
{#IFNDEF D1
   dxmin  = min({rnode(rpdx^D_,particle(ipart)%igrid)})*cfl
}{#IFDEF D1
   dxmin  = rnode(rpdx1_,particle(ipart)%igrid)*cfl
}
   ! initial solution vector:
   y(1:ndir) = particle(ipart)%self%x(1:ndir) !< position of guiding center
   y(ndir+1) = particle(ipart)%self%u(1) !< parallel momentum component (gamma v||)
   y(ndir+2) = particle(ipart)%self%u(2) !< conserved magnetic moment
   ytmp=y
!   y(ndir+3) = particle(ipart)%self%u(3) !< Lorentz factor of guiding centre

   call derivs_gca(particle(ipart)%self%t,y,dydt)
   v0(1:ndir) = dydt(1:ndir)
   ap0        = dydt(ndir+1)

!> guiding center velocity:
   {^C& v(^C)   = abs(dydt(^C))\}
   vp = sqrt({^C& (v(^C))**2|+})

   dt_cfl0    = dxmin/vp
   dt_cfl_ap0 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / ap0)
!  dt_cfl_ap0 = min(dt_cfl_ap0, uparcfl * sqrt(abs(UNIT_LENGTH*dxmin/(ap0+smalldouble))) )

!> make an Euler step with the proposed timestep:
 
  !> new solution vector:
   dt_euler = min(dt_tmp,dt_cfl0,dt_cfl_ap0)
   y(1:ndir+2) = y(1:ndir+2) + dt_euler * dydt(1:ndir+2)
   
   particle(ipart)%self%x(1:ndir) = y(1:ndir) !< position of guiding center
   particle(ipart)%self%u(1)      = y(ndir+1) !< parallel momentum component (gamma v||)
   particle(ipart)%self%u(2)      = y(ndir+2) !< conserved magnetic moment
   
!> first check if the particle is outside the physical domain or in the ghost cells

if (.not. particle_in_igrid(ipart_working,igrid_working)) then
   y(1:ndir+2) = ytmp(1:ndir+2)
end if

   call derivs_gca_rk(particle(ipart)%self%t+dt_euler,y,dydt)
  !call derivs_gca(particle(ipart)%self%t+dt_euler,y,dydt)
   
   v1(1:ndir) = dydt(1:ndir)
   ap1        = dydt(ndir+1)
!> guiding center velocity:
   {^C& v(^C)   = abs(dydt(^C))\}
   vp = sqrt({^C& (v(^C))**2|+})

   dt_cfl1    = dxmin/vp
   dt_cfl_ap1 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / ap1)
  ! dt_cfl_ap1 = min(dt_cfl_ap1, uparcfl * sqrt(abs(UNIT_LENGTH*dxmin/(ap1+smalldouble))) )
  
   dt_tmp = min(dt_euler, dt_cfl1, dt_cfl_ap1)

   particle(ipart)%self%x(1:ndir) = ytmp(1:ndir) !< position of guiding center
   particle(ipart)%self%u(1)      = ytmp(ndir+1) !< parallel momentum component (gamma v||)
   particle(ipart)%self%u(2)      = ytmp(ndir+2) !< conserved magnetic moment
  !      dt_tmp = min(dt_cfl1, dt_cfl_ap1)



! time step due to parallel acceleration:
! The standart thing, dt=sqrt(dx/a) where we comupte a from d(gamma v||)/dt and d(gamma)/dt
!    dt_ap = sqrt(abs(dxmin*UNIT_LENGTH*y(ndir+3)/( dydt(ndir+1) - y(ndir+1)/y(ndir+3)*dydt(ndir+3) ) ) )
!     vp = sqrt({^C& (v(^C)*UNIT_LENGTH)**2|+})
!     gammap = sqrt(1.0d0/(1.0d0-(vp/CONST_c)**2))
!     ap = CONST_c**2/vp*gammap**(-3)*dydt(ndir+3)
!     dt_ap = sqrt(dxmin*UNIT_LENGTH/ap)

!   dt_a = bigdouble
!   if (dt_euler .gt. smalldouble) then 
!      a = sqrt({^C& (v1(^C)-v0(^C))**2 |+})/dt_euler
!      dt_a = min(sqrt(dxmin/a),bigdouble)
!   end if
  

!   particle(ipart)%self%dt = min(dt_tmp , dt_a)
   particle(ipart)%self%dt = dt_tmp


!**************************************************
!> Make sure we don't miss an output or tmax_particles:
!**************************************************
   !> corresponding output slot:
   nout = int(particle(ipart)%self%t/dtsave_ensemble) + 1
   tout = dble(nout) * dtsave_ensemble
   if (particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) then
      dt_max_output = tout - particle(ipart)%self%t
      if (dt_max_output .le. 0.0d0) &
           dt_max_output = smalldouble * tout
   end if

   !> bring to tmax_particles:
   if (particle(ipart)%self%t+particle(ipart)%self%dt .gt. tmax_particles) then
!        particle(ipart)%self%dt = max(tmax_particles - particle(ipart)%self%t , smalldouble * tmax_particles)
        dt_max_time = (tmax_particles - particle(ipart)%self%t)
        if (dt_max_time .le. 0.0d0) &
             dt_max_time = smalldouble * tmax_particles
     end if
!**************************************************
!**************************************************
   particle(ipart)%self%dt = min(particle(ipart)%self%dt,dt_max_time,dt_max_output)

   dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)

   t_min_mype = min(t_min_mype,particle(ipart)%self%t)

end do !< ipart loop

!> keep track of the global minimum:
call MPI_ALLREDUCE(dt_particles_mype,dt_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
     icomm,ierrmpi)

!> keep track of the minimum particle time:
call MPI_ALLREDUCE(t_min_mype,t_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
     icomm,ierrmpi)

end subroutine set_particles_dt_gca
!=============================================================================
}



{#IFDEF PARTICLES_ADVECT
!=============================================================================
subroutine integrate_particles_advect()
!> this solves dx/dt=v for particles

use mod_odeint
use mod_particles
use mod_gridvars, only: get_vec, vp^C_, igrid_working, interpolate_var
include 'amrvacdef.f'
integer                             :: ipart, iipart
double precision                    :: dt_p
double precision, dimension(1:ndir) :: v, x
integer                             :: igrid
double precision                    :: rho, rho1, rho2, td, tloc
double precision, dimension(ixG^T,1:nw)   :: w
!> for odeint:
integer                          :: nok, nbad
double precision                 :: h1
double precision,parameter       :: eps=1.0d-6, hmin=1.0d-8
integer, parameter               :: nvar = ^NC
external derivs_advect
!-----------------------------------------------------------------------------

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);
   
   dt_p = particle(ipart)%self%dt
   igrid = particle(ipart)%igrid
   igrid_working = igrid
   tloc = particle(ipart)%self%t
   x(1:ndir) = particle(ipart)%self%x(1:ndir)

   

   ! **************************************************
   !> Position update
   ! **************************************************

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> Simple forward Euler:
   call get_vec(igrid,x,tloc,v,vp1_,vp^NC_)
   particle(ipart)%self%u(1:ndir) = v(1:ndir)

   particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
        + dt_p * v(1:ndir) &
        / UNIT_LENGTH
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> Adaptive stepwidth RK4:
!   h1 = dt_p/2.0d0
!   call odeint(x,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_advect,rkqs)
!   particle(ipart)%self%x(1:ndir) = x(1:ndir)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




   ! **************************************************
   !> Velocity update
   ! **************************************************

   call get_vec(igrid,x,tloc+dt_p,v,vp1_,vp^NC_)
   particle(ipart)%self%u(1:ndir) = v(1:ndir)





   ! **************************************************
   !> Payload update
   ! **************************************************
   
   !> To give an example, we set the payload to the interpolated density at 
   !> the particles position.  
   !> In general, it would be better to add the auxilary variable to mod_gridvars 
   !> since then you can use the convenient subroutine get_vec() for this.  
   
   
   !if (.not.time_advance) then
   !   
   !   w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
   !   call primitive(ixG^LL,ixG^LL,w,px(igrid)%x)
   !   
   !   call interpolate_var(igrid,ixG^LL,ixM^LL,pw(igrid)%w(ixG^T,rho_),px(igrid)%x(ixG^T,1:ndim),x,rho)

   !else

   !   w(ixG^T,1:nw) = pwold(igrid)%w(ixG^T,1:nw)
   !   call primitive(ixG^LL,ixG^LL,w,px(igrid)%x)
   !   call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,rho_),px(igrid)%x(ixG^T,1:ndim),x,rho1)

   !   w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
   !   call primitive(ixG^LL,ixG^LL,w,px(igrid)%x)
   !   call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,rho_),px(igrid)%x(ixG^T,1:ndim),x,rho2)
   !   
   !   td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt
   !   
   !   rho = rho1 * (1.0d0 - td) + rho2 * td
   !   
   !end if

   !particle(ipart)%self%payload(1) = rho * UNIT_DENSITY


   ! **************************************************
   !> Time update
   ! **************************************************
   particle(ipart)%self%t = particle(ipart)%self%t + dt_p


end do

!=============================================================================
end subroutine integrate_particles_advect
!=============================================================================
subroutine derivs_advect(t_s,x,dxdt)

use mod_gridvars, only: get_vec, vp^C_, igrid_working
include 'amrvacdef.f'

double precision                :: t_s, x(*)
double precision                :: dxdt(*)
! .. local ..
double precision                :: v(1:^NC)
!-----------------------------------------------------------------------------

call get_vec(igrid_working,x(1:^NC),t_s,v,vp1_,vp^NC_)
{^C& dxdt(^C) = v(^C)/ UNIT_LENGTH;}

end subroutine derivs_advect
!=============================================================================
subroutine set_particles_dt_advect()

use mod_particles
include 'amrvacdef.f'

integer                         :: ipart, iipart, nout
double precision                :: t_min_mype, tout, dt_particles_mype, dt_cfl
double precision                :: v(1:ndir)
!-----------------------------------------------------------------------------
dt_particles      = bigdouble
dt_particles_mype = bigdouble
t_min_mype        = bigdouble


do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

!> make sure we step only one cell at a time:
{v(^C)   = abs(particle(ipart)%self%u(^C))\}

{^IFPHI 
!> convert to angular velocity:
if (typeaxial =='cylindrical') v(phi_) = abs(v(phi_)/particle(ipart)%self%x(r_))
}

{#IFNDEF D1
   dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)})
}
{#IFDEF D1
   dt_cfl = rnode(rpdx1_,particle(ipart)%igrid)/v(1)
}

{^IFPHI
if (typeaxial =='cylindrical') then 
!> phi-momentum leads to radial velocity:
if (phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
     sqrt(rnode(rpdx1_,particle(ipart)%igrid)/particle(ipart)%self%x(r_)) &
     / v(phi_))
!> limit the delta phi of the orbit (just for aesthetic reasons):
   dt_cfl = min(dt_cfl,0.1d0/v(phi_))
!> take some care at the axis:
   dt_cfl = min(dt_cfl,(particle(ipart)%self%x(r_)+smalldouble)/v(r_))
end if
}

   particle(ipart)%self%dt = dt_cfl*UNIT_LENGTH


!**************************************************
!> Make sure we don't miss an output or tmax_particles:
!**************************************************
   !> corresponding output slot:
   nout = int(particle(ipart)%self%t/dtsave_ensemble) + 1
   tout = dble(nout) * dtsave_ensemble
   if (particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) &
        particle(ipart)%self%dt = max(tout - particle(ipart)%self%t , smalldouble * tout)

   !> bring to tmax_particles:
   if (particle(ipart)%self%t+particle(ipart)%self%dt .gt. tmax_particles) &
        particle(ipart)%self%dt = max(tmax_particles - particle(ipart)%self%t , smalldouble * tmax_particles)
!**************************************************
!**************************************************


   dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)

   t_min_mype = min(t_min_mype,particle(ipart)%self%t)

end do !< ipart loop

!> keep track of the global minimum:
call MPI_ALLREDUCE(dt_particles_mype,dt_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
     icomm,ierrmpi)

!> keep track of the minimum particle time:
call MPI_ALLREDUCE(t_min_mype,t_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
     icomm,ierrmpi)

end subroutine set_particles_dt_advect
}
!=============================================================================
{#IFDEF PARTICLES_LORENTZ
subroutine integrate_particles_lorentz()
!> this is the relativistic Vay scheme, a leapfrog integrator
use constants
use mod_particles
include 'amrvacdef.f'
integer                             :: ipart, iipart
double precision                    :: lfac, q, m, dt_p, cosphi, sinphi, phi1, phi2, r, re
double precision, dimension(1:ndir) :: b, e, emom, uminus, t_geom, s, udash, tmp, uplus, xcart1, xcart2, ucart2, radmom
double precision, dimension(1:ndir) :: uhalf, tau, uprime, ustar, tovergamma
double precision                    :: lfacprime, sscal, sigma
!-----------------------------------------------------------------------------

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

   q  = particle(ipart)%self%q
   m  = particle(ipart)%self%m
   dt_p = particle(ipart)%self%dt

      ! Push particle over half time step
      call get_lfac(particle(ipart)%self%u,lfac)
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * CONST_C / UNIT_LENGTH

      ! Get E, B at new position
   call get_b(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,b)
   call get_e(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,e)


select case (typeaxial)
! **************************************************
!> CARTESIAN COORDINATES
! **************************************************
case ('slab')
      
! **************************************************
!> Momentum update
! **************************************************


{#IFDEF PARTICLES_VAY
!The Vay mover
        call cross(particle(ipart)%self%u,b,tmp)
        uhalf = particle(ipart)%self%u + &
             q * dt_p /(2.0d0 * m * CONST_C) * (e &
             + tmp/(lfac))
       
	tau = q * dt_p / (2.d0 * m * CONST_C) * b
        uprime = uhalf + q * dt_p / (2.d0 * m * CONST_C) * e
        call get_lfac(uprime,lfacprime)
        sigma = lfacprime**2 - ({tau(^C)*tau(^C)|+})
       
	ustar = ({uprime(^C)*tau(^C)|+}) / CONST_C     
        lfac = sqrt((sigma + sqrt(sigma**2 + 4.d0 * &
             (({tau(^C)*tau(^C)|+}) + ({ustar(^C)*ustar(^C)|+})))) / 2.d0)

        tovergamma = tau / lfac
        sscal = 1.d0 / (1.d0 + ({tovergamma(^C)*tovergamma(^C)|+}))
        call cross(uprime,tovergamma,tmp)
        particle(ipart)%self%u = sscal * (uprime + ({uprime(^C)*tovergamma(^C)|+}) & 
	* tovergamma + tmp)      
}
!The Boris mover
{#IFNDEF PARTICLES_VAY
emom = q * e * dt_p /(2.0d0 * m * CONST_C)

      if (losses) then
         call get_lfac(particle(ipart)%self%u,lfac)
         re = abs(q)**2 / (m * CONST_C**2)
         call cross(particle(ipart)%self%u,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*particle(ipart)%self%u(^C)|+})/lfac)**2) &
              * particle(ipart)%self%u / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if
      
      uminus = particle(ipart)%self%u + emom + radmom
      
      call get_lfac(uminus,lfac)
      call get_t(b,lfac,dt_p,q,m,t_geom)
      call get_s(t_geom,s)
      
      call cross(uminus,t_geom,tmp)
      udash = uminus + tmp
      
      call cross(udash,s,tmp)
      uplus = uminus + tmp

      if (losses) then
         call cross(uplus,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*uplus(^C)|+})/lfac)**2) &
              * uplus / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if

      
      particle(ipart)%self%u = uplus + emom + radmom
}


! **************************************************
!> CYLINDRICAL COORDINATES
! **************************************************
case ('cylindrical')
print*,'This is Boris scheme in cylindrical coordinates'

! **************************************************
!> Momentum update
! **************************************************
      emom = q * e * dt_p /(2.0d0 * m * CONST_C)

      if (losses) then
         call get_lfac(particle(ipart)%self%u,lfac)
         re = abs(q)**2 / (m * CONST_C**2)
         call cross(particle(ipart)%self%u,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*particle(ipart)%self%u(^C)|+})/lfac)**2) &
              * particle(ipart)%self%u / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if
      
      uminus = particle(ipart)%self%u + emom + radmom
      
      call get_lfac(uminus,lfac)
      call get_t(b,lfac,dt_p,q,m,t_geom)
      call get_s(t_geom,s)
      
      call cross(uminus,t_geom,tmp)
      udash = uminus + tmp
      
      call cross(udash,s,tmp)
      uplus = uminus + tmp

      if (losses) then
         call cross(uplus,b,tmp)
         radmom = - third * re**2 * lfac &
              * ( ({(e(^C)+tmp(^C)/lfac)**2|+})  &
              - (({e(^C)*uplus(^C)|+})/lfac)**2) &
              * uplus / m / CONST_C * dt_p
      else
         radmom = 0.0d0
      end if

      
      particle(ipart)%self%u = uplus + emom + radmom

! **************************************************
!> Position update
! **************************************************
      
      !> Get cartesian coordinates:
      phi1        = particle(ipart)%self%x(pphi_)
      cosphi     = cos(phi1)
      sinphi     = sin(phi1)

      xcart1(1)  = particle(ipart)%self%x(r_) * cosphi
      xcart1(2)  = particle(ipart)%self%x(r_) * sinphi
      xcart1(3)  = particle(ipart)%self%x(zz_)

      ucart2(1)   = cosphi * particle(ipart)%self%u(r_) - sinphi * particle(ipart)%self%u(pphi_)
      ucart2(2)   = cosphi * particle(ipart)%self%u(pphi_) + sinphi * particle(ipart)%self%u(r_)
      ucart2(3)   = particle(ipart)%self%u(zz_)

      !> update position
      xcart2(1:ndir) = xcart1(1:ndir) &
           + dt_p * ucart2(1:ndir)/lfac &
           * CONST_C/UNIT_LENGTH


      !> back to cylindrical coordinates
      phi2     = atan2(xcart2(2),xcart2(1))
      if (phi2 .lt. 0.0d0) phi2 = 2.0d0*dpi + phi2
      r       = sqrt(xcart2(1)**2 + xcart2(2)**2)
      particle(ipart)%self%x(r_)   = r
      particle(ipart)%self%x(pphi_) = phi2
      particle(ipart)%self%x(zz_)   = xcart2(3)

! **************************************************
!> Rotate the momentum to the new cooridnates
! **************************************************

      !> rotate velocities
      cosphi     = cos(phi2-phi1)
      sinphi     = sin(phi2-phi1)

      tmp = particle(ipart)%self%u
      particle(ipart)%self%u(r_)   = cosphi * tmp(r_)   + sinphi * tmp(pphi_)
      particle(ipart)%self%u(pphi_) = cosphi * tmp(pphi_) - sinphi * tmp(r_)
      particle(ipart)%self%u(zz_)   = tmp(zz_)

end select

   



! **************************************************
!> Position update over half timestep at the end
! **************************************************
      
      ! update position
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * CONST_C / UNIT_LENGTH

! **************************************************
!> Time update
! **************************************************
   particle(ipart)%self%t = particle(ipart)%self%t + dt_p

! **************************************************
!> Payload update
! **************************************************
      if (npayload > 0) then
      !> current gyroradius
      call cross(particle(ipart)%self%u,b,tmp)
      tmp = tmp / sqrt({b(^C)**2|+})
      particle(ipart)%self%payload(1) = sqrt({tmp(^C)**2|+}) / sqrt({b(^C)**2|+}) * m / abs(q) * 8.9875d+20
      end if
      ! e.b payload
      if (npayload>1) then
      !> e.b (set npayload=2 first):
      particle(ipart)%self%payload(2) = ({e(^C)*b(^C)|+})/ (sqrt({b(^C)**2|+})* sqrt({e(^C)**2|+}))
      end if

end do !< ipart loop

!=============================================================================
contains
!=============================================================================
!> internal procedures
!=============================================================================
subroutine get_t(b,lfac,dt,q,m,t)

use constants
implicit none
double precision, dimension(^NC), intent(in)      :: b
double precision, intent(in)                      :: lfac, dt, q, m
double precision, dimension(^NC), intent(out)     :: t
!-----------------------------------------------------------------------------

t = q * b * dt / (2.0d0 * lfac * m * CONST_C)

end subroutine get_t
!=============================================================================
subroutine get_s(t,s)

implicit none
double precision, dimension(^NC), intent(in)   :: t
double precision, dimension(^NC), intent(out)  :: s
!-----------------------------------------------------------------------------

s = 2.0d0 * t / (1.0d0+{t(^C)**2|+})

end subroutine get_s
!=============================================================================
end subroutine integrate_particles_lorentz
!=============================================================================
subroutine set_particles_dt_lorentz()

use constants
use mod_particles
include 'amrvacdef.f'

integer                         :: ipart, iipart, nout
double precision,dimension(^NC) :: b
double precision                :: lfac,absb,dt_particles_mype,dt_cfl,v(ndir)
double precision                :: t_min_mype, tout
double precision, parameter     :: cfl=0.5d0
!-----------------------------------------------------------------------------
dt_particles      = bigdouble
dt_particles_mype = bigdouble
t_min_mype        = bigdouble

do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);
   
   call get_b(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,b)
   absb = sqrt({b(^C)**2|+})
   call get_lfac(particle(ipart)%self%u,lfac)


!**************************************************
!> CFL timestep
!**************************************************
!> make sure we step only one cell at a time:
{v(^C)   = abs(CONST_C * particle(ipart)%self%u(^C) / lfac)\}

{^IFPHI 
!> convert to angular velocity:
if (typeaxial =='cylindrical') v(phi_) = abs(v(phi_)/particle(ipart)%self%x(r_))
}

{#IFNDEF D1
   dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)})
}
{#IFDEF D1
   dt_cfl = rnode(rpdx1_,particle(ipart)%igrid)/v(1)
}

{^IFPHI
if (typeaxial =='cylindrical') then 
!> phi-momentum leads to radial velocity:
if (phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
     sqrt(rnode(rpdx1_,particle(ipart)%igrid)/particle(ipart)%self%x(r_)) &
     / v(phi_))
!> limit the delta phi of the orbit (just for aesthetic reasons):
   dt_cfl = min(dt_cfl,0.1d0/v(phi_))
!> take some care at the axis:
   dt_cfl = min(dt_cfl,(particle(ipart)%self%x(r_)+smalldouble)/v(r_))
end if
}

dt_cfl = dt_cfl * cfl
!**************************************************
!**************************************************

!> bound by gyro-rotation:
particle(ipart)%self%dt = &
     abs( dtheta * CONST_c/UNIT_LENGTH * particle(ipart)%self%m * lfac &
     / (particle(ipart)%self%q * absb) )

particle(ipart)%self%dt = min(particle(ipart)%self%dt,dt_cfl)*UNIT_LENGTH


!**************************************************
!> Make sure we don't miss an output or tmax_particles:
!**************************************************
   !> corresponding output slot:
   nout = int(particle(ipart)%self%t/dtsave_ensemble) + 1
   tout = dble(nout) * dtsave_ensemble
   if (particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) &
        particle(ipart)%self%dt = max(tout - particle(ipart)%self%t , smalldouble * tout)

   !> bring to tmax_particles:
   if (particle(ipart)%self%t+particle(ipart)%self%dt .gt. tmax_particles) &
        particle(ipart)%self%dt = max(tmax_particles - particle(ipart)%self%t , smalldouble * tmax_particles)
!**************************************************
!**************************************************

   dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)

   t_min_mype = min(t_min_mype,particle(ipart)%self%t)

end do !< ipart loop

!> keep track of the global minimum:
call MPI_ALLREDUCE(dt_particles_mype,dt_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
     icomm,ierrmpi)

!> keep track of the minimum particle time:
call MPI_ALLREDUCE(t_min_mype,t_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
     icomm,ierrmpi)

end subroutine set_particles_dt_lorentz
!=============================================================================
subroutine get_e(igrid,x,tloc,e)

!> Get the electric field in the grid at postion x.
!> For ideal SRMHD, we first interpolate b and u=lfac*v/c
!> The electric field then follows from e = b x beta, where beta=u/lfac.  
!> This ensures for the resulting e that e<b and e.b=0. Interpolating on u 
!> avoids interpolation-errors leading to v>c.  
!>
!> For (non-ideal) MHD, we directly interpolate the electric field as 
!> there is no such constraint.

use mod_particles
use mod_gridvars
include 'amrvacdef.f'

integer,intent(in)                                 :: igrid
double precision,dimension(^NC), intent(in)        :: x
double precision, intent(in)                       :: tloc
double precision,dimension(^NC), intent(out)       :: e

double precision,dimension(^NC)                    :: e1, e2

{^IFSRMHD
double precision,dimension(^NC)                    :: u, u1, u2, beta, b
double precision                                   :: lfac
}
integer                                            :: ic^D
double precision                                   :: td
!-----------------------------------------------------------------------------


{^IFSRMHD

call get_b(igrid,x,tloc,b)

!> get flow four-velocity u:
if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,up^C_),px(igrid)%x(ixG^T,1:ndim),x,u(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,up^C_),px(igrid)%x(ixG^T,1:ndim),x,u1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,up^C_),px(igrid)%x(ixG^T,1:ndim),x,u2(^C))\}

td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt

{^C& u(^C) = u1(^C) * (1.0d0 - td) + u2(^C) * td\}

end if

call get_lfac(u,lfac)

{^C&beta(^C) = u(^C)/lfac\}

!> do the cross-product e = b x beta:

call cross(b,beta,e)
return
}


if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep^C_),px(igrid)%x(ixG^T,1:ndim),x,e2(^C))\}

td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt

{^C& e(^C) = e1(^C) * (1.0d0 - td) + e2(^C) * td\}

end if

end subroutine get_e
!=============================================================================
subroutine get_b(igrid,x,tloc,b)

use mod_particles
use mod_gridvars
include 'amrvacdef.f'

integer,intent(in)                                 :: igrid
double precision,dimension(^NC), intent(in)        :: x
double precision, intent(in)                       :: tloc
double precision,dimension(^NC), intent(out)       :: b
integer                                            :: ic^D
double precision,dimension(^NC)                    :: b1, b2
double precision                                   :: td
!-----------------------------------------------------------------------------

if (.not.time_advance) then

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b(^C))\}

else

{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars_old(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b1(^C))\}
{^C&call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp^C_),px(igrid)%x(ixG^T,1:ndim),x,b2(^C))\}

td = (tloc/(UNIT_LENGTH/UNIT_VELOCITY) - t) / dt

{^C& b(^C) = b1(^C) * (1.0d0 - td) + b2(^C) * td\}

end if


!> flat interpolation:
!{ic^D = int((x(^D)-rnode(rpxmin^D_,igrid))/rnode(rpdx^D_,igrid)) + 1 + dixB \}
!b(1:^NC) = gridvars(igrid)%w(ic^D,bp1_:bp^NC_)

end subroutine get_b
!=============================================================================
subroutine get_lfac(u,lfac)

use constants

double precision,dimension(^NC), intent(in)        :: u
double precision, intent(out)                      :: lfac
!-----------------------------------------------------------------------------

lfac = sqrt(1.0d0 + ({u(^C)**2|+}))

end subroutine get_lfac
}
!=============================================================================
subroutine cross(a,b,c)

include 'amrvacdef.f'
double precision, dimension(^NC), intent(in)   :: a,b
double precision, dimension(^NC), intent(out)  :: c
!-----------------------------------------------------------------------------

!> ndir needs to be three for this to work!!!
{#IFDEF C3
select case (typeaxial)
case ('slab')
   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)
case ('cylindrical')
   c(r_) = a(pphi_)*b(zz_) - a(zz_)*b(pphi_)
   c(pphi_) = a(zz_)*b(r_) - a(r_)*b(zz_)
   c(zz_) = a(r_)*b(pphi_) - a(pphi_)*b(r_)
case default
   call mpistop('geometry not implemented in cross(a,b,c)')
end select
}{#IFNDEF C3
call mpistop("Needs to be run with three components!")
}

end subroutine cross
!=============================================================================
!=============================================================================
}
