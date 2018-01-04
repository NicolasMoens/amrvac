!> Module for including flux limited diffusion in hydrodynamics simulations
module mod_fld
  implicit none

  !> Set FLD mod_variables:
  !double precision :: fld_kappa = 0.34

contains

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine fld_add_source(qdt,ixI^L,ixO^L,wCT,w,x)

    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: radiation_force(ixI^S,1:ndim)
    double precision :: radiation_cooling(ixI^S)
    double precision :: radiation_heating(ixI^S)

    radiation_force(ixI^S,1:ndim) = wCT(ixI^S,iw_r_f(:)) !need c and kappa and r_f
    radiation_cooling(ixI^S) = wCT(ixI^S,iw_r_e) !need boltzman and temp?
    radiation_heating(ixI^S) = wCT(ixI^S,iw_r_e) !need speed  of light and r_e
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!     THESE TERMS ARE JUST TEMPORARY DUMMIES     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! momentum equation source term
    w(ixO^S,iw_mom) = w(ixO^S,iw_mom) &
        + qdt * radiation_force(ixO^S,iw_r_f)

    ! energy equation source terms
    w(ixO^S,iw_e)=w(ixO^S,iw_e) &
        + qdt * radiation_heating(ixI^S) &
        - qdt * radiation_cooling(ixI^S)
  end subroutine fld_add_source

!  subroutine fld_adv_radparam(qdt,ixI^L,ixO^L,wCT,w,x,)

!    use mod_global_parameters
!    use mod_usr_methods

!    integer, intent(in)             :: ixI^L, ixO^L
!    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
!    double precision, intent(in)    :: wCT(ixI^S,1:nw)
!    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! radiation flux

    ! radiation energy

!  end subroutine fld_add_source

!  subroutine fld_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)

!  end subroutine fld_get_dt

end module mod_fld
