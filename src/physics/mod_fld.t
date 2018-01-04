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

    double precision :: radiation_force(ixI^S,ndim)
    double precision :: radiation_cooling(ixI^S)
    double precision :: radiation_heating(ixI^S)

    radiation_force(ixI^S,ndim) = wCT(ixI^S,r_f) !need c and kappa and r_f
    radiation_cooling(ixI^S) = wCT(ixI^S,r_e) !need boltzman and temp?
    radiation_heating(ixI^S) = wCT(ixI^S,r_e) !need speed  of light and r_e

      !do idim = 1, ndim
    ! momentum equation
    w(ixO^S,iw_r_f(:)) = w(ixO^S,,iw_r_f(:)) &
        + qdt * radiation_force(ixO^S,iw_r_f(:))
      !end do
    ! energy equation
    w(ixO^S,iw_e)=w(ixO^S,iw_e) &
        + qdt * radiation_heating(ixI^S,iw_e) &
        - qdt * radiation_cooling(ixO^S,iw_e)
  end subroutine fld_add_source

  subroutine fld_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)

  end subroutine fld_get_dt

end module mod_fld
