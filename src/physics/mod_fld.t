!> Module for including flux limited diffusion in hydrodynamics simulations
module mod_fld
  implicit none

contains

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine fld_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
      do idim = 1, ndim
        w(ixO^S,iw_mom(idim)) = w(ixO^S,iw_mom(idim)) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,iw_rho)
        if(energy .and. .not.block%e_is_internal) then
          w(ixO^S,iw_e)=w(ixO^S,iw_e) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,iw_mom(idim))
        end if
      end do
    end if

  end subroutine fld_add_source

  subroutine fld_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)

  end subroutine fld_get_dt

end module mod_fld
