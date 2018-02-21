module mod_physics_roe

  implicit none
  public

  procedure(sub_average), pointer         :: phys_average => null()
  procedure(sub_get_eigenjump), pointer   :: phys_get_eigenjump => null()
  procedure(sub_rtimes), pointer          :: phys_rtimes => null()

  integer :: nworkroe = -1

  abstract interface
     subroutine sub_average(wL, wR, x, ixmin1,ixmax1, idim, wroe, workroe)
       use mod_global_parameters
       import
       integer, intent(in)             :: ixmin1,ixmax1, idim
       double precision, intent(in)    :: wL(ixGlo1:ixGhi1, nw),&
           wR(ixGlo1:ixGhi1, nw)
       double precision, intent(inout) :: wroe(ixGlo1:ixGhi1, nw)
       double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)
       double precision, intent(in)    :: x(ixGlo1:ixGhi1, 1:1)
     end subroutine sub_average

     subroutine sub_get_eigenjump(wL, wR, wC, x, ixmin1,ixmax1, il, idim,&
         smalla, a, jump, workroe)
       use mod_global_parameters
       import
       integer, intent(in)                          :: ixmin1,ixmax1, il, idim
       double precision, dimension(ixGlo1:ixGhi1, nw)       :: wL, wR, wC
       double precision, dimension(ixGlo1:ixGhi1)           :: smalla, a, jump
       double precision, dimension(ixGlo1:ixGhi1, nworkroe) :: workroe
       double precision, intent(in)                 :: x(ixGlo1:ixGhi1, 1:1)
     end subroutine sub_get_eigenjump

     subroutine sub_rtimes(q, w, ixmin1,ixmax1, iw, il, idim, rq, workroe)
       use mod_global_parameters
       import
       integer, intent(in)             :: ixmin1,ixmax1, iw, il, idim
       double precision, intent(in)    :: w(ixGlo1:ixGhi1, nw),&
           q(ixGlo1:ixGhi1)
       double precision, intent(inout) :: rq(ixGlo1:ixGhi1)
       double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)
     end subroutine sub_rtimes
  end interface

contains

  subroutine phys_roe_check()
    if (.not. associated(phys_average)) call &
       mpistop("Error: no average method has been specified")

    if (.not. associated(phys_get_eigenjump)) call &
       mpistop("Error: no eigenjump method has been specified")

    if (.not. associated(phys_rtimes)) call &
       mpistop("Error: no rtimes method has been specified")
  end subroutine phys_roe_check

end module mod_physics_roe
