!>   The module add viscous source terms and check time step
!>
!>   Viscous forces in the momentum equations:
!>   d m_i/dt +=  - div (vc_mu * PI)
!>   !! Viscous work in the energy equation:
!>   !! de/dt    += - div (v . vc_mu * PI)
!>   where the PI stress tensor is
!>   PI_i,j = - (dv_j/dx_i + dv_i/dx_j) + (2/3)*Sum_k dv_k/dx_k
!>   where vc_mu is the dynamic viscosity coefficient (g cm^-1 s^-1). 
module mod_viscosity
  implicit none

  !> Viscosity coefficient
  double precision, public :: vc_mu = 1.d0

  !> fourth order
  logical :: vc_4th_order = .false.

  !> source split or not
  logical :: vc_split= .false.

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

contains
  !> Read this module"s parameters from a file
  subroutine vc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /vc_list/ vc_mu, vc_4th_order, vc_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, vc_list, end=111)
111    close(unitpar)
    end do

  end subroutine vc_params_read

  !> Initialize the module
  subroutine viscosity_init()
    use mod_global_parameters
    integer :: nwx,idir

    call vc_params_read(par_files)

    ! Determine flux variables
    nwx = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mom(idir) = nwx       ! momentum density
    end do

    nwx = nwx + 1
    e_     = nwx          ! energy density

  end subroutine viscosity_init

  subroutine viscosity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x,energy,qsourcesplit,active)
  ! Add viscosity source to w within ixO 
    use mod_global_parameters
    
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim), wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    
    integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim,idir,jdir,iw
    double precision:: lambda(ixImin1:ixImax1,ixImin2:ixImax2,ndir,ndir),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2),v(ixImin1:ixImax1,ixImin2:ixImax2,ndir)

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources 
      if(.not.vc_4th_order) then
        ! involves second derivatives, two extra layers
        ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;
        if( ixImin1>ixmin1 .or. ixImax1<ixmax1.or. ixImin2>ixmin2 .or. &
           ixImax2<ixmax2)call mpistop&
           ("error for viscous source addition, 2 layers needed")
        ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
      else
        ! involves second derivatives, four extra layers
        ixmin1=ixOmin1-4;ixmin2=ixOmin2-4;ixmax1=ixOmax1+4;ixmax2=ixOmax2+4;
        if( ixImin1>ixmin1 .or. ixImax1<ixmax1.or. ixImin2>ixmin2 .or. &
           ixImax2<ixmax2)&
          call mpistop("error for viscous source addition"//&
          "requested fourth order gradients: 4 layers needed")
        ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;
      end if

      ! get velocity
      do idir=1,ndir
        v(ixImin1:ixImax1,ixImin2:ixImax2,idir)=wCT(ixImin1:ixImax1,&
           ixImin2:ixImax2,mom(idir))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           rho_)
      end do

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      lambda(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir,1:ndir)=zero
      
      !next construct 
      do idim=1,ndim; do idir=1,ndir
      ! Calculate velocity gradient tensor within ixL: gradv= grad v, 
      ! thus gradv_ij=d_j v_i
        tmp(ixImin1:ixImax1,ixImin2:ixImax2)=v(ixImin1:ixImax1,ixImin2:ixImax2,&
           idir)
        select case(typegrad)
        case("central")
          call gradient(tmp,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
             ixmax1,ixmax2,idim,tmp2)
        case("limited")
          call gradientS(tmp,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
             ixmax1,ixmax2,idim,tmp2)
        end select
        lambda(ixmin1:ixmax1,ixmin2:ixmax2,idim,idir)= lambda(ixmin1:ixmax1,&
           ixmin2:ixmax2,idim,idir)+ tmp2(ixmin1:ixmax1,ixmin2:ixmax2)
        lambda(ixmin1:ixmax1,ixmin2:ixmax2,idir,idim)= lambda(ixmin1:ixmax1,&
           ixmin2:ixmax2,idir,idim)+ tmp2(ixmin1:ixmax1,ixmin2:ixmax2)
      enddo; enddo;
      
      ! Multiply lambda with viscosity coefficient and dt
      lambda=lambda*vc_mu*qdt
      
      !calculate div v term through trace action separately
      tmp=0.d0
      do idir=1,ndir
         tmp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp(ixmin1:ixmax1,&
            ixmin2:ixmax2)+lambda(ixmin1:ixmax1,ixmin2:ixmax2,idir,idir)
      end do
      tmp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp(ixmin1:ixmax1,ixmin2:ixmax2)/3.d0
      
      !substract trace from diagonal elements
      do idir=1,ndir
         lambda(ixmin1:ixmax1,ixmin2:ixmax2,idir,idir)=lambda(ixmin1:ixmax1,&
            ixmin2:ixmax2,idir,idir)-tmp(ixmin1:ixmax1,ixmin2:ixmax2)
      enddo
      
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr) 
      ! hence m_j=m_j+d_i tensor_ji
      do idir=1,ndir
        do idim=1,ndim 
              tmp(ixmin1:ixmax1,ixmin2:ixmax2)=lambda(ixmin1:ixmax1,&
                 ixmin2:ixmax2,idir,idim)
              select case(typegrad)
              case("central")
                call gradient(tmp,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,idim,tmp2)
              case("limited")
                call gradientS(tmp,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,idim,tmp2)
              end select
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,mom(idir))+tmp2(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
        enddo
      end do

      if(energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        do idim=1,ndim
          tmp=0.d0
          do idir=1,ndir
             tmp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp(ixmin1:ixmax1,&
                ixmin2:ixmax2)+v(ixmin1:ixmax1,ixmin2:ixmax2,&
                idir)*lambda(ixmin1:ixmax1,ixmin2:ixmax2,idir,idim)
          end do
          call gradient(tmp,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,idim,tmp2)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        enddo
      end if
    end if

  end subroutine viscosity_add_source

  subroutine viscosity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
  
  ! Check diffusion time limit for dt < dtdiffpar * dx**2 / (mu/rho)
  
  use mod_global_parameters
  
  integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2
  double precision, intent(in) :: dx1,dx2, x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(inout) :: dtnew
  
  double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)
  double precision:: dtdiff_visc, dxinv2(1:ndim)
  integer:: idim
  
  ! Calculate the kinematic viscosity tmp=mu/rho
  
  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=vc_mu/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     rho_)
  
  dxinv2(1)=one/dx1**2;dxinv2(2)=one/dx2**2;
  do idim=1,ndim
     dtdiff_visc=dtdiffpar/maxval(tmp(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)*dxinv2(idim))
     ! limit the time step
     dtnew=min(dtnew,dtdiff_visc)
  enddo
  
  end subroutine viscosity_get_dt

end module mod_viscosity
