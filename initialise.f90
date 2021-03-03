module initial_conditions
  use constants, only : dl, twopi
  implicit none

  type BubbleParams
     real(dl) :: r0, meff, phit, phif
  end type BubbleParams
  
contains

  subroutine initialise_vacuum(phi,dphi,phi_fv)
    real(dl), dimension(:,:), intent(out) :: phi, dphi
    real(dl), intent(in) :: phi_fv
    integer :: i
    
    dphi = 0._dl
!    call GRF_2D(phi,(/ (1._dl,i=1,256/),.false.)
    phi = phi + phi_fv
  end subroutine initialise_vacuum
  
  !>@brief
  !> Initialise an initial bubble profile using output from an instanton code.
  subroutine initialize_bubble_from_file(phi,dphi,dx)
    real(dl), dimension(:,:), intent(out) :: phi, dphi
    real(dl), intent(in) :: dx

    phi = 0._dl
    dphi = 0._dl
  end subroutine initialize_bubble_from_file

  !>@brief
  !> Initialise a thin-wall bubble profile using a tanh profile
  subroutine initialize_bubble_thin_wall(phi,dphi,dx,r0,meff,phif,phit)
    real(dl), dimension(:,:), intent(out) :: phi, dphi
    real(dl), intent(in) :: dx, r0, meff, phif, phit

    integer :: n, nn, i,j
    real(dl) :: rx2, ry2, r

    n = size(phi(1,:)); nn=n/2
    do j=1,n; ry2 = (j-nn)**2*dx**2
       do i=1,n; rx2 = (i-nn)**2*dx**2
          r = sqrt(rx2+ry2)
          ! phi(i,j) = 0.5_dl*(phif-phit)*tanh(meff*(r-r0))+0.5_dl*(phif+phit)
          phi(i,j) = (phif-phit)*4._dl/twopi*atan(-0.5*exp(meff*r0)/cosh(meff*r)) + phif 
       enddo
    enddo
    dphi = 0._dl
  end subroutine initialize_bubble_thin_wall

  subroutine initialize_gaussian(phi,dphi,dx,amp,sig)
    real(dl), dimension(:,:), intent(out) :: phi, dphi
    real(dl), intent(in) :: dx,amp,sig

    integer :: n,nn,i,j
    real(dl) :: rx2, ry2

    n = size(phi(1,:)); nn=n/2
    do j=1,n; ry2 = (j-nn)**2*dx**2
       do i=1,n; rx2 = (i-nn)**2*dx**2
          phi(i,j) = amp*exp(-0.5*(rx2+ry2)/sig**2)
       enddo
    enddo
    dphi = 0._dl
  end subroutine initialize_gaussian
  
end module initial_conditions
