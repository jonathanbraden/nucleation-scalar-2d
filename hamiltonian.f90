#include "macros.h"
module Hamiltonian
  use constants, only : dl, twopi
  use fftw3  ! For Laplacian, etc
  use Model  ! Contains potential, etc.
  use simulation
  implicit none

  integer, parameter :: n_Hamiltonian_terms = 2  
 
contains
  
  subroutine Hamiltonian_Split(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call Hamiltonian_kinetic(this,dt)
    case (2)
       call Hamiltonian_potential(this,dt)
    case default
       print*,"Undefined Hamiltonian term"
       stop
    end select
  end subroutine Hamiltonian_Split

  subroutine symp_o2_step(this,dt,w1,w2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2

    integer :: i

    do i=2,n_Hamiltonian_terms-1
       call Hamiltonian_Split(this,0.5_dl*w1*dt,i)
    enddo
    call Hamiltonian_Split(this,w1*dt,n_Hamiltonian_terms)
    do i=n_Hamiltonian_terms-1,2,-1
       call Hamiltonian_Split(this,0.5_dl*w1*dt,i)
    enddo
    call Hamiltonian_Split(this,0.5_dl*(w1+w2)*dt,1)
  end subroutine symp_o2_step
  
  subroutine Hamiltonian_kinetic(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    this%fld = this%fld + dt*this%dfld
  end subroutine Hamiltonian_kinetic

  subroutine Hamiltonian_potential(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i, n
    real(dl), dimension(1:this%nlat,1:this%nfld) :: LAP, DV
    
    n = this%nlat
!    do i=1,this%nFld
!       this%tPair%realSpace(:,:) = this%fld(:,:,i)
!       call laplacian_2d(this%nLat,this%nLat,this%tPair%realSpace,this%tPair%specSpace,this%dk)
!       this%dfld(:,:,i) = this%dfld(:,:,i) + dt*this%tPair%realSpace(:,:)
    !    enddo

#ifdef SPECTRAL
    this%tPair%realSpace(1:n,1:n) = this%fld(1:n,1:n,1)
    call laplacian_2d_wtype(this%tPair,this%dk)
    this%dfld(1:n,1:n,1) = this%dfld(1:n,1:n,1) + dt*this%tPair%realSpace(1:n,1:n)
    this%dfld(1:n,1:n,:) = this%dfld(1:n,1:n,:) - dt*vp(this%fld(1:n,1:n,:))
#else
    call wrap_field(this)
    do i=1,n
       LAP = STENCIL(c,LAPLACIAN)
       this%dfld(:,1:n,:) = this%dfld(1:n,1:n,1) !- dt*STENCIL(c,LAPLACIAN)
       this%dfld(:,1:n,:) = this%dfld(1:n,1:n,1) - dt*vp(this%fld(1:n,1:n,1))
    enddo
#endif
  end subroutine Hamiltonian_potential
  
end module Hamiltonian
