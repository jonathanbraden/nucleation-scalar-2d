#include "macros.h"
module Simulation
  use constants, only : dl, twopi
  use fftw3
  implicit none
  
  type Lattice
     real(dl), dimension(:,:,:), allocatable :: fld, dfld
     real(dl) :: time
     integer :: nlat, nFld
     real(dl) :: dx, lSize, dk
!#ifdef SPECTRAL
     type(transformPair2D) :: tPair
!#endif
  end type Lattice

#ifdef ANISOTROPIC
  real(dl), parameter :: c2=0._dl, c1=1._dl, c0=-4._dl
#else
  real(dl), parameter :: c2=1._dl/6._dl, c1=2._dl/3._dl, c0=-10._dl/3._dl
#endif

  
contains

  !>@brief
  !> Create a lattice object.
  !
  !>@param[out] this - The Lattice object to create
  !>@param[in]  n    - Number of lattice sites per side
  !>@param[in]  len  - Side length of the square
  !>@param[in]  nf   - Number of fields living on the lattice
  subroutine create_lattice(this,n,len,nf)
    type(Lattice), intent(out) :: this
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: len

    this%time = 0._dl
    this%nlat = n; this%lSize = len
    this%dx = len/dble(n); this%dk = twopi/len
    call initialize_transform_2d(this%tPair,(/n,n/))
#ifdef SPECTRAL
    allocate( this%fld(1:n,1:n,1),this%dfld(1:n,1:n,1) )
#else
    allocate( this%fld(0:n+1,0:n+1,1),this%dfld(0:n+1,0:n+1,1) )
#endif
    
  end subroutine create_lattice

  subroutine destroy_lattice(this)
    type(Lattice), intent(inout) :: this
    this%time = -1._dl; this%nlat = -1; this%lSize = -1.
    this%dx = -1.; this%dk = -1.; this%nFld = -1
    if (allocated(this%fld)) deallocate(this%fld)
    if (allocated(this%dfld)) deallocate(this%dfld)
    call destroy_transform_2d(this%tPair)
  end subroutine destroy_lattice

  subroutine write_lattice_header(this,fNum)
    type(Lattice), intent(in) :: this
    integer, intent(in) :: fNum
    logical :: o

    inquire(opened=o,unit=fNum)
    if (o) then
       write(fNum,*) "# Lattice parameters are :"
       write(fNum,*) "#  nlat = ",this%nlat,", side length = ",this%lSize,", num fields = ",this%nfld
#ifdef SPECTRAL
       write(fNum,*) "# Fourier spectral derivatives"
#else
       write(fNum,*) "# Finite-differencing derivatives"
#endif
    endif
  end subroutine write_lattice_header

  !>@brief
  !> Refine the fields in the simulation volume to increase precision.
  !> Note that currently this routine deals with the Nyquist modes by zeroing them out, so it should only be used on resolved fields.
  !>
  !> ONLY TESTED FOR ds = 2^p.
  !
  !>@param[in] this - Lattice object to upsample
  !>@param[in] ds   - The upsample factor (only tested with ds = 2^p)
  subroutine upsample_sim(this,ds)
    type(Lattice), intent(inout) :: this
    integer, intent(in) :: ds

    complex(C_DOUBLE_COMPLEX), dimension(1:this%nlat/2+1,1:this%nlat,1:2) :: Fk
    integer :: i
    real(dl) :: len
    integer :: n, nf, nn, n_new

    n = this%nlat; nf = this%nFld; len = this%lSize
    nn = n/2 ! exclude Nyquist
    n_new = n*ds
    
    this%tPair%realSpace = this%fld(:,:,1)
    call forward_transform_2d_wtype(this%tPair)
    Fk(:,:,1) = this%tPair%specSpace / dble(n)**2

    this%tPair%realSpace = this%dfld(:,:,1)
    call forward_transform_2d_wtype(this%tPair)
    Fk(:,:,2) = this%tPair%specSpace / dble(n)**2
    call destroy_lattice(this)
    call create_lattice(this,n_new,len,nf)

    this%tPair%specSpace = 0._dl
    this%tPair%specSpace(1:nn,1:nn) = Fk(1:nn,1:nn,1)
    this%tPair%specSpace(1:nn,n_new-nn:n_new) = Fk(1:nn,n-nn:n,1)
    call fftw_execute_dft_c2r(this%tPair%planb,this%tPair%specSpace,this%tPair%realSpace)
    this%fld(:,:,1) = this%tPair%realSpace

    this%tPair%specSpace = 0._dl
    this%tPair%specSpace(1:nn,1:nn) = Fk(1:nn,1:nn,2)
    this%tPair%specSpace(1:nn,n_new-nn:n_new) = Fk(1:nn,n-nn:n,2)
    call fftw_execute_dft_c2r(this%tPair%planb,this%tPair%specSpace,this%tPair%realSpace)
    this%dfld(:,:,1) = this%tPair%realSpace
  end subroutine upsample_sim
  
  !>@brief
  !> Compute the gradient energy density averaged over the lattice using spectral derivative stencils.
  !> Comparing the results with direct True and False is a test of integration by parts on the lattice
  !
  !>@param[in] Lattice object containing simulation data
  !>@param[in] direct : If .true. evaluates sum (grad phi)^2 directly
  !>                    If .false. -sum(phi grad^2phi)
  real(dl) function gradient_energy_spectral(this,direct) result(ge)
    type(Lattice), intent(inout) :: this
    logical, intent(in) :: direct
    integer :: n
    
    n = this%nlat
    this%tPair%realSpace(1:n,1:n) = this%fld(1:n,1:n,1)
    if (direct) then
       call grad_squared_2d_wtype(this%tPair,this%dk)
       ge = 0.5_dl*sum(this%tPair%realSpace(1:n,1:n))/dble(n)**2
    else
       call laplacian_2d_wtype(this%tPair,this%dk)
       ge = -0.5_dl*sum(this%fld(1:n,1:n,1)*this%tPair%realSpace(1:n,1:n))/dble(n)**2
    endif
  end function gradient_energy_spectral

  !>@brief
  !> Same as gradient_energy_spectral, except using finite-differencing derivative stencils.
  real(dl) function gradient_energy_discrete(this,direct) result(ge)
    type(Lattice), intent(inout) :: this
    logical, intent(in) :: direct
    integer :: n
    n = this%nlat
    call wrap_field(this)

    if (direct) then
       ge = 0.5_dl
    else
       print*,"discrete gradient energy not yet tested"
!       ge = -0.5_dl*sum(this%fld(:,:,1)*STENCIL(c,LAPLACIAN))
    endif
  end function gradient_energy_discrete
  
  subroutine wrap_field(this)
    type(Lattice), intent(inout) :: this
    integer :: n
#ifdef SPECTRAL
#else
    n = mySim%nlat
    this%fld(0,:,:) = this%fld(n,:,:); this%fld(n+1,:,:) = this%fld(1,:,:)
    this%fld(:,0,:) = this%fld(:,n,:); this%fld(:,n+1,:) = this%fld(:,1,:)
#endif
  end subroutine wrap_field
  
end module Simulation
