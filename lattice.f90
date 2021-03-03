#define SPECTRAL
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
    allocate( this%fld(1:n,1:n,1),this%dfld(1:n,1:n,1) )
!#ifdef SPECTRAL
    call initialize_transform_2d(this%tPair,(/n,n/))
!#endif
  end subroutine create_lattice

  subroutine destroy_lattice(this)
    type(Lattice), intent(inout) :: this
    this%time = -1._dl; this%nlat = -1; this%lSize = -1.
    this%dx = -1.; this%dk = -1.; this%nFld = -1
    if (allocated(this%fld)) deallocate(this%fld)
    if (allocated(this%dfld)) deallocate(this%dfld)
    call destroy_transform_2d(this%tPair)
  end subroutine destroy_lattice

  
end module Simulation
