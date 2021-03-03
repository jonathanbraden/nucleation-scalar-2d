#define SPECTRAL 1
module Hamiltonian
  use constants, only : dl, twopi
  use fftw3  ! For Laplacian, etc
  implicit none

  real(dl), parameter :: lam = 1.2_dl
  
  ! Add macros for finite-differencing stencils
  ! Add options to pad grid
  
  type Lattice
     real(dl), dimension(:,:,:), allocatable :: fld, dfld
     real(dl) :: time
     integer :: nlat, nFld
     real(dl) :: dx, lSize, dk
!#ifdef SPECTRAL
     type(transformPair2D) :: tPair
!#endif
  end type Lattice
  
  real(dl), dimension(:,:,:), allocatable :: fld, dfld
  real(dl) :: time
  integer :: nlat
  real(dl) :: dx, lSize, dk
!#ifdef SPECTRAL  
  type(transformPair2D) :: tPair
!#endif

  integer, parameter :: n_Hamiltonian_terms = 2  ! This is a bit ugly
  
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
    integer :: i

!    do i=1,this%nFld
!       this%tPair%realSpace(:,:) = this%fld(:,:,i)
!       call laplacian_2d(this%nLat,this%nLat,this%tPair%realSpace,this%tPair%specSpace,this%dk)
!       this%dfld(:,:,i) = this%dfld(:,:,i) + dt*this%tPair%realSpace(:,:)
    !    enddo

    this%tPair%realSpace(:,:) = this%fld(:,:,1)
    call laplacian_2d_wtype(this%tPair,this%dk)

    this%dfld(:,:,1) = this%dfld(:,:,1) + dt*this%tPair%realSpace(:,:)
    this%dfld = this%dfld - dt*vp(this%fld)
  end subroutine Hamiltonian_potential

  elemental function v(phi)
    real(dl), intent(in) :: phi
    real(dl) :: v
    !    v = 0.5_dl*phi**2
    v = cos(phi) + 0.5_dl*lam**2*sin(phi)**2 - 1._dl
  end function v

  elemental function vp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vp
    ! vp = phi
    vp = -sin(phi) + 0.5_dl*lam**2*sin(2._dl*phi)
  end function vp

  elemental function vpp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vpp
    ! vpp = 1._dl
    vpp = -cos(phi) + lam**2*cos(2._dl*phi)
  end function vpp
  
end module Hamiltonian
