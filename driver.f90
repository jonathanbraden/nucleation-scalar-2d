program bubble_nucleation_2d
  use constants, only : dl, twopi
  use utils, only : newunit
  use Hamiltonian
  use output
  use initial_conditions
  use yoshida
  implicit none
  type(Lattice) :: mySim

  integer :: u, i, u1, j
  real(dl) :: alph

  call boot_openmp(2)
  call create_lattice(mySim,512,64._dl,1)

#ifdef RESAMPLE_TEST
  mySim%fld = 0._dl; mySim%dfld = 0._dl
  call add_vacuum_flucs(mySim,1._dl,-1._dl+lam**2,twopi/128._dl*32)
  
  open(unit=newunit(u),file='field_base.dat',access='stream',status='replace'); write(u) mySim%fld; close(u)
  open(unit=newunit(u),file='dfield_base.dat',access='stream',status='replace'); write(u) mySim%dfld; close(u)

  call upsample_sim(mySim,2)
  open(unit=newunit(u),file='field_up.dat',access='stream',status='replace'); write(u) mySim%fld; close(u)
  open(unit=newunit(u),file='dfield_up.dat',access='stream',status='replace'); write(u) mySim%dfld; close(u)
#endif

   !call initialize_gaussian(mySim%fld(:,:,1),mySim%dfld(:,:,1),mySim%dx,1._dl,1._dl)
  !call initialize_bubble_thin_wall(mySim%fld(:,:,1),mySim%dfld(:,:,1),mySim%dx,4.,1._dl,0.,2.)

#ifdef ICS  
  open(unit=newunit(u),file='field.dat',access='stream',status='replace')
  open(unit=newunit(u1),file='dfield.dat',access='stream',status='replace')
  do i=1,100
     mySim%fld = 0._dl; mySim%dfld = 0._dl
     call add_vacuum_flucs(mySim,1._dl,-1._dl+lam**2,twopi/128._dl*32)
     write(u) mySim%fld; write(u1) mySim%dfld
  enddo
  close(u); close(u1)
#endif
  
!#ifdef RUN
  call initialize_rand(42,1315)
  do j=1,1000
  mySim%fld = 0._dl; mySim%dfld = 0._dl; mySim%time = 0._dl
  call add_vacuum_flucs(mySim,2.12_dl,-1._dl+lam**2,twopi/64._dl*128)
!  call upsample_sim(mySim,2)
!  open(unit=newunit(u),file='field.dat',access='stream',status='replace')
!  open(unit=newunit(u1),file='dfield.dat',access='stream',status='replace')
!  write(u) mySim%fld; write(u1) mySim%dfld
  call output_lattice_mean(mySim,.true.)
  alph = 4._dl
  do i=1,64
     call step_lattice(mySim,mySim%dx/alph,32)
!     write(u) mySim%fld; write(u1) mySim%dfld
     call output_lattice_mean(mySim)
  enddo
  enddo
!  close(u)
!#endif
  
contains

  subroutine scan_ics(mySim,nSamp,phi0,kcut)
    type(Lattice), intent(inout) :: mySim
    integer, intent(in) :: nSamp
    real(dl), intent(in) :: phi0, kcut

    integer :: i,j
    real(dl) :: alpha
    
    do i=1,nSamp
       mySim%fld = 0._dl; mySim%dfld = 0._dl; mySim%time = 0._dl
       call add_vacuum_flucs(mySim,phi0,-1._dl+lam**2,kcut)
       call output_lattice_mean(mySim,.true.)
       alph = 4._dl
       do j=1,128
          call step_lattice(mySim,mySim%dx/alph,32)
          call output_lattice_mean(mySim)
       enddo
    enddo
  end subroutine scan_ics
  
  subroutine test_integrator()
    integer :: i,j,u
    call create_lattice(mySim,2,16._dl,1)

    open(unit=newunit(u),file='field.dat')
    do j=0,8
       mySim%fld = 1._dl; mySim%dfld = 0._dl; mySim%time = 0._dl
       do i=1,32
          call step_lattice(mySim,twopi/8._dl/2.**j,2**j)
          write(u,*) mySim%time, mySim%fld(1,1,1), mySim%dfld(1,1,1)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine test_integrator
  
  subroutine output_field(fld,dfld)
    real(dl), dimension(:,:), intent(in) :: fld, dfld
    integer :: u, i,j,n

    n = size(fld(1,:))
    open(unit=newunit(u),file='field.dat',access='stream')
    write(u) fld
    close(u)

    open(unit=newunit(u),file='dfld.dat',access='stream')
    write(u) dfld
    close(u)
  end subroutine output_field

  ! This currently deals with the Nyquist by setting it to zero
  subroutine upsample_sim(this,ds)
    type(Lattice), intent(inout) :: this
    integer, intent(in) :: ds

    complex(C_DOUBLE_COMPLEX), dimension(1:mySim%nlat/2+1,1:mySim%nlat,1:2) :: Fk
    integer :: i
    real(dl) :: len
    integer :: n, nf, nn, n_new

    n = this%nlat; nf = this%nFld; len = this%lSize
    nn = n/2 ! exclude Nyquist
    n_new = n*ds
    
    ! Fix this to work with multiple fields
    this%tPair%realSpace = this%fld(:,:,1)
    call forward_transform_2d_wtype(this%tPair)
    Fk(:,:,1) = this%tPair%specSpace / dble(n)**2

    this%tPair%realSpace = this%dfld(:,:,1)
    call forward_transform_2d_wtype(this%tPair)
    Fk(:,:,2) = this%tPair%specSpace / dble(n)**2
    call destroy_lattice(this)
    call create_lattice(this,n_new,len,nf)

    ! Add correct normalization (1/ds**2) I think
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
  
  subroutine add_vacuum_flucs(this,phi0,m2eff,kcut)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: phi0, m2eff,kcut

    real(dl), dimension(1:this%nlat/2+1) :: amp, phase
    real(dl) :: norm, rad2
    integer :: n,nn, i,j,jj
    complex(C_DOUBLE_COMPLEX), parameter :: iImag = (0._dl,1._dl)

    !! Boot random number generator
    n = this%nlat; nn = n/2+1
    
    norm = 1._dl/this%lSize/2.**0.5/phi0
    do j=1,n; if (j <= nn) then; jj=j-1; else; jj=n+1-j; endif
       call random_number(amp(1:nn)); call random_number(phase(1:nn))
       this%tPair%specSpace(:,j) = sqrt(-log(amp))*exp(iImag*twopi*phase)
       do i=1,nn
          rad2 = ((i-1)**2+jj**2)*this%dk**2
          if (rad2 > kcut**2) then
             this%tPair%specSpace(i,j) = 0._dl
          else
             this%tPair%specSpace(i,j) = norm*this%tPair%specSpace(i,j)/(rad2+m2eff)**0.25
          endif
       enddo
    enddo
    this%tPair%specSpace(1,1) = 0._dl
    call fftw_execute_dft_c2r(this%tPair%planb,this%tPair%specSpace,this%tPair%realSpace)
    this%fld(:,:,1) = this%fld(:,:,1) + this%tPair%realSpace(:,:)

    
    do j=1,n; if (j <= nn) then; jj=j-1; else; jj=n+1-j; endif
       call random_number(amp(1:nn)); call random_number(phase(1:nn))
       this%tPair%specSpace(:,j) = sqrt(-log(amp))*exp(iImag*twopi*phase)
       do i=1,nn
          rad2 = ((i-1)**2+jj**2)*this%dk**2
          if (rad2 > kcut**2) then
             this%tPair%specSpace(i,j) = 0._dl
          else
             this%tPair%specSpace(i,j) = norm*this%tPair%specSpace(i,j)*(rad2+m2eff)**0.25
          endif
       enddo
    enddo
    this%tPair%specSpace(1,1) = 0._dl
    call fftw_execute_dft_c2r(this%tPair%planb,this%tPair%specSpace,this%tPair%realSpace)
    this%dfld(:,:,1) = this%tPair%realSpace(:,:)
  end subroutine add_vacuum_flucs

! Warning, this subroutine isn't currently threadsafe
  subroutine initialize_rand(seed, seedfac)
    integer, intent(in) :: seed, seedfac
    integer :: nseed, i
    integer, allocatable, dimension(:) :: seeds

    call random_seed(size=nseed)
    print*,"Seed size is ",nseed
    allocate(seeds(1:nseed))
    seeds = seed + seedfac*(/ (i-1, i=1,nseed) /)
    call random_seed(put=seeds)
    deallocate(seeds)
  end subroutine initialize_rand

  
end program bubble_nucleation_2d
