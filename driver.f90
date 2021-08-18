#include "macros.h"
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

  real(dl) :: l
  integer :: kcut

  l = 64._dl; kcut = 128
  !  call boot_openmp(2)
  !  call set_model_params( (/1.2_dl/) )
  call set_model_params( (/0._dl/) )
  call create_lattice(mySim,256,l,1)

   !call initialize_gaussian(mySim%fld(:,:,1),mySim%dfld(:,:,1),mySim%dx,1._dl,1._dl)
  !call initialize_bubble_thin_wall(mySim%fld(:,:,1),mySim%dfld(:,:,1),mySim%dx,4.,1._dl,0.,2.)

#ifdef ICS  
  open(unit=newunit(u),file='field.bin',access='stream',status='replace')
  open(unit=newunit(u1),file='dfield.bin',access='stream',status='replace')
  do i=1,100
     mySim%fld = 0._dl; mySim%dfld = 0._dl
     call add_vacuum_flucs(mySim,1._dl,-1._dl+lam**2,twopi/128._dl*32)
     write(u) mySim%fld; write(u1) mySim%dfld
  enddo
  close(u); close(u1)
#endif
  
!#ifdef RUN
  call initialize_rand(42,1315)
  open(unit=newunit(u),file='field.bin',access='stream',status='replace')
  open(unit=newunit(u1),file='dfield.bin',access='stream',status='replace')
  do j=1,10
     mySim%fld = -1._dl; mySim%dfld = 0._dl; mySim%time = 0._dl
     !mySim%fld = 0._dl; mySim%dfld = 0._dl; mySim%time = 0._dl
!     call add_vacuum_flucs(mySim,1.5_dl,-1._dl+lam**2,twopi/mySim%lSize*kcut)
!     call add_vacuum_flucs(mySim,2.1_dl,-1._dl+lam**2,twopi/mySim%lSize*kcut)
     call add_vacuum_flucs(mySim,5._dl,2._dl*(1._dl-lam),twopi/mySim%lSize*kcut)
     !call upsample_sim(mySim,2)
     write(u) mySim%fld; write(u1) mySim%dfld
     call output_lattice_mean(mySim,.true.)
     alph = 8._dl
     do i=1,64
        call step_lattice(mySim,mySim%dx/alph,32)
        call output_lattice_mean(mySim)
        write(u) mySim%fld; write(u1) mySim%dfld
     enddo
  enddo
  close(u); close(u1)
!#endif
  
contains

  subroutine scan_ics(mySim,nSamp,phi0,kcut,out_ic)
    type(Lattice), intent(inout) :: mySim
    integer, intent(in) :: nSamp
    real(dl), intent(in) :: phi0, kcut
    logical, intent(in), optional :: out_ic

    real(dl) :: alph  ! move this to an input
    integer :: i,j, u_f, u_df
    logical :: out_

    out_ = .false.; if (present(out_ic)) out_ = out_ic

    if (out_) then
       open(unit=newunit(u_f),file='field_ic.dat',access='stream',status='replace')
       open(unit=newunit(u_df),file='dfield_ic.dat',access='stream',status='replace')
    endif
    do i=1,nSamp
       mySim%fld = 0._dl; mySim%dfld = 0._dl; mySim%time = 0._dl
       call add_vacuum_flucs(mySim,phi0,-1._dl+lam**2,kcut)
       call output_lattice_mean(mySim,.true.)
       if (out_) then; write(u_f) mySim%fld; write(u_df) mySim%dfld; endif
       alph = 4._dl
       do j=1,128
          call step_lattice(mySim,mySim%dx/alph,32)
          call output_lattice_mean(mySim)
          if (out_) then; write(u_f) mySim%fld; write(u_df) mySim%dfld; endif
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
    open(unit=newunit(u),file='field.dat',access='stream',status='replace')
    write(u) fld
    close(u)

    open(unit=newunit(u),file='dfld.dat',access='stream',status='replace')
    write(u) dfld
    close(u)
  end subroutine output_field
  
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
