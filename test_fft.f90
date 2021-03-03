program test_fft
  use constants, only : dl, twopi
  use utils, only : newunit
#ifdef USEOMP
  use omp_lib
#endif
  use fftw3
  implicit none

  type(transformPair2D) :: tPair
  integer :: n

  integer(kind=8) :: ti, tf, rate  ! For timing
  real(dl) :: tr_i, tr_f
  integer :: nt, i

  !call benchmark_fft_vary_threads(2048,100)
  call test_derivatives()
  
contains

  subroutine benchmark_fft_vary_threads(n,nt)
    integer, intent(in) :: n, nt
    integer :: i,j
    type(transformPair2D) :: tPair

    real(dl), dimension(n,n) :: phi
    integer(kind=8) :: tf,ti,rate
    print*,"Testing thread scaling with ",n,"^2 lattice sites and ",nt," transforms"
    call boot_openmp()
#ifdef USEOMP

    call system_clock(ti)
    do j=1,nt
       call random_number(phi)
    enddo
    call system_clock(tf,rate)
    print*,"Time to generate ",nt," random fields: ",dble(tf-ti)/rate

    call random_number(phi)
    
    do i=1,omp_get_max_threads()
       call fftw_plan_with_nthreads(i)
       call initialize_transform_2d(tPair,(/n,n/))

       print*,"Using ",i," threads"
       call system_clock(ti)
       do j=1,nt
          tPair%realSpace = phi
          call laplacian_2d_wtype(tPair,1._dl)
       enddo
       call system_clock(tf,rate)
       print*,"Time for ",nt," Laplacians :",dble(tf-ti)/rate
       call system_clock(ti)
       do j=1,nt
          tPair%realSpace = phi
          call forward_transform_2d_wtype(tPair)
       enddo
       call system_clock(tf,rate)
       print*,"Time to perform",nt," forward transforms: ",dble(tf-ti)/rate
       call destroy_transform_2d(tPair)
    enddo
#else
    print*,"Threads not enabled"
#endif
  end subroutine benchmark_fft_vary_threads
  
  subroutine benchmark_fft_vary_n()
    integer :: nlat
    type(transformPair2D) :: tPair

    integer :: i
    do i=0,4
       call initialize_transform_2d(tPair,(/nlat,nlat/))
    enddo
    
  end subroutine benchmark_fft_vary_n
  
  !>@brief
  !> Test the Laplacian and gradient squared operator
  subroutine test_derivatives()
      call test_2d_laplacian(128,16._dl)
      call test_2d_laplacian(256,16._dl)
      call test_2d_laplacian(256,32._dl)

      call test_2d_grad_squared(128,16._dl)
      call test_2d_grad_squared(256,16._dl)
      call test_2d_grad_squared(256,32._dl)
  end subroutine test_derivatives
  
  subroutine test_2d_laplacian(n,len)
    integer, intent(in) :: n
    real(dl), intent(in) :: len
    type(transformPair2D) :: tPair
    integer :: i,j,u
    real(dl) :: ry2, rx2, dx, dk
    real(dl), dimension(1:n,1:n) :: g_lap
    
    dk = twopi/len; dx = len/dble(n)
    
    do j=1,n; ry2 = (j-n/2)**2*dx**2; do i=1,n; rx2 = (i-n/2)**2*dx**2
       g_lap(i,j) = exp(-0.5_dl*(ry2+rx2))*(rx2+ry2-2._dl)
    enddo; enddo
    call initialize_transform_2d(tPair,(/n,n/))
    do j=1,n; ry2 = (j-n/2)**2*dx**2
       do i=1,n; rx2 = (i-n/2)**2*dx**2
          tPair%realSpace(i,j) = exp(-0.5_dl*(ry2+rx2))
       enddo
    enddo
    call laplacian_2d_wtype(tPair,dk)

    open(unit=newunit(u),file='laplacian_slice.dat')
    j=n/2
    do i=1,n
       write(u,*) (i-n/2)*dx, tPair%realSpace(i,j), g_lap(i,j)
    enddo
    close(u)

    print*,"Maximum error in Laplacian is ",maxval(abs(tPair%realSpace-g_lap))," with dx = ",dx," and L = ",len
  end subroutine test_2d_laplacian

  subroutine test_2d_grad_squared(n,len)
    integer, intent(in) :: n
    real(dl), intent(in) :: len
    type(transformPair2D) :: tPair
    integer :: i,j,u
    real(dl) :: ry2, rx2, dx, dk
    real(dl), dimension(1:n,1:n) :: g_grad2
    
    dk = twopi/len; dx = len/dble(n)
    
    do j=1,n; ry2 = (j-n/2)**2*dx**2; do i=1,n; rx2 = (i-n/2)**2*dx**2
       g_grad2(i,j) = exp(-(ry2+rx2))*(rx2+ry2)
    enddo; enddo
    call initialize_transform_2d(tPair,(/n,n/))
    do j=1,n; ry2 = (j-n/2)**2*dx**2
       do i=1,n; rx2 = (i-n/2)**2*dx**2
          tPair%realSpace(i,j) = exp(-0.5_dl*(ry2+rx2))
       enddo
    enddo
    print*,"dk is ",dk
    call grad_squared_2d_wtype(tPair,dk)

    print*,"Maximum error in (Grad f)^2 is ",maxval(abs(tPair%realSpace-g_grad2))," with dx = ",dx," and L = ",len

    open(unit=newunit(u),file='grad2_slice.dat')
    j=n/2
    do i=1,n
       write(u,*) (i-n/2)*dx, tPair%realSpace(i,j), g_grad2(i,j), tPair%realSpace(n/2,i)
    enddo
    close(u)

  end subroutine test_2d_grad_squared
    
  ! Add subroutine for testing a sine-wave
  
  real(dl) function gaussian(r,amp,sig)
    real(dl), intent(in) :: r,amp,sig

    gaussian = amp*exp(-0.5*r**2/sig**2)
  end function gaussian
  
  real(dl) function gaussian_lap(r,amp,sig)
    real(dl), intent(in) :: r,amp,sig

    gaussian_lap = (amp/sig**2)*exp(-0.5*r**2/sig**2)*(r**2/sig**2-2._dl)
  end function gaussian_lap
    
end program test_fft
