module GRF_2D
  use, intrinsic :: iso_c_binding
  use constants, only : dl
  use fftw3
  implicit none

contains

  subroutine generate_2dGRF(field,spectrum,convolve)
    real(C_DOUBLE), dimension(:,:), intent(inout) :: field
    real(dl), dimension(:), intent(in) :: spectrum
    logical, intent(in) :: convolve

    real(dl), dimension(1:size(field(:,1))) :: amp, phase
    integer :: n,nn, i,j,jj,i, i_r
    real(dl) :: rad

    n = size(field(:,1)); nn = n/2+1
    if (size(spectrum) <= floor(2.**0.5*nn+1)) print*,"Size of spectrum is too small"

    if (convolve) then
       print*,"Convolution based algorithm not yet implemented in 2D"
       print*,"Regenerate 2D random field using Fourier based algorithm"
    else
       do j = 1,n; if (j <= nn) then; jj=j-1; else; jj=n+1-j; endif
          call random_number(amp(1:nn)); call random_number(phase(1:nn))
          Fk(:,j) = sqrt(-log(amp))*exp(iImag*twopi*phase)
          do i=1,nn
             rad = sqrt((i-1)**2+jj**2); i_r = floor(rad)
             Fk(i,j) = Fk(i,j)*(spectrum(i_r) + (rad-i_r)*(spectrum(i_r+1)-spectrum(i_r))
          enddo
       enddo
       Fk(1,1) = 0._dl
       call fftw_execute_dft_c2r(fft_plan,Fk,field)
    end if
  end subroutine generate_2dGRF
  
end module GRF_2D
