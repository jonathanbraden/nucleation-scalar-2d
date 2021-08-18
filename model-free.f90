module Model
  use constants, only : dl, twopi
  implicit none

  real(dl) :: lam = 0._dl
  
contains

  subroutine set_model_params(par)
    real(dl), dimension(1), intent(in) :: par
    lam = par(1)
  end subroutine set_model_params

  subroutine write_model_header(fNum)
    integer, intent(in) :: fNum
    logical :: o

    inquire(opened=o,unit=fNum)
    if (o) then
       write(fNum,*) "# Model is: V = 0.5 f^2"
       write(fNum,*) "# Parameter m^2 = 1"
    endif
  end subroutine write_model_header

  elemental function v(phi)
    real(dl), intent(in) :: phi
    real(dl) :: v
    v = 0.5_dl*phi**2
  end function v

  elemental function vp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vp
    vp = phi
  end function vp

  elemental function vpp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vpp
    vpp = 0._dl
  end function vpp
  
end module Model
