module Model
  use constants, only : dl, twopi
  implicit none

  real(dl) :: lam = 1.2_dl

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
       write(fNum,*) "# Model is: V = cos(f) + 0.5*(",lam,")^2*sin(f)^2 - 1"
       write(fNum,*) "# Parameter lambda = ",lam
    endif
  end subroutine write_model_header

  
  elemental function v(phi)
    real(dl), intent(in) :: phi
    real(dl) :: v
    v = cos(phi) + 0.5_dl*lam**2*sin(phi)**2 - 1._dl
  end function v

  elemental function vp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vp
    vp = -sin(phi) + 0.5_dl*lam**2*sin(2._dl*phi)
  end function vp

  elemental function vpp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vpp
    vpp = -cos(phi) + lam**2*cos(2._dl*phi)
  end function vpp

end module Model
