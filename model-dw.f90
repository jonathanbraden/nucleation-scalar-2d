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
       write(fNum,*) "# Model is: V = 0.25(f^2-1)^2 + eps*(f^3/3 - f + 2/3)"
       write(fNum,*) "# Parameter eps = ",lam
    endif
  end subroutine write_model_header

  elemental function v(phi)
    real(dl), intent(in) :: phi
    real(dl) :: v
    v = 0.25*(phi**2-1._dl)**2 + lam*(phi**3/3._dl-phi+2._dl/3._dl)
  end function v

  elemental function vp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vp
    vp = (phi**2-1._dl)*(phi+lam)
  end function vp

  elemental function vpp(phi)
    real(dl), intent(in) :: phi
    real(dl) :: vpp
    vpp = 3._dl*phi**2 - 1._dl + 2._dl*lam*phi
  end function vpp
  
end module Model
