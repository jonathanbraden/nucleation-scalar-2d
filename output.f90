#include "macros.h"
module output
  use utils, only : newunit
  use simulation
  use Hamiltonian

contains

  subroutine output_lattice(this)
    type(Lattice), intent(inout) :: this
  end subroutine output_lattice

  subroutine output_lattice_mean(this,new_sim)
    type(Lattice), intent(inout) :: this
    logical, intent(in), optional :: new_sim
    
    integer, save :: u_mean
    logical :: o
    real(dl) :: vol
    integer :: n
    real(dl), dimension(2) :: means
    
    inquire(opened=o,file='log.out')
    if (.not.o) then
       open(unit=newunit(u_mean),file='log.out')
       call write_model_header(u_mean)
       call write_lattice_header(this,u_mean)
       write(u_mean,*) "# Time  <Phi>  <DPhi>  <KE>  <GE>  <PE> <GE>_parts <cos(phi)> <V'> <Phi^2> <DPhi^2>"
    endif

    if (present(new_sim).and.new_sim) write(u_mean,*)
    n = this%nlat
    vol = dble(this%nlat)*dble(this%nlat)

    means(1) = sum(this%fld(1:n,1:n,1))/vol; means(2) = sum(this%dfld(1:n,1:n,1))/vol
    
    write(u_mean,*) this%time, sum(this%fld(1:n,1:n,1))/vol, sum(this%dfld(1:n,1:n,1))/vol, 0.5_dl*sum(this%dfld(:,:,1)**2)/vol, gradient_energy_spectral(this,.true.), sum(v(this%fld))/vol, gradient_energy_spectral(this,.false.), sum(cos(this%fld(1:n,1:n,1)))/vol, sum(vp(this%fld(1:n,1:n,1)))/vol, sum((this%fld(1:n,1:n,1)-means(1))**2)/vol, sum((this%dfld(1:n,1:n,1)-means(2))**2)/vol
  end subroutine output_lattice_mean
    
end module output
