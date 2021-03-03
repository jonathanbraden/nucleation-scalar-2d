module output
  use utils, only : newunit
  use Hamiltonian

contains

  subroutine init_output()    
  end subroutine init_output
  
  subroutine output_lattice(this)
    type(Lattice), intent(inout) :: this
  end subroutine output_lattice

  subroutine output_lattice_mean(this,new_sim)
    type(Lattice), intent(inout) :: this
    logical, intent(in), optional :: new_sim
    
    integer, save :: u_mean
    logical :: o
    real(dl) :: vol

    inquire(opened=o,file='log.out')
    if (.not.o) then
       open(unit=newunit(u_mean),file='log.out')
       write(u_mean,*) "# Time  <Phi>  <DPhi>  <KE>  <GE>  <PE>"
    endif

    if (present(new_sim).and.new_sim) write(u_mean,*)
    
    vol = dble(this%nlat)*dble(this%nlat)
    this%tPair%realSpace = this%fld(:,:,1)
    call grad_squared_2d_wtype(this%tPair,this%dk)

    write(u_mean,*) this%time, sum(this%fld(:,:,1))/vol, sum(this%dfld(:,:,1))/vol, 0.5_dl*sum(this%dfld(:,:,1)**2)/vol, 0.5_dl*sum(this%tPair%realSpace)/vol, sum(v(this%fld))/vol, sum(cos(this%fld(:,:,1)))/vol
  end subroutine output_lattice_mean
    
end module output
