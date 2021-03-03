program sample_ics
  use constants, only : dl, twopi
  use utils, only : newunit
  use Hamiltonian
  use output
  implicit none

  type(Lattice) :: mySim
  integer :: u_f, u_df, i, nsamp
  
  call create_lattice(mySim,512,32._dl,1)
  nsamp = 1000
  
  open(unit=newunit(u_f),file='field.dat',access='stream')
  open(unit=newunit(u_df),file='dfield.dat',access='stream')
  do i=1,nsamp
     mySim%fld = 0._dl; mySim%dfld = 0._dl
     call add_vacuum_flucs(mySim,1._dl,-1._dl+lam**2)
     write(u_f) mySim%fld
     write(u_df) mySym%dfld
  enddo
  close(u_f); close(u_df)

contains

  ! Copy add_vauum_flucs here
  
end program sample_ics
