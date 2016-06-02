! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module sdci_driver_mod
  use initial_davidson_guess
  use blocked_locks_mod
  use sgga_struct_mod
  use davidson_mod
  use davidson_io_mod
  use natural_mod
  use sgga   
  use aposteriori_mod
  use get_basis_data
  
  implicit none

  private
  public :: sdci_driver 

contains

  !> \brief Drives any (MR)SDCI calculation with a Davidson-Liu Solver
  !> \param roots      The number of roots to solve for
  !> \param energy     On exit contains the energy of each root
  !> \param sgga_data  Scratch space/data for the SGGA routine 
  subroutine sdci_driver(roots, energy, sgga_data)
    ! inputs
    integer, intent(in) :: roots
    real(kind=real8), dimension(:) :: energy
    type(sgga_struct) :: sgga_data

    ! local variables 
    type(initial_dav_guess) :: guess
    type(blockedLockVectorType) :: civec, sigmavec
    procedure(htimesc), pointer :: mult
    logical :: ref_ci
    real(kind=real8) :: nuc, ref_energy

    ! get the initial wavefunction guess(es). use reference CI calculation for multireference cases
    ref_ci = .false.
    if (num_ref>1) ref_ci = .true.
    call davidson_setup(guess, roots, restart_flag, ref_ci, sgga_data)       

    ! print out the reference energy (useful for debugging)
    ref_energy = get_reference_energy(guess, 1)
    call get_nuclear_repulsion(nuc)
    write(*,*) "the reference energy is (elec + nuclear) = " , ref_energy +nuc

    ! produce output 
    write(*,*) " "
    write(*,*) "         ****************************************"
    write(*,*) " "
    write(*,*) "                    (MR)SDCI Driver                      "
    write(*,*) " "
    write(*,*) "         **************************************** "     

    ! allocate the needed memory
    call allocateBlockedLockVector(civec, total_csfs, blockSizeCI)
    call allocateBlockedLockVector(sigmavec, total_csfs, blockSizeSigma)
    civec%v = 0
    sigmavec%v = 0 

    ! set up the io (files and what not) for the solver
    call davidson_io_setup()

    ! set my H*c multiplier
    mult => htimesc

    ! run the solver
    call davidson_liu(civec, sigmavec, roots, energy, guess, sgga_data, mult)

    ! calculate a posteriori corrections
    if (roots==1) call calculate_a_posteriori_corrections(civec%v, guess, energy(1), get_reference_energy(guess,1))

    ! if we want natural orbitals calculate them now
    if ( (nat_orb_flag == 1) .and. (roots==1) ) then
       call natural_orbital_driver(civec)
    end if
  end subroutine sdci_driver

end module sdci_driver_mod

