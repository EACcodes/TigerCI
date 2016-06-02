! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module acpf_driver_mod
  use initial_davidson_guess
  use blocked_locks_mod
  use sgga_struct_mod
  use davidson_mod
  use davidson_io_mod
  use natural_mod
  use sgga
  use gtimesc_mod
  use initial_davidson_struct
  use get_basis_data
  implicit none

  private
  public :: acpf_driver

contains



  !> \brief Runs a standard (MR)ACPF/ACPF2/AQCC calculation
  !> \brief Drives any (MR)SDCI calculation with a Davidson-Liu Solver
  !> \param roots      The number of roots to solve for
  !> \param energy     On exit contains the energy of each root
  !> \param sgga_data  Scratch space/data for the SGGA routine 
  subroutine acpf_driver(roots, energy, sgga_data)
    ! inputs
    integer, intent(in) :: roots
    real(kind=real8), dimension(:) :: energy
    type(sgga_struct) :: sgga_data

    ! local variables 
    integer :: ierr
    logical :: ref_ci
    real(kind=real8) :: ref_energy, nuc
    real(kind=real8), dimension(:), allocatable :: work1, work2         
    type(initial_dav_guess) :: guess
    type(blockedLockVectorType) :: civec, sigmavec
    procedure(htimesc), pointer :: hmultc
    procedure(gtimesc), pointer :: gmultc

    ! get the initial wavefunction guess(es). use reference CI calculation for multireference cases
    ref_ci = .false.
    if (num_ref>1) ref_ci = .true.
    call davidson_setup(guess, roots, restart_flag, ref_ci, sgga_data)


    ! produce output 
    write(*,*) " "
    write(*,*) "         ****************************************"
    write(*,*) " "
    write(*,*) "                  (MR)ACPF/ACPF2/AQCC Driver                      "
    write(*,*) " "
    write(*,*) "         **************************************** "     

    ! get the energy for the initial guess
    ! note that when the user supplies the reference energy is comes with the nuc potential, 
    ! when we calculate the reference energy we don't include the nuc potential
    ref_energy = get_reference_energy(guess, 1)
    call get_nuclear_repulsion(nuc)
    if (using_user_reference_energy()) then
       write(*,*) "the reference energy is (elec + nuclear) = " , ref_energy
       ref_energy = ref_energy - nuc
    else
       write(*,*) "the reference energy is (elec + nuclear) = " , ref_energy + nuc
    end if

    ! allocate the needed memory
    call allocateBlockedLockVector(civec, total_csfs, blockSizeCI)
    call allocateBlockedLockVector(sigmavec, total_csfs, blockSizeSigma)
    allocate(work1(total_csfs), work2(total_csfs), stat=ierr)
    call allocatecheck(ierr, "work array in acpf driver")
    civec%v = 0
    sigmavec%v = 0 
    work1 = 0.0
    work2 = 0.0

    ! set up the io (files and what not) for the solver
    call davidson_io_setup()

    ! set my H*c multiplier
    hmultc => htimesc
    gmultc => gtimesc

    ! set up the g values 
    call build_compressed_g(work1)

    ! run the solver
    call davidson_liu_gen_eig(civec, sigmavec, work1, work2, roots, energy, ref_energy, guess, sgga_data, hmultc, gmultc)


  end subroutine acpf_driver



end module acpf_driver_mod
