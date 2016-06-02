! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief Sets up the one electron integrals in the MO basis on disk.
!> \author David Krisiloff
!> \date Jan 2012

module one_electron_integrals
    use global_var_mod
    use molecule_var_mod
    use utilities_mod
    use get_integrals 
    use math_utils

    contains

    !> \brief Handles the one electron integrals for the CI calculation.
    !>
    !> First reads the integrals in the AO basis off of disk. Then it transforms them to the MO basis. 
    !> Finally they are written to disk again for later use.
    subroutine one_electron_driver( molecular_orbitals )
        implicit none

        real(kind=real8), dimension(:,:), intent(in) :: molecular_orbitals

        real(kind=real8), dimension(:), allocatable :: one_ints_AO         ! one electron integrals in the AO basis
        real(kind=real8), dimension(:), allocatable :: one_ints_MO         ! one electron integrals in the MO basis
        real(kind=real8), dimension(:), allocatable :: mol_orb_1D          ! molecular_orbitals packed into a 1D array
        real(kind=real8), dimension(:), allocatable :: workspace           ! a workspace for the AO -> MO transformation

        integer :: nao,nmo                                                 ! Size of the one electron integrals                    
        integer :: ierr                                                    ! return code
        integer :: iOrb, iBas, k                                           ! loop variables

        ! First step allocate the arrays
        nmo = num_orbitals * (num_orbitals + 1) / 2
        nao = number_bas * (number_bas + 1) / 2
        allocate(one_ints_AO(nao), one_ints_MO(nmo), &
                 mol_orb_1D(num_orbitals*number_bas), workspace(number_bas), stat = ierr)
        call allocatecheck(ierr, "one_ints ")
        one_ints_AO = 0.0
        one_ints_MO = 0.0
        mol_orb_1D  = 0.0
        workspace   = 0.0
        
        ! Ask the underlying qchem program for the one electron integrals (AO basis)
        call get_one_ints(one_ints_AO, number_bas, nao)

        ! Now we transform from the AO to MO basis using some old f77 code
        ! (this could always be changed later)
        ! first we pack the molecular_orbitals into a 1D array
        k = 0
        do iOrb = 1 , num_orbitals
            do iBas = 1, number_bas
                k = k + 1
                mol_orb_1D(k) = molecular_orbitals(iBas, iOrb)
            end do
        end do
        

        ! now call spxfrm
        call spxfrm( mol_orb_1D, one_ints_AO, one_ints_MO, workspace, number_bas, num_orbitals)

        ! write the one electron integrals (MO basis) to disk
        write(ioii) one_ints_MO

        ! Clean up
        deallocate(one_ints_AO, one_ints_MO, mol_orb_1D, workspace,stat=ierr)
        call deallocatecheck(ierr, "one_ints ")
    end subroutine one_electron_driver

end module one_electron_integrals
