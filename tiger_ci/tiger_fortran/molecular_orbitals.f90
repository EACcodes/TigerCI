! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief The Molecular Orbital module. Handles exclusively the MOs
!> \author David Krisiloff
!> \date 9/1/2011
!>
!> Designed in the hopes of making the code just slightly more modular.
!> This should reduced the number of copies of the MOs (we need just 1 !!)
!>


! A couple key points:
! (1) This is in no way compatabile with anything PAO ... mostly because of the headache
!     that rewriting so much of the code ( which I'm not familiar with ) would be. I'm
!     terribly sorry if this causes someone else a headache
!
! (2) Read but DO NOT WRITE TO molecular_orbitals. Unless you really know what you are
!     doing. If you randomly alter the molecular_orbitals, weird things are going to
!     happen. 

module molecular_orbital_mod
    use fortran_timing
    use global_var_mod
    use molecule_var_mod
    use utilities_mod
    use get_integrals 
    use math_utils
    
    implicit none

    real(kind=real8), allocatable, dimension(:,:) :: molecular_orbitals

    contains

    !> \brief Does all the molecular orbital setup, reading them off this and possibly
    !> calculating the LOVOs
    subroutine get_molecular_orbitals()
        implicit none

        ! local variables
        integer :: nOrb                         ! number of molecular orbitals
        integer :: nBas                         ! number of atom centered basis functions
        integer :: ierr
        integer :: iCoeff                                   ! a pointer for MOLCAS

        ! Initialize just for safety
        ierr = 0
        iCoeff = 0
        
        write(ioOutput,130)
    130 format(/,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/&
                ,1x,"!//                              ",/&
                ,1x,"!// ENTERED MOLECULAR ORBITALS   ",/&
                ,1X,"!//                              ",/&
                ,1x,"!//                              ",/&
                ,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/)
        
        
        !> \todo We should have as a user input the name of the orbital file
        nOrb = num_orbitals
        nBas = number_bas

        ! Allocate space for the molecular_orbitals. Notice the typical FORTRAN column ordering.
        ! Each column is a separate orbital
        allocate( molecular_orbitals(nBas, nOrb), stat=ierr )
        call allocatecheck(ierr,'molecular_orbitals')
        molecular_orbitals = 0
        ! Read the MOs from disk (doesn't matter if it's local or nonlocal anymore)
        call get_coefficients(nBas, nOrb, molecular_orbitals)

        ! This next part shouldn't be necessary. However because ever once in a while
        ! something bad happens to (i.e. I screw up) the orbitals file ... lets just double check that
        ! everything is normalized.
        call check_norms(nOrb, nBas)

    end subroutine get_molecular_orbitals

    !> \brief Double checks that all molecular orbitals are normalized to 1
    subroutine check_norms(nOrb, nBas)
        implicit none

        ! input variables
        integer, intent(in) :: nOrb                         ! number of molecular orbitals
        integer, intent(in) :: nBas                         ! number of atom centered basis functions

        ! local variables
        real(kind=real8), allocatable, dimension(:,:) :: AOoverlap
        real(kind=real8) :: norm
        integer :: iOrb, jOrb, ierr
        logical :: error

        ! First allocate and obtain the AOoverlap matrix
        allocate( AOoverlap(nBas, nBas), stat=ierr)
        call allocatecheck(ierr, "AOoverlap")
        call get_overlap(nBas, AOoverlap)
        
        ! Go through each orbital and check the norm
        error = .false.
        write(*,*) "Checking orbital norms ..."
        do iOrb = 1 , nOrb
            norm = calc_norm(molecular_orbitals(:,iOrb), molecular_orbitals(:,iOrb), AOoverlap, nBas )
            if ( abs(norm -1.0) > 5.0D-6 ) then 
               100 format("Orbital ", i5, " has an unusual norm " , E20.10)
               write(*,100) iOrb, norm
               error = .true.
            end if
            do jOrb = iOrb + 1, nOrb
               norm = calc_norm(molecular_orbitals(:,iOrb), molecular_orbitals(:,jOrb), AOoverlap, nBas )
               if ( abs(norm) > 5.0D-6 ) then 
101               format("Orbitals ", i5, " and ", i5, " are not orthogonal, with overlap " , E20.10)
                  write(*,101) iOrb, jOrb, norm
                  error = .true.
               end if
            end do
        end do
        write(*,*) 

        if (error) then 
           write(*,*) "******************************************"
           write(*,*) "                                          "
           write(*,*) "SERIOUS WARNING FROM MOLECULAR_ORBITAL_MOD"
           write(*,*) "NORMALIZATION ISSUES !!!! "
           write(*,*) "                                          "
           write(*,*) "******************************************"
           stop
        else
           write(*,*) "Molecular orbitals are normalized and orthogonal!"
        end if
             
    end subroutine check_norms

    !> \brief Calculates the norm of a molecular orbital, (AO basis is not orthogonal)
    real function calc_norm (vec1, vec2, overlap, nBas )
        implicit none

        ! input variables
        real(kind=real8), dimension(:), intent(in)   :: vec1, vec2
        real(kind=real8), dimension(:,:), intent(in) :: overlap
        integer, intent(in) :: nBas

        ! local variables
        real(kind=real8) :: tmp
        integer :: i,j

        tmp = 0.0
        do i = 1, nBas
            do j = 1, nBas
                tmp = tmp + vec1(i) * vec2(j) * overlap(i,j)
            end do
        end do
        calc_norm = tmp

    end function calc_norm

!#ifndef TIGER_GAMESS
!    !> \brief Writes out the orbital to disk in MOLCAS' format (useful for visualization)
!    !> \param nBas Number of AO orbitals
!    !> \param nOrb Number of MO orbitals to write to disk
!    subroutine write_orbs_MOLCAS(nBas,nOrb)
!        implicit none
!
!        ! inputs
!        integer, intent(in) :: nBas, nOrb
!
!        ! local variables
!        character(len=8) :: fName, Title
!        integer :: unit_num, nSym, ierr
!        integer, dimension(:), allocatable :: indt
!        real(kind=real8), dimension(:), allocatable :: eOrb, Occ
!
!        ! allocate some memory
!        allocate(eOrb(nOrb), Occ(nOrb), indt(nOrb), stat=ierr)
!        call allocatecheck(ierr, "in write_orbs_MOLCAS")
!        eOrb = 0 
!        Occ  = 0 
!        indt = 0 
!
!        !call MOLCAS routine to write the orbitals
!        fName = 'LOVOs'
!        Title = 'debug'
!        unit_num = 768
!        call WrVec( fName, unit_num, 'C', nSym, nBas, nOrb, molecular_orbitals, Occ, eOrb, indt, Title)
!
!        deallocate(Occ, eOrb, indt)
!
!    end subroutine write_orbs_MOLCAs      
!#endif

end module molecular_orbital_mod



