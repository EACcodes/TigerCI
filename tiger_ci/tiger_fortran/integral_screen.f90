! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
! all the components need to screen the two-electron integrals using 
! a variety of approaches. 
!
! moved into a module of its own because I need these routines in different 
! cholesky modules and I'm getting tired of dealing with circular dependencies

module integral_screen
  use global_var_mod
  use utilities_mod
  use function_map
  implicit none

  ! Data structure for handling integral screening
  type integral_screen_data
     real(kind=real8), allocatable, dimension(:,:) :: Q_int         ! The Q quantity for CS screening, sqrt[(ij|ij)]
     real(kind=real8), allocatable, dimension(:,:) :: Q_shell      ! The largest sqrt[(ij|ij)] in the Shell pair IJ (i in I, j in J) 
     real(kind=real8) :: tol
     integer :: nSCreenCS
  end type integral_screen_data


contains

  !> \brief This subroutine sets up the data structure used for all integral screening
  !> \param screen The memory for the CS screening quantities
  !> \param nBas  Number of basis functions
  !> \param nShell Number of AO shells
  !> \param diag  The diagonal two electron integrals (ij|ij)
  !> \param bas2shell Maps each AO basis function to its AO shell
  subroutine setup_screening(screen, nBas, nShell, diag, bas2shell)
    implicit none

    ! inputs
    integer, intent(in) :: nBas, nShell
    integer, dimension(:), intent(in) :: bas2shell
    real(kind=real8), dimension(:), intent(in) :: diag
    type(integral_screen_data) :: screen

    ! Hopefully these will get moved to input parameters
    screen%tol   = ao_integral_threshold

    ! Variable initialization
    screen%nSCreenCS  = 0

    call setup_CS_screen(screen, nBas, nShell, diag, bas2shell)
  end subroutine setup_screening




  !> \brief This subroutine sets up the data structure used for Cauchy-Schwarz (CS) screening
  !> \param screen The memory for the CS screening quantities
  !> \param nBas  Number of basis functions
  !> \param nShell Number of AO shells
  !> \param diag  The diagonal two electron integrals (ij|ij)
  !> \param bas2shell Maps each AO basis function to its AO shell
  subroutine setup_CS_screen(screen, nBas, nShell, diag, bas2shell)
    implicit none

    ! inputs
    integer, intent(in) :: nBas, nShell
    integer, dimension(:), intent(in) :: bas2shell
    real(kind=real8), dimension(:), intent(in) :: diag
    type(integral_screen_data) :: screen
    
    ! local variables
    integer :: ierr
    integer :: iBas, jBas
    integer :: iDiag
    integer :: iShell, jShell
    real(kind=real8) :: tmp

    ! allocate memory, if not already allocated 
    if ( .not. allocated(screen%Q_int) ) then
       allocate( screen%Q_int(nBas,nBas), screen%Q_shell(nShell, nShell), stat=ierr)
       call allocatecheck(ierr,"CS integral_screen_data")
    end if
    screen%Q_int = 0.0
    screen%Q_shell = 0.0


    ! First step, building the Q screening quantities at the individual basis set level
    iDiag = 0
    do iBas = 1, nBas
       do jBas = 1, iBas
          iDiag = diag_address(iBas,jBas)
          tmp = sqrt( diag(iDiag) )
          
          screen%Q_int( iBas, jBas ) = tmp
          screen%Q_int( jBas, iBas ) = tmp
       end do
    end do

    ! Now we build the Q screening quantities at the shell level
    do iBas = 1, nBas
       iShell = bas2shell(iBas)
       do jBas = 1, nBas
          jShell = bas2shell(jBas)
          tmp = screen%Q_int(iBas,jBas)
          if ( tmp > screen%Q_shell(iShell,jShell) ) then 
             screen%Q_shell(iShell,jShell) = tmp 
             screen%Q_shell(jShell,iShell) = tmp
          end if
       end do
    end do


  end subroutine setup_CS_screen

  !> \brief Answers the question is this integral shell block screened 
  !> \param data The integral screening data structure
  !> \param iShell  The I in the integral block (IJ|KL) to be screened
  !> \param jShell  The J in the integral block (IJ|KL) to be screened
  !> \param kShell  The K in the integral block (IJ|KL) to be screened
  !> \param lShell  The L in the integral block (IJ|KL) to be screened
  !> \param nFuncInShell  The number of AO functions in each shell
  logical function shellIsScreened(data, iShell, jShell, kShell, lShell, nFuncInShell)
    implicit none
    integer, intent(in) :: iShell, jShell, kShell, lShell
    integer, intent(in) :: nFuncInShell(:)
    type(integral_screen_data) :: data

    shellIsScreened = .false.
    shellIsScreened = isScreenCS(data, iShell, jShell, kShell, lShell, nFuncInShell)

  end function shellIsScreened


  !> \brief Answers the question is this integral shell block screened under CS
  !> \param data The integral screening data structure
  !> \param iShell  The I in the integral block (IJ|KL) to be screened
  !> \param jShell  The J in the integral block (IJ|KL) to be screened
  !> \param kShell  The K in the integral block (IJ|KL) to be screened
  !> \param lShell  The L in the integral block (IJ|KL) to be screened
  !> \param nFuncInShell  The number of AO functions in each shell
  logical function isScreenCS(data, iShell, jShell, kShell, lShell, nFuncInShell)
    implicit none

    ! inputs
    integer, intent(in) :: iShell, jShell, kShell, lShell, nFuncInShell(:)
    type(integral_screen_data) :: data

    ! local variables
    real(kind=real8) :: bound
    integer :: n

    bound = data%Q_shell(iShell,jShell) * data%Q_shell(kShell,lShell)
    isScreenCS = .false.
    if ( bound < data%tol ) then
       iSScreenCS = .true.
       n = nFuncInShell(iShell) * nFuncInShell(jShell) * nFuncInShell(kShell) * nFuncInShell(lShell)
       data%nScreenCS = data%nScreenCS + n
    end if

  end function isScreenCS

  !> \brief Prints out some stats on the number of integrals which weren't calculated due to different screenings
  !> \param data The integral screening data structure
  subroutine print_integral_screen_stats(data)
    implicit none
    type(integral_screen_data) :: data
    write(*,*) "Integrals screened by CS  = ", data%nScreenCS
  end subroutine print_integral_screen_stats

end module integral_screen
