! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! At the moment this is just the old ao_int_eval_fly.f code
! placed into a F90 module. Hopefully this is a first step 
! in cleaning up this code. 
!

module cholesky_decomposition
  use global_var_mod
  use utilities_mod
  use math_utils
  use molecule_var_mod
  use io_unit_numbers
  use cholesky_structs
  use integral_screen
  use cholesky_out_of_core
  use get_integrals
  use two_electron_integrals
  use function_map
  use c_io_finterface
  implicit none


contains
  !> \brief Our new AO cholesky decomposition routine. Performs the decomposition using LAPACK for the actual math
  !> \param nShell The number of AO shells
  !> \param nFuncInShell The number of AO basis functions in each AO Shell
  !> \param bas2shell   Maps each basis function (i) to its shell (I)
  !> \param nBas The number of AO basis functions
  !> \param nCho  On exit contains the number of cholesky vectors ( effective rank of the two-electron integrals )
  !> \param threshold  Numerical tolerance (stopping criteria) for the incomplete cholesky decomposition
  !> \param mem      The maximum amount of memory to use (roughly)
  !> \param activePairs     The number of basis products (ij) that are included in the decomposition after prescreening
  subroutine new_cholesky(nShell, nFuncInShell, bas2shell, nBas, nCho, threshold, mem, activePairs)
    implicit none

    ! LAPACK function for the pivoted CD
    external :: dpstrf

    ! inputs
    integer, intent(in) :: nBas, nShell, mem
    integer, intent(out) :: nCho, activePairs
    real(kind=real8), intent(in) :: threshold
    integer, intent(in), dimension(:) :: nFuncInShell, bas2shell

    ! local variables
    integer :: nPairs, disk
    integer :: ierr, columns
    integer, dimension(:), allocatable :: P, disk_pointers 

    real(kind=real8), dimension(:), allocatable :: diag, ooc_diag

    logical, dimension(:,:), allocatable :: prescreened

    type(pairIndex)            :: pairMap
    type(clock)                :: decompTimer
    type(integral_screen_data) :: screen

    ! This cholesky routine is designed to be as simple as possible
    ! We are going to build the two-electron integral matrix (or those
    ! pieces which survive screening) and then let LAPACK deal with 
    ! doing the heavy math.
  
    write(ioOutput,130)
130 format(/,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/&
         ,1x,"!//                              ",/&
         ,1x,"!// ENTERED TWO ELECTRON INTEGRAL",/&
         ,1X,"!//     CHOLESKY DECOMPOSITION   ",/&
         ,1x,"!//                              ",/&
         ,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/)
    call start_clock(decompTimer)

    ! First step we need to get a copy of the diagonal integrals for screening
    nPairs = nBas * (nBas+1) / 2
    allocate(diag(nPairs), stat=ierr)
    call allocatecheck(ierr, "diagonal integrals")
    call calculate_diagonal_integrals(diag, nFuncInShell, nShell)

    ! Next we are going to prescreen. This is going to also require us to setup 
    ! our (ij) --> i,j mapping
    allocate(pairMap%pair2ij(nPairs), pairMap%ij2pair(nBas,nBas), prescreened(nBas,nBas), stat=ierr)
    call allocatecheck(ierr, "pair map and screening")
    call prescreen(diag, pairMap, threshold, activePairs, prescreened, nBas)

    ! Setup for integral screening during calculation of the two-electron integral matrix
    call setup_screening(screen, nBas, nShell, diag, bas2shell)
    
    ! Allocate space for recording the pivoting we'll use
    allocate(P(activePairs), stat=ierr)
    call allocatecheck(ierr, "permutation p")

    ! At this point we know how much memory we are going to need. Lets figure out if we 
    ! can get LAPACK to do the math for us or if we need to get our hands dirty with the ooc method

    ! how many columns can we store in memory
    columns = ooc_mem(activePairs, mem)
    write(*,*) " "
    write(*,*) "=== mem considerations ==="
    write(*,*) "available memory in # reals =  " , mem
    write(*,*) "dimension of the two-electron integrals = " , activePairs
    write(*,*) "number of columns I can store at one time = " , columns

        
    if ( columns < activePairs ) then

       write(*,*) "Insufficient memory to store the prescreened two-electron integral matrix in memory"
       write(*,*) " --- switching to the out-of-core decomposition algorithm --- "

       ! Calculate the two-electron integrals (storing them on disk) 
       allocate(disk_pointers(nBas*nBas), stat=ierr)
       call allocatecheck(ierr, "disk pointers")
       disk = 178       
       disk_pointers = -1 
       call ooc_calculate_two_electron_integrals(nShell, pairMap, prescreened, nFuncInShell, screen, & 
                                                           disk, activePairs, disk_pointers)
                                                           
       ! at this point we could just run the decomposition. However the diag array contains (ij|ij) integrals
       ! we don't need (they were prescreened) we can reduce memory footprint and bookkeeping headache by removing them
       allocate(ooc_diag(activePairs),stat=ierr)
       call allocatecheck(ierr, "ooc_diag ") 
       call get_active_diag(diag, ooc_diag, nBas, prescreened, pairMap)
       deallocate(diag)
                                     
       ! perform the ooc CD
       call ooc_cholesky(columns, ooc_diag, activePairs, threshold, P, disk, disk_pointers, nCho, activePairs)
       
       ! Last step, write the permutation to disk 
       call ooc_write_index(P, pairMap, activePairs)
       
       ! clean up memory
       deallocate(P,pairMap%pair2ij,pairMap%ij2pair,prescreened)
       deallocate(disk_pointers, ooc_diag, screen%Q_int, screen%Q_shell )
       close(unit=disk,status='delete')

    else

       write(*,*) "Already did the decomposition in C++, so I'm just fetching some data now!"
       
       call get_cholesky_data(nCho)
       return
    end if

    ! finish the timing
    write(*,*) " "
    call print_clock(decompTimer, "total two-electron integral cholesky decomposition")

    
  end subroutine new_cholesky


  !> \brief Writes out the cholesky decomposition to disk
  !> \param V The matrix containing the cholesky decomposition
  !> \param P Contains the permutation matrix used for the pivoted decomposition
  !> \param pairMap The data structure mapping the (ij) pairs to i,j
  !> \param nCho The effective rank of A after decomposition
  subroutine write_factorization(V, P, pairMap, activePairs, nCho)

    implicit none

    ! inputs
    integer, intent(in) :: activePairs, nCho, P(:)
    real(kind=real8), intent(in) :: V(:,:)
    integer,allocatable :: pair(:)
    type(pairIndex), intent(in) :: pairMap

    ! local variables
    integer :: i
    integer :: iPair
    
    allocate(pair(2))
    
    ! Write out the permuted pair (ij) basis
    call setup_vector_io_with_c(scratch_directory // 'cvec_index.dat')
    do i = 1, activePairs
       iPair = P(i)
       pair(1)  = pairMap%pair2ij(iPair)%i
       pair(2)  = pairMap%pair2ij(iPair)%j
       call write_ints_with_c(scratch_directory // 'cvec_index.dat', pair, 2, i)
    end do

    deallocate(pair)
    
    ! Write out the cholesky vectors
    call setup_vector_io_with_c(scratch_directory // 'cvec.dat')
    do i = 1, nCho
        call write_vector_with_c(scratch_directory // 'cvec.dat', V(:,i), activePairs,i)
    end do
    
  end subroutine write_factorization

  !> \brief Calculates the lower triangular part of the two-electron integral matrix
  !> \param V The two-electron integral matrix
  !> \param nShell The number of AO shells
  !> \param pairMap The data structure mapping the (ij) pairs to i,j
  !> \param prescreened Data structure containing he (i,j) pairs removed during prescreening
  !> \param nFuncInShell The number of AO basis functions in each AO Shell
  !> \param screen Data structure containing information for screening integrals
  subroutine calculate_two_electron_integrals(V, nShell, pairMap, prescreened, nFuncInShell, screen)
    implicit none

    ! inputs
    integer, intent(in) ::  nFuncInShell(:)
    type(pairIndex) :: pairMap
    real(kind=real8), dimension(:,:), intent(out) :: V
    logical, intent(in) :: prescreened(:,:)
    integer :: nShell
    type(integral_screen_data) :: screen                    ! for integral screening
    
    ! local variables
    integer :: iShell, jShell, kShell, lShell
    integer :: iFunc, jFunc, kFunc, lFunc
    integer :: ij, kl, i, j, k, l
    integer :: ip, jp, kp, lp
    integer :: ijkl, ijShell, klShell
    integer :: ierr
    integer :: intSize
    integer :: tmp

    integer :: my_iShell, my_jShell, my_kShell, my_lShell
    integer :: am_i, am_j, am_k, am_l
    integer :: p_i, p_ij, p_ijk

    real(kind=real8), dimension(:), allocatable :: scratch

    type(clock) :: timer                                    ! timing variable

    ! Start the timer
    call start_clock(timer)

    ! Memory allocations (we need scratch space for the calculations)
    intSize = maxval(nFuncInShell) ** 4
    allocate( scratch(intSize), stat=ierr)
    call allocatecheck(ierr, "scratch space for (ij|kl)")

    ! Alright lets do the hard work. Note that this, without (pre)screening
    ! scales as O(nBas**4)

    ! This is necessary once we start integral screening 
    V = 0.0 

    ! Integral evaluations 
    !***********************      
    ! loop over shells for the eval_ijkl call at each iteration 
    do iShell = 1, nShell
       do jShell= 1, iShell
          do kShell = 1, nShell
             do lShell= 1, kShell
                ijShell = iShell*(iShell+1)/2 + jShell
                klShell = kShell*(kShell+1)/2 + lShell

                ! By permutational symmetry I only need ijShell <= klShell (or the opposite .. just not both)
                if ( ijShell <= klShell ) then

                   ! Integral screening
                   if ( shellIsScreened(screen, iShell, jShell, kShell, lShell, nFuncInShell) ) cycle

                   ! reorder shells for proper erkale angular momentum ordering
                   my_iShell = iShell
                   my_jShell = jShell
                   my_kShell = kShell
                   my_lShell = lShell

                   call funcs_to_am(nFuncInShell(my_iShell), am_i)
                   call funcs_to_am(nFuncInShell(my_jShell), am_j)
                   call funcs_to_am(nFuncInShell(my_kShell), am_k)
                   call funcs_to_am(nFuncInShell(my_lShell), am_l)
                   
                   if (am_i < am_j) then
                      tmp = my_iShell 
                      my_iShell = my_jShell
                      my_jShell = tmp
                   end if
                   
                   if (am_k < am_l) then
                      tmp = my_kShell 
                      my_kShell = my_lShell
                      my_lShell = tmp
                   end if
                   
                   ! evaluate integrals
                   scratch = 0.0
                   call eval_ijkl(my_iShell,my_jShell,my_kShell,my_lShell,scratch,intSize)
                   
                   ! the following are the number of basis functions before this shell
                   ! you can consider them conceptually to be "pointers" to the current position
                   ! in the basis set
                   ip = sum( nFuncInShell(1:my_iShell-1))
                   jp = sum( nFuncInShell(1:my_jShell-1))
                   kp = sum( nFuncInShell(1:my_kShell-1))
                   lp = sum( nFuncInShell(1:my_lShell-1))

                   ! pointer to the TInt array for the (ij|kl) integral
                   ijkl = 0
                   
                   ! loop over all basis functions in the shells we evaluated
                   do iFunc = 1, nFuncInShell(my_iShell)
                      p_i = (iFunc-1) * nFuncInShell(my_jShell)
                      do jFunc = 1, nFuncInShell(my_jShell)
                         p_ij = (p_i + jFunc -1) * nFuncInShell(my_kShell)
                         do kFunc = 1, nFuncInShell(my_kShell)
                            p_ijk = (p_ij + kFunc - 1) * nFuncInShell(my_lShell)
                            do lFunc = 1, nFuncInShell(my_lShell)

                               ! Move the pointer to the next integral
                               ijkl = p_ijk + lfunc
                               
                               ! These are the absolute basis function numbers for this integral
                               i = ip + iFunc
                               j = jp + jFunc
                               k = kp + kFunc
                               l = lp + lfunc
                               
                               ! Avoid (ij) or (kl) if we prescreened them !
                               if ( prescreened(i,j) .or. prescreened(k,l) ) cycle

                               ! Get the positions in the V matrix for (ij) and (kl) 
                               ! this is non-trivial due to the prescreening
                               ij = pairMap%ij2pair(i,j)
                               kl = pairMap%ij2pair(k,l)

                               ! I'm screwing this up slightly ... I need a lower triangular V
                               if ( ij < kl ) then
                                  tmp = kl
                                  kl  = ij 
                                  ij  = tmp
                               end if
                               V(ij,kl) = scratch(ijkl)
                            end do
                         end do
                      end do
                   end do
                end if

             end do
          end do
       end do
    end do
    
    ! clean up 
    deallocate(scratch)
    ! produce a bit of output
    write(*,*) " "
    write(*,*) "=== Calculated the necessary (ij|kl) ==="
    call print_integral_screen_stats(screen)
    call print_clock(timer, "(ij|kl) evaluations")

  end subroutine calculate_two_electron_integrals

  !> \brief Prescreens the product basis (i,j) using the formula from Pedersen. This step is extremely important for performance
  !> \param diag The diagonal integrals (ij|ij)
  !> \param pairMap The data structure mapping the (ij) pairs to i,j
  !> \param threshold  Numerical tolerance (stopping criteria) for the incomplete cholesky decomposition
  !> \param activePairs     The number of basis products (ij) that are included in the decomposition after prescreening
  !> \param prescreened Data structure containing he (i,j) pairs removed during prescreening
  !> \param nBas The number of AO basis functions
  subroutine prescreen(diag, pairMap, threshold, activePairs, prescreened, nBas)
    implicit none

    ! inputs
    integer, intent(in) :: nBas
    real(kind=real8), intent(in) :: threshold, diag(:)
    type(pairIndex) :: pairMap
    integer, intent(out) :: activePairs
    logical, intent(out) :: prescreened(:,:)

    ! local variables
    integer :: iBas, jBas, iDiag, ij
    real(kind=real8) :: max_ijij, test

    !error condition flag
    pairMap%ij2pair = -10000

    ij = 0 
    iDiag = 0
    max_ijij = maxval(diag)
    activePairs = 1
    prescreened = .false.

    test = (threshold ** 2) / max_ijij          ! Prescreening advocated by pedersen et al.

    ! Now we loop over all the diagonal elements and figure out which (if any) we can safely ignore
    ! iDiag points to the diagonal elements 
    ! ij    keeps track of how many diagonal elements we have kept so far
    ! If an element passes screening we add it to the pairMap
    ! otherwise we add it to prescreened
    do iBas = 1, nBas
       do jBas = 1, iBas
          iDiag = diag_address(iBas, jBas)
          if ( diag(iDiag) > test) then 
             ij = ij + 1
             pairMap%pair2ij(ij)%i = iBas
             pairMap%pair2ij(ij)%j = jBas
             pairMap%ij2pair(iBas,jBas) = ij
             pairMap%ij2pair(jBas,iBas) = ij
          else
             prescreened(iBas,jBas) = .true.
             prescreened(jBas,iBas) = .true.
          end if
       end do
    end do

    ! we call the number of diagonal elements (ij|ij) we kept the number of active pairs
    activePairs = ij 

    write(*,*) "=== prescreen results ==="
    write(*,*) "total number of basis pairs    " , nBas * (nBas+1)/2
    write(*,*) "threshold for prescreen        " , test
    write(*,*) "size of the reduced pair basis " , activePairs

  end subroutine prescreen

  !> \brief Calculates the diagonal integrals (ij|ij)
  !> \param diag The diagonal integrals (ij|ij)
  !> \param nFuncInShell The number of AO basis functions in each AO Shell
  !> \param nShell The number of AO shells
  !> \param sortFunc  An external routine required by MOLCAS for sorting integrals after they are calculated
  subroutine calculate_diagonal_integrals(diag, nFuncInShell, nShell)
    implicit none

    ! inputs
    real(kind=8), dimension(:), intent(out) :: diag
    integer, intent(in) :: nShell
    integer, intent(in), dimension(:) :: nFuncInShell

    ! local variables
    integer :: iShell, jShell
    integer :: iBas, jBas
    integer :: intSize
    integer :: ierr
    integer :: i,j
    integer :: ij                  ! pointer for the pair (ij)
    integer :: ijij                ! pointer for the integral (ij|ij)
    integer :: jLimit
    integer :: my_iShell, my_jShell, am_i, am_j, tmp

    real(kind=real8), dimension(:), allocatable :: integral_scratch

    ! Alright to get all the diagonal integrals we calculate all the 
    ! (IJ|IJ) shell batches and grab the individual diagonal integrals (ij|ij)
    ! from those

    intSize = ( maxval(nFuncInShell) ) ** 4
    allocate( integral_scratch(intSize), stat=ierr)

    do iShell = 1, nShell
       do jShell = 1, iShell

          my_iShell = iShell
          my_jShell = jShell

          call funcs_to_am(nFuncInShell(my_iShell), am_i)
          call funcs_to_am(nFuncInShell(my_jShell), am_j)
                   
          if (am_i < am_j) then
             tmp = my_iShell 
             my_iShell = my_jShell
             my_jShell = tmp
          end if
          
          integral_scratch = 0.0 
          call eval_ijkl(my_iShell, my_jShell, my_iShell, my_jShell, integral_scratch, intSize)

          ! We need to place the integrals in integral_scratch into diag
          ! we do this according to Jeremy's original formulas
          do iBas = 0, nFuncInShell(my_iShell)-1
             jLimit = nFuncInShell(my_jShell) -1
             if (my_iShell == my_jShell) jLimit = iBas
             do jBas = 0, jLimit
                i = sum(nFuncInShell(1:my_iShell-1)) + iBas +1
                j = sum(nFuncInShell(1:my_jShell-1)) + jBas +1
                ij = diag_address(i,j)
                ijij = ((iBas*nFuncInShell(my_jShell)+jBas)*nFuncInShell(my_iShell)+iBas)*nFuncInShell(my_jShell)+jBas
                diag(ij) = integral_scratch(ijij+1)
             end do
          end do

       end do ! iShell
    end do     ! jShell
    deallocate(integral_scratch)
  end subroutine calculate_diagonal_integrals

end module cholesky_decomposition

