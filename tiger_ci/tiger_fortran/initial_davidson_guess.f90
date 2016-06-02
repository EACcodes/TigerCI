! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief This module computes an initial guess for the eigenvalue solver

module initial_davidson_guess
  use global_var_mod
  use molecule_var_mod
  use utilities_mod
  use sort_utils
  use c_sort_finterface
  use sgga_struct_mod
  use io_unit_numbers
  use ci_utilities_mod
  use blocked_locks_mod
  use davidson_io_mod
  use sgga
  use davidson_mod
  use initial_davidson_struct
  use fortran_timing  
  use locist_var_mod
  use graph_var_mod
  use allowed_virtuals_mod
  use graph_mod
  use diag_element_mod
  use find_references_mod
  implicit none

  private
  ! note that some of these are actually in initial_davidson_struct ...
  public :: initial_dav_guess, davidson_setup, pull_initial_guess, get_reference_energy, &
       deallocate_initial_dav_guess

contains

  !> \brief Sets up an initial guess for a Davidson-type solver
  !> \param guess The initial guess data structure
  !> \param roots The number of roots which require an initial guess
  !> \param restart  Logical. If true generates an initial guess from an old wavefunction on disk 
  !> \param ref_ci   Logical. If true generates an initial guess from a reference CI calculation
  !> \param sgga_data Scratch data for the SGGA code
  !> \param filter (Optional) A function pointer to a filter if the pool of diagonal H elements should be filtered before sampling for initial guesses 
  !>
  !> If (restart==.false.) and (ref_ci==.false.) we get the initial guess by looking
  !> at the CSFs with the lowest energies (\f$ \langle i | H | i \rangle \f$)
  recursive subroutine davidson_setup(guess, roots, restart, ref_ci, sgga_data, filter)
    ! inputs
    type(initial_dav_guess) :: guess
    integer, intent(in) :: roots
    logical, intent(in) :: restart, ref_ci 
    type(sgga_struct) :: sgga_data
    procedure(filter_interface), pointer, optional :: filter

    ! local variables
    integer :: i, ierr


    write(*,*) " "
    write(*,*) "  Starting initial vector selection for a Davidson type algorithm"

    ! first handle restarting from disk 
    if (restart) then 
       guess%from_disk = .true.
       guess%unit_num = ioCIFinal
       write(*,*) "  We are starting from a previous calculation ... initial vectors are on disk"
    else
       guess%from_disk = .false.
       guess%unit_num = kill_number
       ! we need to generate the initial guess ourselves
       ! first find the memory
       allocate(guess%vecs(roots), stat=ierr)
       call allocatecheck(ierr, "guess%vecs")
       do i = 1, roots
          call create_sparse_vector(guess%vecs(i),roots)
       end do
       write(*,*) "  New calculation. Packing initial vector into sparse storage"
       allocate(guess%energies(roots), stat=ierr)
       call allocatecheck(ierr, "guess%energies")
       guess%energies = 0

       ! Don't run a ref_ci guess for a ref_ci calculation 
       ! (besides being overkill it won't work ... I manually reset reference_ci_flag == 0 after
       ! the ref_ci guess)
       if (ref_ci .and. .not. (reference_ci_flag==1) ) then 
          write(*,*) "  Grabbing initial vector/energies from a reference CI"
          call reference_ci(guess%vecs, guess%energies, roots, sgga_data)
       else
          write(*,*) "  Grabbing initial vector/energies from lowest few <i|H|i> "
          call init_from_diag(guess%vecs, guess%energies, roots, filter)
       end if
    end if
  end subroutine davidson_setup

  !> \brief Produces \f$n\f$ initial guess from the lowest few diagonal elements of the Hamiltonian
  !> \param vectors  The initial vectors (sparse)
  !> \param energies The energies of each initial vector
  !> \param n        The number of guess to produce
  !> \param filter   (Optional) A function pointer to filter out a subspace of diagonal elements before sorting
  subroutine init_from_diag(vectors, energies, n, filter)
    ! inputs
    type(sparse_vector), dimension(:) :: vectors
    real(kind=real8), dimension(:) :: energies
    integer, intent(in) :: n 
    procedure(filter_interface), pointer :: filter

    ! local variables
    integer :: i, ierr
    integer, dimension(:), allocatable :: index
    real(kind=real8), dimension(:), allocatable :: diag, lowest_few

    ! read in the diagonal elements 
    allocate(diag(total_csfs), stat=ierr)
    call allocatecheck(ierr, "diag elements in init_from_diag")
    rewind(ioDiagStore)
    call read_big_vector(diag, total_csfs, ioDiagStore)

    ! build an index for sorting and a secondary array to store just the lowest few
    allocate(lowest_few(n), index(n), stat=ierr)
    call allocatecheck(ierr, "lowest_few and index")
    
    ! find the lowest few elements
    call find_lowest_n(diag, lowest_few, index, filter)

    ! set up the initial guesses
    do i = 1, n 
       call add_to_sparse_vector(1.0, index(i), vectors(i))
       energies(i) = lowest_few(i)
    end do

    ! clean up 
    deallocate(diag, lowest_few, index)
  end subroutine init_from_diag

  !> \brief Runs a CI on just the reference space to build a few initial wavefunction guesses
  !> \param vectors  On exit contains the initial guess vectors (in sparse format)
  !> \param energies The energies of the initial vectors
  !> \param roots    The number of requested initial vectors
  !> \param sgga_data  Scratch space for the SGGA
  subroutine reference_ci(vectors, energies, roots, sgga_data)
    ! inputs
    type(sparse_vector), dimension(:) :: vectors
    real(kind=real8), dimension(:) :: energies
    integer, intent(in) :: roots
    type(sgga_struct) :: sgga_data

    ! local variables
    integer :: i, j, n, ierr
    integer :: nRefCSFs, addr_in_ref, addr_in_full
    integer, allocatable :: ref_addr_full_graph(:), ref_addr_ref_graph(:)  ! location of the reference configs on both graphs

    type(blockedLockVectorType) :: civec, sigmavec
    type(initial_dav_guess) :: ref_ci_guess
    procedure(htimesc), pointer :: mult
    procedure(filter_interface), pointer :: filter
    type(clock) :: clock_ref_ci

    
    write(*,*) " "
    write(*,*) "         ****************************************"
    write(*,*) " "
    write(*,*) "                    Reference CI Calculation                     "
    write(*,*) " "
    write(*,*) "         **************************************** "     
    call start_clock(clock_ref_ci)

   ! find where the references are in the wavefunction with just the references    
    call find_all_reference_CSFs(ref_addr_full_graph, nRefCSFs)

    ! Completely scrap the current SGGA graph and everything related to it
    ! rebuild it using the reference_ci_flag 
    write(*,*) 
    write(*,*) 
    write(*,*) "  Rebuilding the SGGA graph to include just the references "
    write(*,*) 
    write(*,*) 
    reference_ci_flag = 1
    call deallocLocScratch(sgga_data%loc_scr)
    call delete_graph_variables()
    call delete_allowed_virtuals()  
    call graph_driver(sgga_data%loc_scr)
    call delete_diagonal_elements()
    call diag_element(sgga_data%cho_data,sgga_data%loc_scr)
    rewind(iodiagstore)
    call write_big_vector(diagonal_elements,total_csfs,iodiagstore)
    deallocate(diagonal_elements, stat=ierr)
    
    ! setup for the davidson-liu algorithm 
    call allocateBlockedLockVector(civec, total_csfs, blockSizeCI)
    call allocateBlockedLockVector(sigmavec, total_csfs, blockSizeSigma)
    !mult => ref_ci_htimesc
    mult => htimesc
    call davidson_io_setup()

    ! handle the case where we have roots>num_ref
    n = min(num_ref, roots)

    ! get the initial guess for the reference ci 
    ! from the lowest <i|H|i> where i is a reference
    !filter => filter_everything_but_refs
    call davidson_setup(ref_ci_guess, n, .false., .false., sgga_data)

    ! run davidson-liu
    reference_ci_flag = 1
    call davidson_liu(civec, sigmavec, n, energies, ref_ci_guess, sgga_data, mult)
    reference_ci_flag = 0 
    call SGGA_reboot() 
    
    ! read in each final wavefunction and packed them into spare storage
    !
    ! this is probably the most dangerous part of this entire routine !!!!!!!!!
    ! I'm assuming that because we use the same algorithm for determining the 
    ! index for both the reference_ci and "normal" calculation that the references
    ! are in the same relative positions
    !
    ! For instance if the reference ci has references in positions: 1, 2, 3
    ! and the "normal" calculation has references in positions: 1, 19, 27
    ! I'm assuming 1->1 , 2->19, 3->27
  

    ! find where the references are in the wavefunction with just the references    
    call find_all_reference_CSFs(ref_addr_ref_graph, nRefCSFs)
    do i = 1, n
       call davidson_io_read(civec%v, ci_final, i)      
       do j = 1, nRefCSFs
          ! locations of this CSFs in both wavefunctions
          addr_in_ref = ref_addr_ref_graph(j)
          addr_in_full = ref_addr_full_graph(j)
          if ( abs(civec%v(addr_in_ref)) > 10E-10 ) then 
             call add_to_sparse_vector( civec%v(addr_in_ref) , addr_in_full, vectors(i))
          end if
       end do
    end do
    
    ! deallocate memory
    call deallocateBlockedLockVector(civec)
    call deallocateBlockedLockVector(sigmavec)
    call deallocate_initial_dav_guess(ref_ci_guess)   
    
    ! Scrap and rebuild the entire SGGA graph
    write(*,*) 
    write(*,*) 
    write(*,*) "  Rebuilding the SGGA graph to include all the desired configurations "
    write(*,*) 
    write(*,*) 
    rewind(unit=307)          ! don't ask me why, but if you forget this ... the next SGGA calc won't work
    reference_ci_flag = 0
    call deallocLocScratch(sgga_data%loc_scr)
    call delete_graph_variables()
    call delete_allowed_virtuals()
    call graph_driver(sgga_data%loc_scr)
    call delete_diagonal_elements()
    call diag_element(sgga_data%cho_data,sgga_data%loc_scr)
    rewind(iodiagstore)
    call write_big_vector(diagonal_elements,total_csfs,iodiagstore)
    deallocate(diagonal_elements, stat=ierr)
    
    ! if we need even more guesses ...
    if (roots > nRefCSFs) then 
       write(*,*) 
       write(*,*) "YOU ARE SOLVING FOR MORE ROOTS THAN REFERENCES"
       write(*,*) "INITIAL GUESSES /  REFERENCE ENERGIES GET AMBIGUOUS"
       write(*,*) "TREAT ALL RESULTS WITH CAUTION"
       write(*,*) 
       call handle_extra_ref_ci_guesses(vectors, energies, roots)
    end if
    
    write(*,*) 
    call print_clock(clock_ref_ci, "reference ci ")
  end subroutine reference_ci

!   !> \brief The matrix vector product H*c for a reference CI. Just a wrapper around htimesc
!   !> \param  x c 
!   !> \param  y Hc
!   !> \param data      SGGA scratch space
!   !>
!   !> We calculate the H*c product on the reference space by \f$ P^{-1}HPc\f$ where P is a projection
!   !> operator onto the reference space. This destroys everything in x
!   subroutine ref_ci_htimesc(x, y, data)
! 111    ! inputs
!     type(blockedLockVectorType) :: x, y
!     type(sgga_struct) :: data

!     ! local variables 
!     integer :: i, p 

!     ! first step, multiply P*x
!     y%v = 0.0
!     do i = 1, projector%n_ref_csfs
!        p = projector%ref(i)
!        y%v (p) = x%v (p)
!     end do
!     x%v = 0.0

!     ! now multiply H*(Px)
!     call htimesc(y, x, data)

!     ! now finish it with P^(-1) (HPx)
!     y%v  = 0.0
!     do i = 1, projector%n_ref_csfs
!        p = projector%ref(i)
!        y%v (p) = x%v (p)
!     end do
!   end subroutine ref_ci_htimesc

  !> \brief Generates extra guesses if I need more roots than references for an MR calculation 
  !> \param vectors  On exit contains the initial guess vectors (in sparse format)                                                   
  !> \param energies The energies of the initial vectors                                                                                              
  !> \param roots    The number of requested initial vectors                                                                                           
  !> \param sgga_data  Scratch space for the SGGA   
  !>
  !> We are going to use the lowest few diagonal elements (excluding the reference space)
  !> for the remaining guesses. I'm not sure if this is always the best idea -- be careful
  subroutine  handle_extra_ref_ci_guesses(vectors, energies, roots)
    ! inputs
    type(sparse_vector), dimension(:) :: vectors
    real(kind=real8), dimension(:) :: energies
    integer, intent(in) :: roots

    ! local variables
    integer :: i, ierr
    integer, dimension(:), allocatable :: index
    real(kind=real8), dimension(:), allocatable :: diag, lowest_few
    procedure(filter_interface), pointer :: filter

    ! read in the diagonal elements 
    allocate(diag(total_csfs), stat=ierr)
    call allocatecheck(ierr, "diag elements in init_from_diag")
    rewind(ioDiagStore)
    call read_big_vector(diag, total_csfs, ioDiagStore)

    ! build an index for sorting and a secondary array to store just the lowest few
    allocate(lowest_few(roots), index(roots), stat=ierr)
    call allocatecheck(ierr, "lowest_few and index")

    ! find the lowest few elements excluding the active space
    filter => filter_out_ref_space
    call find_lowest_n(diag, lowest_few, index, filter)

    ! set up the initial guesses
    do i = num_ref+1, roots
       call add_to_sparse_vector(1.0, index(i-num_ref), vectors(i))
       energies(i) = lowest_few(i-num_ref)
    end do

    ! clean up 
    deallocate(diag, lowest_few, index)
  end subroutine handle_extra_ref_ci_guesses

  !> \brief Find the lowest n elements of an array
  !> \param vector The vector being searched
  !> \param lowest The lowest n elements of the vector, in ascending order
  !> \param index  The position of the lowest elements in the vector
  !> \param filter_func  (Optional) A function pointer to a function to filter the elements of the array
  !>
  !> Note that the number of elements to find (n) is implicitly defined in the size of lowest and index
  !> which is a bit weird ... but n really doesn't appear (explicitly) anywhere in this subroutine
  !> 
  !> The filter function allows the user to remove some elements before sorting. For instance if you want 
  !> the lowest n elements of a vector, excluding a subspace, you can provide a function to remove those 
  !> unwanted elements before sorting.
  subroutine find_lowest_n(vector, lowest, index, filter_func)
    ! inputs
    integer, dimension(:) :: index
    real(kind=real8), dimension(:) :: vector, lowest
    procedure(filter_interface), pointer, optional :: filter_func

    ! local variables
    integer :: i,p
    integer :: length
    real(kind=real8) :: largest
    integer, dimension(1) :: position

    ! initialize the lowest and index arrays
    index = 0 
    lowest = maxval(vector) + 1

    ! if we need to filter, do it now
    if ( present(filter_func) ) call filter_func(vector)

    ! We search through the vector while maintaining a list
    ! of the lowest elements so far. For each element if
    ! it is smaller than the largest element on our lowest 
    ! list we replace the largest element with the new one.
    ! This means our lowest list, while accurate, isn't sorted.
    length = size(vector)
    do i = 1, length
       largest = maxval(lowest)
       position = maxloc(lowest)
       p = position(1)
       if ( vector(i) < largest) then
          lowest(p) = vector(i)
          index(p)  = i
       end if
    end do

    ! Now that we have a list of the n lowest elements, lets sort it into ascending order
    call  insertion_sort_realint( lowest, index )    
  end subroutine find_lowest_n

  !> \brief Defines an interface for the filter functions. Doesn't actually do anything
  subroutine filter_interface(x)
    real(kind=real8), dimension(:) :: x
    ! do nothing 
    x = x + 1
    write(*,*) "This is just a generic description of a filter function. You shouldn't be calling this."
    stop
  end subroutine filter_interface

  !> \brief Removes reference contributions from the vector (by setting them to a very large number)
  !> \param vector  The vector 
  subroutine filter_out_ref_space(vector)
    ! inputs
    real(kind=real8), dimension(:) :: vector

    ! local variables
    integer :: i, p, nRefCSFs
    integer, parameter :: large_number = 1000000000
    integer, allocatable :: references(:)
    
    call find_all_reference_CSFs(references, nRefCSFs)
    do i = 1, nRefCSFs
       p = references(i)
       vector(p) = 1000000000                                      
    end do
  end subroutine filter_out_ref_space

  ! !> \brief Removes everything except the reference space contributions from the vector (by setting them to a very large number)
  ! !> \param vector The vector
  ! subroutine filter_everything_but_refs(vector)
  !   ! inputs
  !   real(kind=real8), dimension(:) :: vector

  !   ! local variables
  !   integer :: i, p, ierr, nRefCSFs
  !   integer, allocatable :: references(:)
  !   integer, parameter :: large_number = 1000000000
  !   real(kind=real8), dimension(:), allocatable :: tmp

  !   call find_all_reference_CSFs(references, nRefCSFs)

  !   ! allocate some temp space and place the reference values here
  !   allocate(tmp(nRefCSFs), stat=ierr)
  !   call allocatecheck(ierr, "tmp space in filter_everything_but_refs")
  !   do i = 1, nRefCSFs
  !      p = references(i)
  !      tmp(i) = vector(p)
  !   end do

  !   ! set the entire vector
  !   vector = large_number

  !   ! restore the old values in the reference space
  !   do i = 1, nRefCSFs
  !      p = references(i)
  !      vector(p) = tmp(i)
  !   end do
  ! end subroutine filter_everything_but_refs


end module initial_davidson_guess








  
