! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! ! Out of core cholesky decomposition of the AO two electron integrals
! !
! ! note I'm abbreviating out-of-core as ooc in some function names/variables
! !
! ! note that when I store columns of the matrix on disk the format is column(:), i
! ! where i is an integer denoting the last transposition(pivot) applied to that column



module cholesky_out_of_core
  use global_var_mod
  use molecule_var_mod
  use utilities_mod
  use math_utils
  use cholesky_structs
  use integral_screen
  use fortran_timing
  use get_integrals
  use function_map
  use c_io_finterface
  !debug
  use sort_utils
  implicit none

  type transposition
     integer :: a,b
  end type transposition

  !> When storing a portion of a matrix in a buffer the logic for figuring what column of the matrix
  !> corresponds to which column of the buffer can get ugly (and buggy). This simplifies it a bit by
  !> having the buffer keep track for us
  type indexed_buffer
      real(kind=real8), dimension(:,:), allocatable :: buff
      integer, dimension(:), allocatable :: actual_index
  end type indexed_buffer
  
contains


  subroutine ooc_calculate_two_electron_integrals(nShell, pairMap, prescreened, nFuncInShell, screen, disk, &
       activePairs, disk_pointers )
    implicit none

    ! inputs
    integer, intent(in) ::  nFuncInShell(:)
    type(pairIndex) :: pairMap
    logical, intent(in) :: prescreened(:,:) 
    integer :: nShell, disk, activePairs
    integer, dimension(:) :: disk_pointers
    type(integral_screen_data) :: screen                    ! for integral screening

    ! local variables
    integer :: iShell, jShell, kShell, lShell
    integer :: iFunc, jFunc, kFunc, lFunc
    integer :: ij, kl, i, j, k, l
    integer :: ip, jp, kp, lp
    integer :: ijkl, ijShell, klShell
    integer :: ierr, iBuff
    integer :: intSize
    integer :: my_iShell, my_jShell, my_kShell, my_lShell, tmp
    integer :: am_i, am_j, am_k, am_l
    integer :: p_i, p_ij, p_ijk

    integer, dimension(:), allocatable :: buffer_map
    real(kind=real8), dimension(:), allocatable :: scratch
    real(kind=real8), dimension(:,:), allocatable :: buffer

    type(clock) :: timer                                    ! timing variable

    integer :: rec_size, info

    ! Start the timer
    call start_clock(timer)

    ! Memory allocations (we need scratch space for the calculations)
    intSize = maxval(nFuncInShell) ** 2 
    allocate(buffer(activePairs,intSize), buffer_map( sum(nFuncInShell) ** 2 ), stat=ierr)
    call allocatecheck(ierr, "more scratch space for (ij|kl")

    intSize = maxval(nFuncInShell) ** 4
    allocate( scratch(intSize), stat=ierr)
    call allocatecheck(ierr, "scratch space for (ij|kl)")

    ! Lets open the file we are going to write the integrals to
    inquire(iolength=rec_size) buffer(:,1), info
    open(unit=disk,file=scratch_directory // 'two_electron_integrals.dat',form='unformatted',access='direct',recl=rec_size)

    ! Alright lets do the hard work. Note that this, without (pre)screening
    ! scales as O(nBas**4)

    ! Integral evaluations 
    !***********************      
    ! loop over shells for the eval_ijkl call at each iteration 
    do iShell = 1, nShell
       do jShell= 1, iShell

          ! reorder shells for proper erkale angular momentum ordering
          my_iShell = iShell
          my_jShell = jShell
          
          call funcs_to_am(nFuncInShell(my_iShell), am_i)
          call funcs_to_am(nFuncInShell(my_jShell), am_j)
         
          if (am_i < am_j) then
             tmp = my_iShell 
             my_iShell = my_jShell
             my_jShell = tmp
          end if

          call set_buffer_map(buffer_map, pairMap, prescreened, nFuncInShell, my_iShell, my_jShell)
          buffer = 0
          ! Calculate the columns (**|IJ) ... note that I need the whole thing, not just the
          ! lower triangular part

          do kShell = 1, nShell
             do lShell= 1, kShell
                ijShell = iShell*(iShell+1)/2 + jShell
                klShell = kShell*(kShell+1)/2 + lShell

                ! Integral screening
                if ( shellIsScreened(screen, iShell, jShell, kShell, lShell, nFuncInShell) ) cycle

                ! second half of shell reordering
                my_kShell = kShell
                my_lShell = lShell

                call funcs_to_am(nFuncInShell(my_kShell), am_k)
                call funcs_to_am(nFuncInShell(my_lShell), am_l)
                
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
                         p_ijk = (p_ij + kFunc -1) * nFuncInShell(my_lShell)
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

                            iBuff = buffer_map(ij)
                            buffer(kl,iBuff) = scratch(ijkl)

                         end do
                      end do
                   end do
                end do


             end do
          end do

          ! We finished with the (**|IJ) columns
          ! write each column out to disk 
          call write_integral_buffer(buffer, buffer_map, pairMap, prescreened, nFuncInShell, iShell, jShell, disk, disk_pointers)

       end do
    end do

    ! clean up 
    deallocate(scratch)
    ! produce a bit of output
    write(*,*) " "
    write(*,*) "=== Calculated the necessary (ij|kl) ==="
    call print_integral_screen_stats(screen)
    call print_clock(timer, "(ij|kl) evaluations")

  end subroutine ooc_calculate_two_electron_integrals



  !> \brief Writes out a set of columns (**|IJ) to disk 
  !> \param buffer     Contains the columns (**|IJ)
  !> \param buffer_map The mapping between the buffer and the (ij|kl)
  !> \param pairMap    The data structure mapping the (ij) pairs to i,j
  !> \param prescreened
  !> \param nFuncInShell The number of AO basis functions in each AO Shell
  !> \param iShell     The I from the set of columns (IJ) 
  !> \param jShell     The J from the set of columns (IJ) 
  !> \param disk       Unit number to store the columns
  !> \param disk_pointers Maps column number to record number 
  subroutine write_integral_buffer(buffer, buffer_map, pairMap, prescreened, nFuncInShell, iShell, jShell, disk, disk_pointers)
    implicit none

    ! inputs
    integer, intent(in) :: iShell, jShell, nFuncInShell(:), disk, buffer_map(:)
    integer, intent(out) :: disk_pointers(:)
    type(pairIndex) :: pairMap   
    real(kind=real8), dimension(:,:) :: buffer
    logical , intent(in) :: prescreened(:,:)

    ! local variables 
    integer :: i,j
    integer :: iFunc, jFunc
    integer :: ip,jp
    integer :: ij

    do iFunc = 1, nFuncInShell(iShell)
       ip = sum(nFuncInShell(1:iShell-1))
       do jFunc = 1, nFuncInShell(jShell)
          jp = sum(nFuncInShell(1:jShell-1))
          i = iFunc + ip
          j = jFunc + jp 
          ij = pairMap%ij2pair(i,j)      
          ! skip this pair if it isn't on the map .. it was either prescreened or we already finished this pair
          if (prescreened(i,j)) cycle
          if ( buffer_map(ij) == 0  ) cycle
          disk_pointers(ij) = ij
          write(unit=disk, rec=ij) buffer(:,buffer_map(ij)),0
       end do
    end do

  end subroutine write_integral_buffer



  !> \brief Determines the mapping between columns stored in the buffer and the columns of the two-electron integral matrix
  !> \param buffer_map The mapping between the buffer and the (ij|kl)
  !> \param pairMap    The data structure mapping the (ij) pairs to i,j
  !> \param prescreened
  !> \param nFuncInShell The number of AO basis functions in each AO Shell
  !> \param iShell     The I from the set of columns (IJ) we are going to store in the buffer
  !> \param jShell     The J from the set of columns (IJ) we are going to store in the buffer
  subroutine set_buffer_map(buffer_map, pairMap, prescreened, nFuncInShell, iShell, jShell)
    implicit none

    ! inputs 
    integer, intent(in) :: iShell, jShell, nFuncInShell(:)
    integer, intent(out) :: buffer_map(:)
    type(pairIndex) :: pairMap
    logical , intent(in) :: prescreened(:,:)

    ! local variables
    integer :: i,j
    integer :: iFunc, jFunc
    integer :: ip,jp
    integer :: ij
    integer :: iBuff

    buffer_map = 0
    iBuff = 1

    do iFunc = 1, nFuncInShell(iShell)
       ip = sum(nFuncInShell(1:iShell-1))
       do jFunc = 1, nFuncInShell(jShell)
          jp = sum(nFuncInShell(1:jShell-1))
          i = iFunc + ip
          j = jFunc + jp 
          ij = pairMap%ij2pair(i,j)      

          ! skip if we removed this pair or if it equivalent to one we already did ( (ij) and (ji))
          ! note that the order of these statements is important !!!! (if (ij) is prescreened then ij has a weird value
          if (prescreened(i,j)) cycle
          if (buffer_map(ij) /= 0 ) cycle
          buffer_map(ij) = iBuff 
          iBuff = iBuff + 1      


       end do
    end do

  end subroutine set_buffer_map

  !> \brief This routine takes a vector already pivoted to iteration a and updates it with pivots up to (and including) iteration b
  !> \param vec The vector needing pivoting
  !> \param a   The last element pivoted in vec
  !> \param b   The extent of new pivoting
  !> \param pivot_steps  A list of transpositions that record the pivot steps
  subroutine pivot_update( vec, a, b, pivot_steps)
    implicit none

    ! inputs
    real(kind=real8), dimension(:) :: vec
    integer, intent(in) :: a,b
    type(transposition), dimension(:), intent(in) :: pivot_steps

    ! local variables
    integer :: i, start, finish

    start = a + 1
    finish = b
    do i = start, finish
       call swap_element_real(vec, pivot_steps(i)%a, pivot_steps(i)%b)
    end do
  end subroutine pivot_update

  !> \brief Swaps two elements of a real vector
  !> \param vec The vector whose elements I'm swapping
  !> \param j   Index of one element to swap
  !> \param q   Index of the other element to swap
  subroutine swap_element_real(vec,j,q)
    implicit none
    real(kind=real8)::vec(:),tmp
    integer::j,q
    tmp = vec(q)
    vec(q) = vec(j)
    vec(j) = tmp
  end subroutine swap_element_real


  !> \brief Swaps two elements of an integer vector
  !> \param vec The vector whose elements I'm swapping
  !> \param j   Index of one element to swap
  !> \param q   Index of the other element to swap
  subroutine swap_element_int(vec,j,q)
    implicit none
    integer::vec(:),tmp
    integer::j,q
    tmp = vec(q)
    vec(q) = vec(j)
    vec(j) = tmp
  end subroutine swap_element_int


  !> \brief Swaps two rows of a matrix
  !> \param A The matrix whose rows I'm swapping
  !> \param j   Index of the first row to swap
  !> \param q   Index of the other row to swap
  !> \param n The size of the matrix A(nxn)
  subroutine swap_row(A,j,q,n)
    implicit none
    real::A(:,:), tmp
    integer::i,j,q,n
    do i = 1, n
       tmp = A(q,i)
       A(q,i) = A(j,i)
       A(j,i) = tmp
    end do
  end subroutine swap_row

  !> \brief Swaps two columns of a matrix
  !> \param A The matrix whose columns I'm swapping
  !> \param j   Index of the first column to swap
  !> \param q   Index of the other column to swap
  !> \param n The size of the matrix A(nxn)
  subroutine swap_col(A,j,q,n)
    implicit none
    real(kind=real8)::A(:,:), tmp
    integer::i,j,q,n
    do i = 1, n
       tmp = A(i,q)
       A(i,q) = A(i,j)
       A(i,j) = tmp
    end do
  end subroutine swap_col

  !> \brief Finds the next pivot for the cholesky decomposition by finding the largest remaining diagonal element
  !> \param diag The current values (updated in the outer product formulation) of the diagonal elements
  !> \param j The current iteration
  !> \param q On exit the column to pivot to
  !> \param n The size of the matrix A(nxn)
  subroutine find_pivot(diag, j, q, n)
    implicit none
    real(kind=real8) :: diag(:), max
    integer :: j,q,n, max_pos,i
    max_pos=j
    max=diag(j)
    do i = j, n
       if (diag(i) > max) then
          max = diag(i)
          max_pos=i
       end if
    end do
    q = max_pos
  end subroutine find_pivot



  !> \brief Performs the out-of-core (ooc) cholesky. 
  !> \param C The number of columns of the matrix that we can store in memory at one time 
  !> \param diag The diagonal elements of the matrix to be decomposed
  !> \param n The size of the matrix (nxn) to be decomposed
  !> \param tol The numerical tolerance (cut-off point) for the incomplete cholesky
  !> \param P   On exit contains the permutation used during pivoting. P(i) = j indicates that the ith column is now the jth column
  !> \param disk The unit number where the vectors are stored. File should already be open
  !> \param disk_pointers Maps column number to record number 
  !> \param rank The effective rank of the matrix being decomposed (on exit)
  !> \param activePairs     The number of basis products (ij) that are included in the decomposition after prescreening
  !>
  !> I'm making no claims on the efficency of this routine, if you can think of something better than go for it :)
  !> This routine assumes you have already written the columns of the matrix to be decomposed to disk. 
  subroutine ooc_cholesky(C, diag, n, tol, P, disk, disk_pointers, rank, activePairs)
    implicit none
    ! inputs
    integer, intent(in) :: n, activePairs, disk, C, disk_pointers(:)
    real(kind=real8), intent(in) :: tol
    real(kind=real8), dimension(:) :: diag
    integer, dimension(:) :: P
    integer, intent(out) :: rank

    ! local variables
    integer :: ierr
    integer :: i,j,k,q
    integer :: piv_stop
    integer :: p_buff
    real(kind=real8), allocatable :: buffer2(:)

    type(transposition), allocatable :: pivots(:)
    type(indexed_buffer) :: buff

    ! Initial Setup
    ! **************
    allocate(buff%buff(n,C), buff%actual_index(C), buffer2(n),pivots(n), stat = ierr)
    call allocatecheck(ierr, "initialization in ooc_cholesky")
    do i = 1, n 
       P(i) = i 
       pivots(i)%a = 0
       pivots(i)%b = 0
    end do
    
    
    ! Decomposition
    ! **************
    p_buff = 1
    buff%buff = 0  
    buff%actual_index = -1
    do j = 1 , n

       ! find the next pivot
       call find_pivot(diag,j,q,n)

       ! check to see if we are done
       if ( diag(q) < tol) then
          ! we finished the decomposition 
          rank = j - 1
          exit
       end if

       ! start pivoting
       call swap_element_real(diag,j,q)
       call swap_element_int(P,j,q)
       call swap_element_int(disk_pointers,j,q)
       pivots(j)%a = j
       pivots(j)%b = q
       call swap_row(buff%buff(:,1:C),j,q,C)

       ! read in the new column
       call load_into_index_buff(buff, j, p_buff, n, disk , disk_pointers, piv_stop)
       call pivot_update(buff%buff(:,p_buff), piv_stop, j, pivots)

       ! now we decompose using the inner product formulation
       buff%buff(j,p_buff) = sqrt( diag(j) )
       call dgemv('N',(n-j),(p_buff-1),-1.0,buff%buff(j+1,1),n,buff%buff(j,1:p_buff-1),1,1.0,buff%buff(j+1:n,p_buff),1)
       buff%buff(j+1:n,p_buff) = buff%buff(j+1:n,p_buff) / buff%buff(j,p_buff)

       ! update the diagonal elements
       do i = j+1, n
          diag(i) = diag(i) - buff%buff(i,p_buff) ** 2
       end do

       ! move to the next column in the buffer
       p_buff = p_buff + 1 

       ! check to see if the buffer is full. if it is, perform the outer product update
       ! using the decomposed columns in the buffer
       if ( p_buff == C + 1 ) then 
          ! we are going to "retire" everything currently in the buffer
          ! this means we are going to add their contributions to everything on disk
          ! and we won't need these vectors again 
          do k = j + 1, n 

             ! read in a vector I haven't decomposed yet
             read(unit=disk, rec=disk_pointers(k)) buffer2(:), piv_stop
             call pivot_update(buffer2, piv_stop, j, pivots)          

             ! update the vector with the contributions from everything in buff%buff
             ! Note that while this approach originates with the outer-product formulation
             ! of the cholesky decomposition, that doing it this way (with each not-decomposed
             ! column on disk) ends up organizing the update in matrix-vector product instead of 
             ! an outer product
             call dgemv('N',(n-j+1),C,-1.0,buff%buff(j,1),n,buff%buff(k,1),n,1.0,buffer2(j:n),1)

             ! write the new vector back to disk 
             write(unit=disk, rec=disk_pointers(k)) buffer2(:), j

          end do

          ! dump the vectors in buff to disk 
          call dump_buffer(buff,disk,disk_pointers,j)

          ! reset buff and p_buff
          p_buff = 1
          buff%buff = 0 
          buff%actual_index = -1
       end if

    end do ! cholesky decomposition

    ! dump whatever vectors remain in the buff to disk 
    call dump_buffer(buff,disk,disk_pointers,n)


    ! run through a final set of pivot updates for all vectors on disk
    ! then place the final vectors in the correct file (cvec_no) for the transform code
    ! this also unscrambles the vectors on disk ( column i is now in record number i )
    call setup_vector_io_with_c(scratch_directory // 'cvec.dat')
    do j = 1 , n 
       read(unit=disk, rec=disk_pointers(j)) buffer2(:), piv_stop
       call pivot_update(buffer2, piv_stop, rank, pivots)
       call write_vector_with_c(scratch_directory // 'cvec.dat', buffer2(:), activePairs,j)
    end do

    ! Produce a bit of output about what happened
    write(*,*)  " "
    write(*,*) "=== cholesky factorization ==="
    write(*,*) "CD threshold = " , tol
    write(*,*) "Finished the OOC-CD "
    write(*,*) "effective rank of the two-electron integrals = ", rank

  end subroutine ooc_cholesky


  !> \brief Writes out the permutation from our CD to disk
  !> \param P Contains the permutation matrix used for the pivoted decomposition
  !> \param pairMap The data structure mapping the (ij) pairs to i,j
  !> \param activePairs     The number of basis products (ij) that are included in the decomposition after prescreening
  subroutine ooc_write_index(P, pairMap, activePairs)

    implicit none

    ! inputs
    integer, intent(in) :: activePairs, P(:)
    type(pairIndex), intent(in) :: pairMap

    ! local variables
    integer :: i
    integer :: iPair
    integer, allocatable :: pair(:)

    ! Write out the permuted pair (ij) basis
    allocate(pair(2))
    call setup_vector_io_with_c(scratch_directory // 'cvec_index.dat')
    do i = 1, activePairs
       iPair = P(i)
       pair(1)  = pairMap%pair2ij(iPair)%i
       pair(2)  = pairMap%pair2ij(iPair)%j
       call write_ints_with_c(scratch_directory // 'cvec_index.dat', pair, 2, i)
    end do
    deallocate(pair)
  end subroutine ooc_write_index

  !> \brief Figures out how many colummns of the two-electron integral matrix we can store at one time
  !> \param n Size of the two electron integral matrix (nxn)
  !> \param mem Amount of memory available in # of reals
  !> 
  integer function ooc_mem(n, mem)
    ! implicit none
    integer, intent(in) :: n, mem
    integer :: C

    if ( n*n < mem ) then
       ! We can store them all
       C = n*n + 1        
    else
       C = int(mem/n)
       ! remove one column for the second buffer
       ! remove one column for the diagonal integrals
       ! remove one column for the pairMap
       ! remove one column for the integral screening
       ! remove one column for the disk pointers
       C = C - 5
    end if
    
    if ( C < 0 ) then
        write(*,*) "Bad news .... you don't have enough memory to do the ooc CD ... sorry"
        write(*,*) "FATAL ERROR"
        stop
    end if
    ooc_mem = C
  end function ooc_mem

  
  !> \brief This routine takes the diagonal integrals and produces an array of just the (ij|ij) integrals for active (not prescreened) (ij)
  !> \param all_diag   All the diagonal integrals
  !> \param activePair_diag  On exit contains the diagonal integrals for just the active (ij) pairs
  !> \param nBas         The number of atomic orbitals 
  !> \param prescreened  A logical array describing which pairs are prescreened
  !> \param pairMap      A mapping (of the not prescreened) (ij) pairs
  !> 
  !> This routine is nice because (a) it reduces the memory footprint and (b) makes life less complicated when doing the ooc ( I don't
  !>  have the headache of keeping track of which parts of the diag array aren't elements in the matrix I'm decomposing)
  subroutine get_active_diag( all_diag, activePair_diag, nBas, prescreened, pairMap)
      implicit none
      
      ! inputs
      real(kind=real8), intent(in) :: all_diag(:)
      real(kind=real8), intent(out) :: activePair_diag(:)
      integer, intent(in) :: nBas
      logical, intent(in) :: prescreened(:,:)
      type(pairIndex)   :: pairMap
      
      ! local variables
      integer :: iBas, jBas
      integer :: iDiag
      integer :: iPair
      
      iDiag = 0
      do iBas = 1, nBas
          do jBas = 1, iBas
              iDiag = diag_address(iBas, jBas)
              if ( prescreened(iBas,jBas) ) cycle
              iPair = pairMap%ij2Pair(iBas,jBas)
              activePair_diag(iPair) = all_diag(iDiag)
          end do
      end do
  end subroutine get_active_diag
  
  !> \brief Loads a vector from disk into an indexed buffer
  !> \param buff The buffer we load into
  !> \param j    The column of the matrix we load into the buffer
  !> \param p_buff The position in the buffer where column j is stored
  !> \param n    The size of the column
  !> \param disk  The unit number we are reading the column from 
  !> \param mapping  A mapping of column index (j) to record number on disk
  !> \param piv_stop  The last element of the vector to be pivoted
  subroutine load_into_index_buff(buff, j, p_buff, n, disk , mapping, piv_stop)
      implicit none
      
      ! inputs
      integer, intent(in) :: n, disk, j, p_buff, mapping(:)
      integer, intent(out) :: piv_stop
      type(indexed_buffer) :: buff
       
      read(unit=disk, rec=mapping(j)) buff%buff(1:n,p_buff), piv_stop
      buff%actual_index(p_buff) = j 
  end subroutine load_into_index_buff


  !> \brief Writes out columns from an indexed buffer to disk
  !> \param buff   The indexed buffer
  !> \param disk   The unit number we should write the columns to
  !> \param mapping  Maps column number to record number 
  !> \param piv     The last pivot step applied to vectors [a,b] (this information is also written to disk)
  subroutine dump_buffer( buff, disk, mapping, piv)
    implicit none
    ! inputs
    integer, intent(in) :: mapping(:), disk, piv
    type(indexed_buffer) :: buff

    ! local variables
    integer :: i
    integer :: buf_length
    integer :: record, col
    
    buf_length = size(buff%actual_index)
    do i = 1, buf_length
        col = buff%actual_index(i)             ! what column of the matrix this column of the buffer corresponds to
        if ( col > 0 ) then
            record = mapping(col)              ! where this column of the whole matrix is located on disk
            write(unit=disk, rec=record) buff%buff(:,i), piv
        end if
    end do
  end subroutine dump_buffer
  
end module cholesky_out_of_core
