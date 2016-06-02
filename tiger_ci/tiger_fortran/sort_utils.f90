! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief Contains a few basic sorting algorithms
!> \author David Krisiloff
!> \date 12/23/2010
!>
!> Sorting utilities heavily based on discussions in
!> Numerical Recipies and Introduction to Algorithms
!> by Cormen. (And also some notes from an Algorithms
!> and data structure class at princeton)
!>
!> I really wish Fortran has some basic sorting utilities
!> as part of the language ... sigh


module sort_utils
  use random
  use global_var_mod
  implicit none
  private 
  public insertion_sort_real, insertion_sort_int, shellsort_real, hybrid_quicksort_real, &
       iquicksrt, insertion_sort_realint

contains

  !> \brief Insertion sort for an array of real(kind=real8)s
  !>
  !> Note that insertion sort is fairly
  !> fast for a only a few elements. Its \f$O(N^2)\f$
  !> though so avoid using it for anything more than say
  !> 10 elements.
  !> \param array An array of real(kind=real8)s to be sorted
  subroutine insertion_sort_real( array )
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: array
    integer :: i,j
    do i = 2 , size(array)             ! loop over the array
       do j = i , 2 , -1              ! start with j=i and decrease by 1 until j=1
          if ( array(j) < array(j-1) ) then
             call swap_real(array,j,j-1)
          else
             exit
          end if
       end do
    end do
  end subroutine insertion_sort_real

  !> \brief Insertion sort for an array of ints
  !>
  !> Note that insertion sort is fairly
  !> fast for a only a few elements. Its \f$O(N^2)\f$
  !> though so avoid using it for anything more than say
  !> 10 elements.
  !> \param array An array of real(kind=real8)s to be sorted
  subroutine insertion_sort_int( array )
    implicit none
    integer, dimension(:), intent(inout) :: array
    integer :: i,j
    do i = 2 , size(array)             ! loop over the array
       do j = i , 2 , -1              ! start with j=i and decrease by 1 until j=1
          if ( array(j) < array(j-1) ) then
             call swap_int(array,j,j-1)
          else
             exit
          end if
       end do
    end do
  end subroutine insertion_sort_int


  !> \brief Shellshort an array of real(kind=real8)s
  !> Implements a shellsort using the increment sequence recommend by Knuth
  !> \param array An array of real(kind=real8)s to be sorted
  subroutine shellsort_real( array )
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: array
    integer                        :: length
    integer                        :: h,i,j

    length = size(array)
    ! First lets build the largest value of h we are going to use
    h = 1
    do while ( h < length/3 )
       h = 3*h +1
    end do

    do while ( h > 0 )
       ! h sort the array
       ! essentially an insertion sort of every hth element
       do i = h  , length
          j = i
          do while ( j > h .and. array(j) < array(j-h) )
             call swap_real(array,j,j-h)
             j = j - h
          end do
       end do
       h = h / 3       ! get the next value in the sequence
    end do

  end subroutine shellsort_real

  !> \brief An fancy recursive version of quicksort
  !>
  !>  Quicksorts into groups of 10 elements then finishes it
  !>  off with insertion sort. Uses a random shuffle to minimize
  !>  chances of finding a pathological case. 
  !> \param array An array of real(kind=real8)s to be sorted
  !> \param low   Lower bound of the array (inclusive)
  !> \param high  Upper bound of the array (inclusive)
  subroutine hybrid_quicksort_real(array, low, high)
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: array
    integer, intent(in)              :: low, high

    ! start up the random number generator
    call setup_rand()

    ! the hybrid method uses quicksort to sort everything into
    ! sort subgroups of 10 elements. In side the subgroup the elements
    ! are not sorted. it also uses a random shuffle of the elements
    ! to avoid the worse case scenario of the partitioning element being
    ! the highest or lowest element
    call quicksort_real(array, low, high)

    ! now we need to sort in each of the groups
    ! do this quickly with insertion sort
    call insertion_sort_real(array)


  end subroutine hybrid_quicksort_real


  !> \brief An extremely simple recursive version of quicksort
  !>
  !> The very basics ... no fancy bells or whistles yet
  !> \param array An array of real(kind=real8)s to be sorted
  !> \param low   Lower bound of the array (inclusive)
  !> \param high  Upper bound of the array (inclusive)
  recursive subroutine quicksort_real( array, low, high )
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: array
    integer, intent(in)              :: low, high
    integer :: i
    ! START AUTOGENERATED INITIALIZATION 
    i = 0
    ! END AUTOGENERATED INITIALIZATION 

    ! end of recursion
    ! if we only have ten elements left than ignore it
    ! and we will finish it off with insertion sort later !
    if ( (high-low) <= 10 ) then
       return
    end if

    ! Quicksort works by setting a partitioning element i
    ! and placing it in its correct position in the array
    ! all elements below i are < i and all elements above i
    ! are > i. This is called partitioning
    !
    ! This step is then repeated recurisvely for each subarray,
    ! the elements below i and the elements above i

    call partition(array,low,high,i)
    call quicksort_real(array,low,i-1)
    call quicksort_real(array,i+1,high)

  end subroutine quicksort_real

  subroutine partition(array, low, high, q)
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: array
    integer, intent(in)              :: low, high
    integer, intent(out)             :: q

    real(kind=real8) :: val                   ! the value of the partitioning element
    integer :: left, right

    ! Ok we need to decide on the partitioning element
    ! arbitrarily we choose the last element array(high) 
    ! after the array has been randomly shuffled
    ! left and right are pointers to the left and right side of the array
    call shuffle_real( array(low:high) )
    val = array(high)
    left = low - 1
    right = high

    !        write(*,*) "Starting up the paritioning ... partition value is = " , val , 'bounds are = ' , low, high
    !        write(*,*) "array = " , array
    !        call flush(6)

    ! Now we partition the array by finding one element on the left
    ! which is greater than the partition and one element on the right which is
    ! less than the partition
    do 
       ! stop condition
       if ( left >= right ) then
          exit
       end if

       ! search the left side for a large element
       ! starting one after where we stopped last time
       do left = left + 1 , right - 1 
          if ( array(left) > val ) then
             exit
          end if
       end do

       ! search the right side for a small element
       ! notice we loop in a decreasing manner
       do right = right -1  , left + 1, -1
          if ( array(right) <= val ) then
             exit
          end if
       end do

       !            write(*,*) "Pointers are now "
       !            write(*,*) "left = ", left, "val = " , array(left)
       !            write(*,*) "Right = " , right, "val = " , array(right)

       ! ok now we have found our two element which we can switch
       if ( left < right ) then
          call swap_real(array,left,right)
          !                write(*,*) "swapped array = " , array
       end if
    end do

    ! now we need to swap the partitioning element into place
    ! its currently at left.
    call swap_real(array, left, high )
    q = left

    !        write(*,*) "swapping in the partitioning element "
    !        write(*,*) "array = " , array
    !        write(*,*) " "
    !        write(*,*) " "

  end subroutine partition


  !> \brief Insertion sort for an array of real(kind=real8)s
  !>
  !> Note that insertion sort is fairly
  !> fast for a only a few elements. Its \f$O(N^2)\f$
  !> though so avoid using it for anything more than say
  !> 10 elements.
  !> \param real(kind=real8)_array An array of real(kind=real8)s to be sorted
  !> \param int_array  An array of ints used to index real(kind=real8)_array
  subroutine insertion_sort_realint( real_array, int_array )
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: real_array
    integer, dimension(:), intent(inout) :: int_array
    integer :: i,j

    i = 0
    j = 0

    do i = 2 , size(real_array)             ! loop over the array
       do j = i , 2 , -1              ! start with j=i and decrease by 1 until j=1
          if ( real_array(j) < real_array(j-1) ) then
             call swap2_realint(real_array, int_array,j,j-1)
          else
             exit
          end if
       end do
    end do

  end subroutine insertion_sort_realint


  !> \brief A little helper function swaps elements i,j in array
  !> \param array An array of real(kind=real8)s
  !> \param i     the position of the first element to be swapped
  !> \param j     the position of the second element to be swapped
  subroutine swap_real( array, i, j)
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: array
    integer , intent(in) :: i,j
    real(kind=real8) :: tmp
    tmp = array(i)
    array(i) = array(j) 
    array(j) = tmp
  end subroutine swap_real


  !> \brief A little helper function swaps elements i,j in array
  !> \param array An array of ints
  !> \param i     the position of the first element to be swapped
  !> \param j     the position of the second element to be swapped
  subroutine swap_int( array, i, j)
    implicit none
    integer, dimension(:), intent(inout) :: array
    integer , intent(in) :: i,j
    integer:: tmp
    tmp = array(i)
    array(i) = array(j) 
    array(j) = tmp
  end subroutine swap_int


  !> \brief A little helper function swaps elements i,j in array and intarray
  !> \param array An array of real(kind=real8)s
  !> \param intarray An array of ints
  !> \param i     the position of the first element to be swapped
  !> \param j     the position of the second element to be swapped
  !>  
  !>   This is real(kind=real8)ly a private method of this module ...
  !>   you probably don't want to be calling this from anywhere else
  !>   (though I suppose technically you could).
  subroutine swap2_realint( array, intarray, i, j)
    implicit none
    real(kind=real8), dimension(:), intent(inout) :: array
    integer,  dimension(:), intent(inout) :: intarray
    integer , intent(in) :: i,j
    integer :: inttmp
    real(kind=real8) :: tmp

    tmp = array(i)
    array(i) = array(j)
    array(j) = tmp

    inttmp = intarray(i)
    intarray(i) = intarray(j)
    intarray(j) = inttmp        
  end subroutine swap2_realint



  !> \brief Not sure where this integer quicksort comes from. Try to avoid using it. Keep here only to retain compatability 
  !> with the CD code
  !> \todo REPLACE THIS CODE
  SubRoutine Iquicksrt(iVec,iPos,nVec,iStack)
    ! 
    !     non-recursive quicksort for an integer vector ...
    !     sorts vector in ascending order ...
    !     iVec(nVec)     : vector to sort, sorted on return
    !     iStack(2,nVec) : scratch space
    !     Hacked to also return position indices in iPos ...

    Implicit None

    !     declare subroutine paramemetrs...
    Integer nVec
    Integer iVec(nVec),iStack(2,nVec)
    Integer iPos(nVec)
    !     declare some local variables...
    Integer iLeft,iRight,iS,iX,iW
    Integer i,j
    Integer idum
    !
    iS = 1
    iStack(1,1)=1
    iStack(2,1)=nVec
    !     REPEAT...UNTIL
100 Continue
    iLeft = iStack(1,iS)
    iRight = iStack(2,iS)
    iS = iS-1
    !       REPEAT...UNTIL
110 Continue
    i = iLeft
    j = iRight
    iX = iVec((iLeft+iRight)/2)
    !         REPEAT...UNTIL
120 Continue
    !           WHILE...DO
200 If (iVec(i).ge.iX) Go To 201
    i = i+1
    Go To 200
201 Continue
    !           END WHILE...DO
    !           WHILE...DO
300 If (iX.ge.iVec(j)) Go To 301
    j = j-1
    Go To 300
301 Continue
    If (i.le.j) Then
       !             swap...
       iW = iVec(i)
       idum = iPos(i)
       iVec(i) = iVec(j)
       iPos(i) = iPos(j)
       iVec(j) = iW
       iPos(j) = idum
       i = i+1
       j = j-1
    End If
    If (i.le.j) Go To 120
    !         END REPEAT...UNTIL
    If (i.lt.iRight) Then
       iS = iS+1
       iStack(1,iS) = i
       iStack(2,iS) = iRight
    End If
    iRight = j
    If (iLeft.lt.iRight) Go To 110
    !       END REPEAT...UNTIL
    If (iS.ne.0) Go To 100
    !     END REPEAT...UNTIL
    Return

  end subroutine Iquicksrt

end module sort_utils
