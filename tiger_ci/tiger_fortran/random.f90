! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief This module contains utility functions involving pseudo-random numbers

module random
    use global_var_mod
    private
    public :: shuffle_real, shuffle_int, setup_rand
    contains
    
  !> \brief Randomly shuffles an array of real(kind=real8)s
  !>
  !>  Needed for a quicksort algorithm which doesn't hit \f$O(N^2)\f$ on a sorted array.
  !>  Uses the Fisher-Yates algorithm.
  !>  \param array An array of real(kind=real8)s to be shuffled randomly
  subroutine shuffle_real ( array )
    implicit none

    real(kind=real8), dimension(:), intent(inout) :: array
    integer :: i,j,n
    real(kind=real8) :: random, tmp
    
    ! Note we are assuming that the random number generator has been seeded correctly!
    ! Now lets shuffle by randomly picking a number on [i,n] and placing it at position array(i)
    n = size(array)
    do i = 1 , n-1
       call random_number(random)
       j = int( (n+1-i) * random + i )
       ! swap i and j 
       tmp = array(i)
       array(i) = array(j)
       array(j) = tmp
    end do
  end subroutine shuffle_real
  
    !> \brief Randomly shuffles an array of ints
  !>
  !>  Needed for a quicksort algorithm which doesn't hit \f$O(N^2)\f$ on a sorted array.
  !>  Uses the Fisher-Yates algorithm.
  !>  \param array An array of ints to be shuffled randomly
  subroutine shuffle_int ( array )
    implicit none

    integer, dimension(:), intent(inout) :: array
    integer :: i,j,n,tmp
    real(kind=real8) :: random
    
    ! Note we are assuming that the random number generator has been seeded correctly!
    ! Now lets shuffle by randomly picking a number on [i,n] and placing it at position array(i)
    n = size(array)
    do i = 1 , n-1
       call random_number(random)
       j = int( (n+1-i) * random + i )
       ! swap i and j 
       tmp = array(i)
       array(i) = array(j)
       array(j) = tmp
    end do
  end subroutine shuffle_int


  !> \brief Sets up the random number generator
  !>
  !> This is based on the gfortran documentation
  subroutine setup_rand()
    implicit none
    integer :: seed_size, clock, i
    integer, dimension(:), allocatable :: seed

    ! get the size of the seed
    call random_seed( size = seed_size )

    ! now allocate the seed
    allocate(seed(seed_size),stat=i)
    if ( i /= 0 ) then
       write(*,*) "Memory allocation error in setup_rand"
       stop
    end if

    ! now get the system clock time and setup the seed
    call system_clock( count = clock )
    do i = 1 , seed_size
       seed(i) = clock + 37 * i
    end do

    ! now set the seed
    call random_seed(put=seed)

    deallocate(seed)        
  end subroutine setup_rand
    
end module random
