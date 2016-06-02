! Copyright (c) 1998-2004, University of California, Los Angeles, Emily A. Carter
!                       2004-2016, Princeton University, Emily A. Carter
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of
!     conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!     of conditions and the following disclaimer in the documentation and/or other
!     materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!     used to endorse or promote products derived from this software without specific
!     prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
! CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
! OF THE POSSIBILITY OF SUCH DAMAGE.
!> \brief This module stores lists of allowed excitations from internal orbitals (i,j) to virtual orbitals (a)
!> \author Derek Walter
!> \author David Krisiloff
!>
!> This module is the heart of the TOV approximation, allowing us to truncate virtual excitations. An internal configuration maps
!> to the allowed_virtuals list by allowed_virtuals(N1_index(weight)) (or N2 if its an N-2 path). The mapping from N1/N2_index to
!> allowed virtuals is not 1 to 1. Each array in allowed_virtuals is unique: multiple paths map to the same
!> allowed_virtuals array. For the nonlocal case there is just one allowed_virtual array (with all the virtual orbitals in it). 
!>
!> Note that this scheme is very different from the previous. We are no longer using the internal holes (i,j) to index anything.
!> The previous scheme breaks down in complicated MR cases were one configuration could be viewed as an excitation from multiple 
!> references, each with their own (i,j) hole pair.

module allowed_virtuals_mod

  use global_var_mod
  use molecule_var_mod
  use graph_var_mod
  use utilities_mod
  use ci_utilities_mod
  use sort_utils
  use c_sort_finterface

  implicit none

  ! default everything to private ... you must use one of the public APIs
  private
  public :: num_allowed_virtuals, get_virtuals, get_virtuals_and_index, init_allowed_virtuals, &
            add_new_virtual_list, debug_allowed_virtuals, zero_allowed_virtuals, delete_allowed_virtuals
          
  ! A list of virtual orbitals
  type virt_list
     integer, dimension(:), allocatable:: virtuals
     integer :: num_virt
  end type virt_list

  ! Module global (but private!) variables which store the allowed virtuals
  type(virt_list), pointer, dimension(:) :: allowed_virtuals 
  integer :: curr_size
  integer :: max_size
  logical :: was_init = .false.
    
  ! Module global (but private!) indices for allowed virtuals
  integer, dimension(:), allocatable, target :: N1_index              ! maps weight -> position on allowed virtuals for (N-1) paths
  integer, dimension(:), allocatable, target :: N2_index              ! maps weight -> position on allowed virtuals for (N-2) paths
  

contains

!*************************************************************************************
!> \brief Deallocates the allowed virtuals list
subroutine delete_allowed_virtuals()
    implicit none
    integer :: ierr, i 
    do i =1, max_size 
        deallocate(allowed_virtuals(i)%virtuals, stat=ierr)
        call deallocatecheck(ierr, "delete_allowed_virtuals")
    end do
    curr_size = 0
    max_size = 0 
    was_init = .false.
    
    deallocate(N1_index, N2_index, stat=ierr)
    call deallocatecheck(ierr, "delete_allowed_virtuals 2")
end subroutine

 !*************************************************************************************
  !> \brief Allocates space for an allowed_virtuals list
  subroutine init_allowed_virtuals()
    implicit none

    ! local variables
    integer :: i, ierr
    integer, parameter :: start_size = 2

    ! Skip this if we are already allocated
    if ( was_init ) return
    was_init = .true.

    ! Allocate the allowed virtuals
    allocate(allowed_virtuals(start_size), stat=ierr)
    call allocatecheck(ierr, "allowed_virtuals")
    do i = 1, start_size
       allocate(allowed_virtuals(i)%virtuals(num_external),stat=ierr)
       call allocatecheck(ierr, "member of allowed_virtuals")
       allowed_virtuals(i)%virtuals = 0 
       allowed_virtuals(i)%num_virt = 0 
    end do
    
    ! Allocate the index
    allocate(N1_index(abs(y0(vertex(num_internal,num_elec-1)))), stat=ierr)
    call allocatecheck(ierr, "N1_index")
    allocate(N2_index(abs(y0(vertex(num_internal,num_elec-2)))), stat=ierr)
    call allocatecheck(ierr, "N1_index")
    N1_index = -1
    N2_index = -1

    ! Set the size of our data structures
    curr_size = 0 
    max_size  = 2
  end subroutine init_allowed_virtuals

  !*************************************************************************************
  !> \brief Zeros the allowed virtuals (doesn't deallocate their memory!)
  subroutine zero_allowed_virtuals()
    implicit none

    ! local variables
    integer ::i 

    do i = 1, curr_size
       allowed_virtuals(i)%virtuals = 0
       allowed_virtuals(i)%num_virt = 0 
    end do
    
    N1_index = 0
    N2_index = 0

    curr_size = 0

  end subroutine zero_allowed_virtuals
  

   !*************************************************************************************
  !> Adds a new list of allowed virtuals for an internal configuration
  !> \param weight The weight of the internal configurations
  !> \param elec   The number of electrons in the internal configuration
  !> \param num_virt The number of allowed virtuals
  !> \param virtuals The allowed virtuals
  !>
  !> Important Note! These list of orbitals must be sorted! (the SGGA assumes it)
  subroutine add_new_virtual_list(weight, elec, num_virt, virtuals)
      implicit none
      
      ! inputs
      integer, intent(in) :: weight, elec, num_virt
      integer, intent(in) :: virtuals(:)
      
      ! local variables
      integer ::  i
      integer, dimension(:), pointer :: index
      logical :: unique_list
      integer, dimension(:), allocatable :: tmp_virtuals

    
    ! Set the correct index
    if ( elec == num_elec-2 ) then 
        index => N2_index
    else if ( elec == num_elec-1 ) then
        index => N1_index
    else
        write(*,*) "FATAL ERROR ... add_new_virtual_list was called with an internal path of N electrons"
        write(*,*) "If all the electrons are in the internal space there are no virtual orbitals."
        stop
    end if

    ! sort the list of virtuals into ascending order
    allocate(tmp_virtuals(size(virtuals)))
    tmp_virtuals = virtuals
    ! this call creates a temporary, but it is necessary. The ends of this array may be padded with 0s
    ! (on [num_virt+1:size(virtuals)]) and when we sort those get placed at the beginning of the array,
    ! where we don't want them
    
    !call insertion_sort_int(tmp_virtuals(1:num_virt))
    !
    !
    call sort_int_array(tmp_virtuals(1:num_virt))
    
    
    ! check to see if this list already exists on the allowed_virtuals list
    unique_list = .true. 
    do i = 1 , curr_size
       if ( arraysAreEqual(tmp_virtuals, allowed_virtuals(i)%virtuals) ) then
          unique_list = .false.
          exit
       end if
    end do

    if ( .not. unique_list ) then
        index(weight) = i 
        return
    end if
    
    ! ok we know that this is a unique list of virtuals ... lets add it
    if ( curr_size == max_size ) call double_allowed_virtuals()

    curr_size = curr_size + 1
    allowed_virtuals(curr_size)%virtuals(1:num_virt) = tmp_virtuals(1:num_virt)
    allowed_virtuals(curr_size)%num_virt = num_virt
    index(weight) = curr_size

    deallocate(tmp_virtuals)
  end subroutine add_new_virtual_list

  
  !*************************************************************************************
  !> \brief Doubles the size of the current allowed_virtuals list 
  subroutine double_allowed_virtuals()
      implicit none
      
      ! local variables
      integer :: old_size, new_size
      integer :: i, ierr
      type(virt_list), pointer, dimension(:) :: new_list
      
      old_size = curr_size
      new_size = 2 * old_size

      ! allocate the new space
      allocate(new_list(new_size), stat=ierr)
      call allocatecheck(ierr, "new allowed virtuals list")
      do i = 1, new_size
          allocate(new_list(i)%virtuals(num_external),stat=ierr)
          call allocatecheck(ierr, "member of new_list")
          new_list(i)%virtuals = 0 
          new_list(i)%num_virt = 0 
      end do
      
      ! copy the old data to the new space
      do i = 1, old_size
          new_list(i)%virtuals(:) = allowed_virtuals(i)%virtuals(:)
          new_list(i)%num_virt    = allowed_virtuals(i)%num_virt
      end do
      
      ! deallocate the old struct and point allowed_virtuals to the new memory
      do i = 1, old_size
          deallocate(allowed_virtuals(i)%virtuals, stat=ierr)
          call deallocatecheck(ierr, ",member of allowed_virtuals")
      end do
      deallocate(allowed_virtuals, stat=ierr)
      call deallocatecheck(ierr, "allowed_virtuals")
      allowed_virtuals => new_list
      
      ! update max_size
      max_size = new_size
      
  end subroutine double_allowed_virtuals

  !*************************************************************************************
  !> \brief Given the weight of a configuration and the number of elec (N-1) or (N-2) returns the number of allowed excitations
  !> \param weight  The weight (from the SGGA graph) of an internal configuration
  !> \param state   "D","d" for N-2, "S","s" for N-1   (Yes I agree this could be confusing ... sigh)
  function num_allowed_virtuals(weight,state)
      implicit none
      
      ! inputs
    character::state      
    integer::weight   
    
    ! local variables 
    integer::i
    integer::num_allowed_virtuals

    i = 0  
    if (state == "S" .or. state == "s") then
       i = N1_index(weight)
    elseif (state == "D" .or. state == "d") then
       i = N2_index(weight)
    endif

    if (i < 0) then
       num_allowed_virtuals = 0
       return
    else 
       num_allowed_virtuals = allowed_virtuals(i)%num_virt
       return
    endif

    return
  end function num_allowed_virtuals

  !*************************************************************************************
  !>\brief Given an internal configuration return the allowable virtuals to excite to
  !>\param weight  The weight of the internal configuration
  !>\param state   "S" or "D" -> one or two electrons in the virtual space
  !>\param virtuals An integer array on exit containg the allowed virtuals
  subroutine get_virtuals(weight, state, virtuals)
      implicit none
      
      ! inputs
      character, intent(in)::state      !// "D","d" for N-2, "S","s" for N-1
      integer, dimension(:),intent(inout) :: virtuals
      integer, intent(in)::weight       !// WEIGHT OF THE INTERNAL CONFIG
    
      ! local variables
      integer::i        !// INDEX OF THE HOLE PAIR
      integer::num_virt     !// NUMBER OF VIRTUALS FOR THIS INTERNAL CONFIG

      i = 0 
      if (state == "S" .or. state == "s") then
         i = N1_index(weight)
      elseif (state == "D" .or. state == "d") then
         i = N2_index(weight)
      endif
    
      ! weird error condition ... not sure why this isn't caught earlier ....
      if ( i < 0 ) return 

      num_virt = allowed_virtuals(i)%num_virt

      virtuals(1:num_virt) = allowed_virtuals(i)%virtuals(1:num_virt)
  end subroutine get_virtuals

  !*************************************************************************************
  !>\brief Given an internal configuration return the allowable virtuals to excite to
  !>\param weight  The weight of the internal configuration
  !>\param state   "S" or "D" -> one or two electrons in the virtual space
  !>\param virtuals An integer array on exit containg the allowed virtuals
  !>\param i The index for the allowed_virtuals. Please note that the client code must check that it is both S or D!
  subroutine get_virtuals_and_index(weight, state, virtuals, i)
      implicit none
      ! inputs
    integer, intent(in)::weight      
    integer, intent(out)::i        
 integer, dimension(:),intent(out) :: virtuals
    character, intent(in)::state      
    
    ! local variables
    integer::num_virt     

    i = 0 
    if (state == "S" .or. state == "s") then
       i = N1_index(weight)
    elseif (state == "D" .or. state == "d") then
       i = N2_index(weight)
    endif
    
    ! weird error condition ... not sure why this isn't caught earlier ....
    if ( i < 0 ) return 

    num_virt = allowed_virtuals(i)%num_virt
    virtuals(1:num_virt) = allowed_virtuals(i)%virtuals(1:num_virt)
  end subroutine get_virtuals_and_index
  
  !> \brief Returns the average allowed virtual length
  integer function average_num_virt()
      implicit none
      
       ! local variables
      integer :: i, ierr, p, total 
      integer, dimension(:), allocatable :: count
      
      allocate(count(curr_size),stat=ierr)
      call allocatecheck(ierr, "count in average_num_virt")
      count = 0
      
      ! count how many paths point to each entry in allowed_virtuals
      do i = 1, size(N1_index)
          p = N1_index(i)
          if (p<=0) cycle
          count(p) = count(p) + 1
      end do
      do i = 1, size(N2_index)
          p = N2_index(i)
          if (p<=0) cycle
          count(p) = count(p) + 1
      end do

      total = sum(count)
      
      do i = 1, curr_size
          average_num_virt = average_num_virt + allowed_virtuals(i)%num_virt * count(i) 
      end do
      
      average_num_virt = int( real(average_num_virt,real8)/ total )
  end function  average_num_virt

  !> \brief Writes out the allowed virtual data structures for debugging                                                                                       
  subroutine debug_allowed_virtuals
    implicit none

    ! local variables                                                                                                                                          
    integer :: i 
    character(len=200) :: format_str
    integer :: empty(num_external)

    ! build the format for the table
    write(format_str,'(a,I20,a)') "('    Weight =',I5,1x,'Index =',I5,1x,'Virtuals =',1x,",num_external,'I5)'
    empty = -1

    ! write out the tables
    write(*,*)
    write(*,*)
    write(*,*) "    Allowed Virtual Debug Requested ...."
    write(*,*)
    write(*,*) "    Table 1 : (N-1) paths"
    write(*,*) "    ---------------------"
    do i = 1 , size(N1_index)
        if ( N1_index(i) > 0 ) then
        write(*,format_str) i, N1_index(i), allowed_virtuals(N1_index(i))%virtuals
        else
            write(*,format_str) i, N1_index(i), empty
        end if
    end do
    write(*,*) 
    write(*,*) "    Table 2 : (N-2) paths"
    write(*,*) "    ---------------------"
    do i = 1 , size(N1_index)
        if ( N2_index(i) > 0 ) then
       write(*,format_str) i, N2_index(i), allowed_virtuals(N2_index(i))%virtuals
        else
            write(*,format_str) i, N2_index(i), empty
        end if
    end do
  end subroutine debug_allowed_virtuals

end module allowed_virtuals_mod

