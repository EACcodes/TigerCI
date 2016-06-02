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
!**************************************************************
  !//  UTILITIES_MOD - THIS MODULE CONTAINS UTILITY FILES WE MAY
  !//  USE.  SILLY HOUSEKEEPING STUFF GETS PUT HERE.
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************

module utilities_mod

  use global_var_mod
  use io_unit_numbers


contains

  subroutine copyArrayBool(original,copy)
    implicit none
    integer::size0,i
    logical,allocatable::original(:),copy(:)

    size0 = size(original,1)
    if(.not.allocated(copy)) then
       allocate(copy(size0))
    endif
    do i=1,size0
       copy(i) = original(i)
    enddo

  end subroutine copyArrayBool

  subroutine copyArrayInt(original,copy)
    implicit none
    integer::size0,i
    integer,allocatable::original(:),copy(:)

    size0 = size(original,1)
    if(.not.allocated(copy)) then
       allocate(copy(size0))
    endif
    do i=1,size0
       copy(i) = original(i)
    enddo

  end subroutine copyArrayInt

  subroutine copyArray(original,copy)
    implicit none
    integer::size0,i
    real(real8),allocatable::original(:),copy(:)

    size0 = size(original,1)
    if(.not.allocated(copy)) then
       allocate(copy(size0))
    endif
    do i=1,size0
       copy(i) = original(i)
    enddo

  end subroutine copyArray

  subroutine copyMatrix(original,copy)
    implicit none
    integer::size0,size1,i,j
    real(real8),allocatable::original(:,:),copy(:,:)

    size0 = size(original,1)
    size1 = size(original,2)
    if(.not.allocated(copy)) then
       allocate(copy(size0,size1))
    endif

    do j=1,size1
       do i=1,size0
          copy(i,j) = original(i,j)
       enddo
    enddo
  end subroutine copyMatrix

  subroutine copyMatrixInt(original,copy)
    implicit none
    integer::size0,size1,i,j
    integer,allocatable::original(:,:),copy(:,:)

    size0 = size(original,1)
    size1 = size(original,2)
    if(.not.allocated(copy)) then
       allocate(copy(size0,size1))
    endif

    do j=1,size1
       do i=1,size0
          copy(i,j) = original(i,j)
       enddo
    enddo
  end subroutine copyMatrixInt

  subroutine copyTensor(original,copy)
    implicit none
    integer::size0,size1,size2,i,j,k
    real(real8),allocatable::original(:,:,:),copy(:,:,:)

    size0 = size(original,1)
    size1 = size(original,2)
    size2 = size(original,3)
    if(.not.allocated(copy)) then
       allocate(copy(size0,size1,size2))
    endif

    do k=1,size2
       do j=1,size1
          do i=1,size0
             copy(i,j,k) = original(i,j,k)
          enddo
       enddo
    enddo
  end subroutine copyTensor

  subroutine allocatecheck(allocatestatus,array_string)

    !// THIS SUBROUTINE JUST CHECKS TO MAKE SURE DYNAMIC MEMORY ALLOCATION
    !// WENT O.K.  IF IT DIDN'T WE STOP THE PROGRAM AND PRINT A LITTLE 
    !// MESSAGE TELLING THE LOSER WHERE HE LOST.  THE ARRAY NAME SHOULD
    !// BE A 8 CHARACTER STRING
    !// JMD: changed that

    implicit none

    integer,intent(in)::allocatestatus
    character*(*),intent(in)::array_string

    if (allocatestatus /= 0) then
       write(ioError,*) "PROBLEM ALLOCATING ARRAY ", ARRAY_STRING
       write(ioError,*) "allocatestatus: ", allocatestatus
       stop
    endif

  end subroutine allocatecheck


  !**************************************************************
  subroutine deallocatecheck(deallocatestatus,array_string)

    !// THIS SUBROUTINE JUST CHECKS TO MAKE SURE DYNAMIC MEMORY DEALLOCATION
    !// WENT O.K.  IF IT DIDN'T WE STOP THE PROGRAM AND PRINT A LITTLE 
    !// MESSAGE TELLING THE LOSER WHERE HE LOST.  THE ARRAY NAME SHOULD
    !// BE A 8 CHARACTER STRING
    !// JMD: changed that

    implicit none

    integer,intent(in)::deallocatestatus
    character*(*),intent(in)::array_string

    if (deallocatestatus /= 0) then
       write(ioError,*) "PROBLEM DEALLOCATING ARRAY ", ARRAY_STRING
       write(ioError,*) "deallocatestatus: ", deallocatestatus
       stop
    endif

  end subroutine deallocatecheck


 
  !*****************************************************************
  !// THE FOLLOWING 4 FUNCTIONS ARE JUST USED FOR INDEXING

  integer function index2(i,j)

    !// THIS IS JUST A SIMPLE FUNCTION TO COMPUTE INDICES.
    !// I DON'T INCLUDE IT IN A SEPARATE ROUTINE BECAUSE IT IS A SMALL THING AND IT
    !// WILL BE CALLED OVER AND OVER AGAIN, AND I DON'T WANT IT TO SLOW THE CODE DOWN.
    integer,intent(in)::i,j
    index2 = 0 
    index2 = max(i,j)*(max(i,j)-1)/2 + min(i,j)

  end function index2

  integer function index2_test(i,j)

    !// THIS IS JUST A SIMPLE FUNCTION TO COMPUTE INDICES.
    !// I DON'T INCLUDE IT IN A SEPARATE ROUTINE BECAUSE IT IS A SMALL THING AND IT
    !// WILL BE CALLED OVER AND OVER AGAIN, AND I DON'T WANT IT TO SLOW THE CODE DOWN.
    integer,intent(in)::i,j
    index2_test = 0 
    index2_test = max(i,j)*(max(i,j)-1)/2 + min(i,j)

  end function index2_test

  integer function index4(i,j,k,l)

    !// THIS IS JUST A SIMPLE FUNCTION WHICH MAPS 4 DIFFERENT INTEGERS ONTO
    !// AN INDEX WHICH RUNS FROM 0 TO NC4.  THIS FUNCTION ASSUMES THAT
    !//  I > J > K > L
    integer,intent(in)::i,j,k,l
    index4 = 0 
    index4 = (i-1)*(i-2)*(i-3)*(i-4)/24 + (j-1)*(j-2)*(j-3)/6 + (k-1)*(k-2)/2 + l - 1

  end function index4

  integer function index2m1(i,j)

    !// GIVEN TWO UNEQUL INTEGERS, THIS COMPUTES AND INDEX WHERE THE PAIR
    !// (2,1) RETURNS 1, (3,1) RETURNS 2, AND SO ON.
    integer, intent(in)::i,j
    index2m1 = 0
    index2m1 = (max(i,j)-1)*(max(i,j)-2)/2 + min(i,j)

  end function index2m1


  !*****************************************************************
  function identity_matrix(size)

    integer,intent(in)::size
    real(real8), dimension(size,size)::identity_matrix
    integer::j  !// LOOP CONTROL VARIABLE

    identity_matrix = real(0.0,real8)

    do j = 1, size
       identity_matrix(j,j) = real(1.0,real8)
    enddo

  end function identity_matrix

  !*****************************************************************
  subroutine write_big_vector(vector, length, unit)

    implicit none

    !// WHAT DO YOU THINK THE ROUTINE DOES, DUMBASS?

    integer::length, unit
    real(real8)::vector(length)

    write(unit) vector

  end subroutine write_big_vector

  !*****************************************************************
  subroutine read_big_vector(vector, length, unit)

    implicit none

    !// WHAT DO YOU THINK THE ROUTINE DOES, DUMBASS?

    integer::length, unit
    real(real8)::vector(length)

    read(unit) vector

  end subroutine read_big_vector


  !*****************************************************************
  subroutine writemat(matrix,dim1,dim2,iounit,matrix_string)

    !// THIS SUBROUTINE WRITES A MATRIX TO DISK.  SHOULDN'T BE BIGGER THAN....
    !// SAY 100 X 100

    implicit none

    integer, parameter::biggest=100
    integer::dim1, dim2, iounit
    integer::i,j

    real(real8)::matrix(dim1,dim2)
    character(8)::matrix_string

    if (dim1 > biggest .or. dim2 > biggest) then
       write(ioOutput,*) "Matrix to big in writemat; adjust biggest"
       write(ioOutput,*) dim1, dim2, biggest
       stop
    endif

    write(iounit,*)
    write(iounit,*) matrix_string

    do i = 1, dim1
       write(iounit,110) (matrix(i,j), j = 1, dim2)
    enddo
    write(iounit,*)

110 format(100f12.5)    


  end subroutine writemat



  !> \brief Returns true if the two arrays are of equal length and contain the same values
  !> \param a The first array
  !> \param b The second array
  logical function arraysAreEqual(a,b)
    implicit none
    ! inputs
    integer, dimension(:) ,intent(in) :: a,b
    ! local variables
    integer :: a_size, b_size, n, i
    a_size = size(a)
    b_size = size(b)
    if ( .not. a_size == b_size ) then
       arraysAreEqual = .False.
       return
    end if
    n = a_size
    arraysAreEqual = .True.
    do i = 1, n
       if ( a(i) /= b(i) ) then
          arraysAreEqual = .False.
          return
       end if
    end do
  end function arraysAreEqual

  !> because sometimes you really want the code to stop 
  !> and (hopefully) produce a stack trace
  subroutine cause_floating_point_exception()
    implicit none
    integer :: i
    call flush(6)
    i = i / 0
    i = i / i
    write(*,*) "boom ?", i
  end subroutine cause_floating_point_exception
  
  
  ! the following routines are simple deallocators for arrays which may or
  ! may not be allocated
  
    subroutine try_deallocate_int_1D(array, str)
      implicit none
      integer, allocatable, intent(inout) :: array(:)
      integer :: ierr
      character(len=*), intent(in) :: str
      if (allocated(array)) then 
        deallocate(array, stat=ierr)
        call deallocatecheck(ierr, str)
      end if
  end subroutine
  
  subroutine try_deallocate_real_1D(array, str)
      implicit none
      real(kind=real8), allocatable, intent(inout) :: array(:)
      integer :: ierr
      character(len=*), intent(in) :: str
      if (allocated(array)) then 
        deallocate(array, stat=ierr)
        call deallocatecheck(ierr, str)
      end if
  end subroutine
  
    subroutine try_deallocate_int_2D(array, str)
      implicit none
      integer, allocatable, intent(inout) :: array(:,:)
      integer :: ierr
      character(len=*), intent(in) :: str
      if (allocated(array)) then 
        deallocate(array, stat=ierr)
        call deallocatecheck(ierr, str)
      end if
  end subroutine
  
  subroutine try_deallocate_real_2D(array, str)
      implicit none
      real(kind=real8), allocatable, intent(inout) :: array(:,:)
      integer :: ierr
      character(len=*), intent(in) :: str
      if (allocated(array)) then 
        deallocate(array, stat=ierr)
        call deallocatecheck(ierr, str)
      end if
  end subroutine
  
    subroutine try_deallocate_real_3D(array, str)
      implicit none
      real(kind=real8), allocatable, intent(inout) :: array(:,:,:)
      integer :: ierr
      character(len=*), intent(in) :: str
      if (allocated(array)) then 
        deallocate(array, stat=ierr)
        call deallocatecheck(ierr, str)
      end if
  end subroutine

end module utilities_mod
