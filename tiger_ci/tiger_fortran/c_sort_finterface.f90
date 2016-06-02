! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> We are going to link to the C++ STL for sorting capabilities
!> (because apparently nobody on the fortran standards committee 
!>  thinks that FORTRAN developers ever need to sort things ....)
!>
!> For references on combining C/Fortran see:
!> http://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html 
!> http://stackoverflow.com/tags/fortran-iso-c-binding/info
!> http://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Subroutines-and-Functions.html
!>
!> All sorts sort into ascending order

module c_sort_finterface
  use global_var_mod
  use iso_c_binding
  implicit none
  
  private
  public :: sort_real_array, sort_int_array, sort_real_array_with_index, sort_int_array_with_index

  interface 
     subroutine c_sort_real_array(x,n) bind(c)
       use iso_c_binding
       implicit none
       real(kind=C_DOUBLE), intent(inout) :: x(*)
       integer(kind=C_INT64_T), intent(in) :: n
     end subroutine 
  end interface
  
    interface 
     subroutine c_sort_int_array(x,n) bind(c)
       use iso_c_binding
       implicit none
       integer(kind=C_INT64_T), intent(inout) :: x(*)
       integer(kind=C_INT64_T), intent(in) :: n
     end subroutine 
  end interface
  
  interface 
      subroutine c_sort_real_array_with_index(x,index,n) bind(c)
          use iso_c_binding
          implicit none
          real(kind=C_DOUBLE), intent(inout) :: x(*)
          integer(kind=C_INT64_T), intent(inout) :: index(*)
          integer(kind=C_INT64_T), intent(in) :: n
      end subroutine
  end interface
  
  interface 
      subroutine c_sort_int_array_with_index(x,index,n) bind(c)
          use iso_c_binding
          implicit none
          integer(kind=C_INT64_T), intent(inout) :: x(*)
          integer(kind=C_INT64_T), intent(inout) :: index(*)
          integer(kind=C_INT64_T), intent(in) :: n
      end subroutine
  end interface

  contains
  
  !> \brief Sort real array x
  !> \param x array to be sorted
  subroutine sort_real_array(x)
      implicit none
      real(kind=real8), intent(inout), target :: x(:)
      integer(kind=C_INT64_T) :: n
      n = size(x)
      call c_sort_real_array(x,n)
  end subroutine

  !> \brief Sort integer array x
  !> \param x array to be sorted
  subroutine sort_int_array(x)
      implicit none
      integer(kind=C_INT64_T), intent(inout), target :: x(:)
      integer(kind=C_INT64_T) :: n
      n = size(x)
      call c_sort_int_array(x,n)
  end subroutine
  
  !> \brief Sorts a real array x and it's index vector (the index is assumed to be base 1)
  !> \param x array to be sorted
  !> \param index array to be permuted according to the sort of x
  !>
  !> This routine is very specific. The index vector needs to be 
  !> 1, 2, 3, 4, 5, 6, ...
  !> on entry. When its returned it will indicate the permutation created by sorting
  !> for instance 
  !> 2, 1, 5, ...
  !> means that the lowest element (we sort into ascending order) was originally the second 
  !> element. If the index vector is anything other than a sequence starting with 1 then  
  !> the sort is going to fail. Horribly Fail. 
  subroutine sort_real_array_with_index(x, index)
      implicit none
      real(kind=real8), intent(inout), target :: x(:)
      integer(kind=C_INT64_T), intent(inout), target :: index(:)
      integer(kind=C_INT64_T) :: n
      n = size(x)
      if ( n /= size(index) ) then 
          write(*,*) "Call to sort with permute failed"
          write(*,*) "because the array to be sorted " 
          write(*,*) "and the permutate array are different sizes!!!"
          stop
      end if
      call c_sort_real_array_with_index(x,index,n)
  end subroutine 

  !> \brief Sorts a integer array x and it's index vector (the index is assumed to be base 1)
  !> \param x array to be sorted
  !> \param index array to be permuted according to the sort of x
  !>
  !> This routine is very specific. The index vector needs to be 
  !> 1, 2, 3, 4, 5, 6, ...
  !> on entry. When its returned it will indicate the permutation created by sorting
  !> for instance 
  !> 2, 1, 5, ...
  !> means that the lowest element (we sort into ascending order) was originally the second 
  !> element. If the index vector is anything other than a sequence starting with 1 then  
  !> the sort is going to fail. Horribly Fail. 
  subroutine sort_int_array_with_index(x, index)
      implicit none
      integer(kind=C_INT64_T), intent(inout), target :: x(:)
      integer(kind=C_INT64_T), intent(inout), target :: index(:)
      integer(kind=C_INT64_T) :: n

      n = size(x)
      if ( n /= size(index) ) then 
          write(*,*) "Call to sort with permute failed"
          write(*,*) "because the array to be sorted " 
          write(*,*) "and the permutate array are different sizes!!!"
          stop
      end if
      call c_sort_int_array_with_index(x,index,n)
  end subroutine 
end module  
