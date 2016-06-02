! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> Basic reading and writing of real(kind=real8) vectors using C routines
!>
!> For references on combining C/Fortran see:
!> http://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html 
!> http://stackoverflow.com/tags/fortran-iso-c-binding/info
!> http://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Subroutines-and-Functions.html
!>
!> Before you read/write a vector you need to run the setup_vector_io routine
!> otherwise weird things can happen!

module c_io_finterface
  use global_var_mod
  use iso_c_binding
  implicit none


  interface 
     subroutine c_setup_vector_io(filename) bind(c)
       use iso_c_binding
       implicit none
       character(kind=c_char), intent(in) :: filename(*)
     end subroutine c_setup_vector_io
  end interface

  interface
     subroutine c_close_files() bind(c)
     end subroutine c_close_files
  end interface

  interface 
     subroutine c_write_vector(filename, position, vector, length) bind(c)
       use iso_c_binding
       implicit none
       character(kind=c_char), intent(in) :: filename(*)
       integer(kind=C_INT64_T), intent(in) :: position
       real(kind=C_DOUBLE), intent(in) :: vector(*)
       integer(kind=C_INT64_T), intent(in) :: length
     end subroutine c_write_vector
  end interface

  interface 
     subroutine c_read_vector(filename, position, vector, length) bind(c)
       use iso_c_binding
       implicit none
       character(kind=c_char), intent(in) :: filename(*)
       integer(kind=C_INT64_T), intent(in) :: position
       real(kind=C_DOUBLE), intent(out) :: vector(*)
       integer(kind=C_INT64_T), intent(in) :: length
     end subroutine c_read_vector
  end interface

  
  interface 
     subroutine c_write_ints(filename, position, vector, length) bind(c)
       use iso_c_binding
       implicit none
       character(kind=c_char), intent(in) :: filename(*)
       integer(kind=C_INT64_T), intent(in) :: position
       integer(kind=C_INT64_T), intent(in) :: vector(*)
       integer(kind=C_INT64_T), intent(in) :: length
     end subroutine c_write_ints
  end interface

  interface 
     subroutine c_read_ints(filename, position, vector, length) bind(c)
       use iso_c_binding
       implicit none
       character(kind=c_char), intent(in) :: filename(*)
       integer(kind=C_INT64_T), intent(in) :: position
       integer(kind=C_INT64_T), intent(out) :: vector(*)
       integer(kind=C_INT64_T), intent(in) :: length
     end subroutine c_read_ints
  end interface

  
contains

  subroutine clean_up_io_with_c()
    call c_close_files()
  end subroutine clean_up_io_with_c

  subroutine setup_vector_io_with_c(name)
    character(len=*) :: name
    ! null terminate C strings   
    call c_setup_vector_io(name//C_NULL_CHAR)
  end subroutine setup_vector_io_with_c

  subroutine read_vector_with_c(name, vector, vector_length, position)
    integer :: vector_length
    character(len=*) :: name 
    real(kind=real8), intent(out), target :: vector(vector_length)
    integer, intent(in) :: position

    integer(kind=C_INT64_T) :: c_length, c_position

    c_length = int(vector_length, C_INT64_T)
    c_position = int(position, C_INT64_T)

    ! null terminate C strings
    call c_read_vector(name//C_NULL_CHAR, c_position, vector, c_length)
  end subroutine read_vector_with_c

  subroutine write_vector_with_c(name, vector, vector_length, position)
    integer :: vector_length
    character(len=*) :: name 
    real(kind=real8), intent(in), target :: vector(vector_length)
    integer, intent(in) :: position

    integer(kind=C_INT64_T) :: c_length, c_position

    c_length = int(vector_length, C_INT64_T)
    c_position = int(position, C_INT64_T)

    ! null terminate C strings
    call c_write_vector(name//C_NULL_CHAR, c_position, vector, c_length)
  end subroutine write_vector_with_c

  subroutine read_ints_with_c(name, vector, vector_length, position)
    integer :: vector_length
    character(len=*) :: name 
    integer(kind=C_INT64_T), intent(out), target :: vector(vector_length)
    integer, intent(in) :: position

    integer(kind=C_INT64_T) :: c_length, c_position

    c_length = int(vector_length, C_INT64_T)
    c_position = int(position, C_INT64_T)

    ! null terminate C strings
    call c_read_ints(name//C_NULL_CHAR, c_position, vector, c_length)
  end subroutine read_ints_with_c

  subroutine write_ints_with_c(name, vector, vector_length, position)
    integer :: vector_length
    character(len=*) :: name 
    integer(kind=C_INT64_T), intent(in), target :: vector(vector_length)
    integer, intent(in) :: position

    integer(kind=C_INT64_T) :: c_length, c_position

    c_length = int(vector_length, C_INT64_T)
    c_position = int(position, C_INT64_T)

    ! null terminate C strings
    call c_write_ints(name//C_NULL_CHAR, c_position, vector, c_length)
  end subroutine write_ints_with_c

  
end module c_io_finterface
