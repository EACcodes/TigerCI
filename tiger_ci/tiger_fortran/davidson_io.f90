! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> Fortran before 2003 didn't technically allow anything other
!> than a 4 byte integer for the size of a record (recl)
!>
!> As of 8/2013 ifort supports larger integers, but gfortran
!> doesn't. As a result we need this code write out blocks of 
!> data larger than ~2Gb
!>
!>
!> Note (as requested by Johannes): if the day comes to pass 
!> that gfortran supports the necessary part of the 2003 standard
!> (8 byte/64 bit interfaces for 'direct' mode binary I/O) then 
!> feel free to toss this and all the C++ code out and go back
!> to using fortran


module davidson_io_mod
  use global_var_mod
  use utilities_mod
  use c_io_finterface
  implicit none

    character(len=2), parameter :: cfile='c', Hcfile='Hc', Gcfile='Gc'
    character(len=8), parameter :: ci_final="ci_final"

contains

  subroutine davidson_io_setup()   
    call setup_vector_io_with_c(scratch_directory // ci_final)
    call setup_vector_io_with_c(scratch_directory // cfile)
    call setup_vector_io_with_c(scratch_directory // Hcfile)
    call setup_vector_io_with_c(scratch_directory // Gcfile)
  end subroutine davidson_io_setup
  
  subroutine davidson_io_write(vec, filename, vec_num)
    ! inputs
    real(kind=real8), dimension(:) :: vec
    character(len=*), intent(in) :: filename 
    integer, intent(in) :: vec_num
    call write_vector_with_c(scratch_directory // filename, vec, size(vec), vec_num)
  end subroutine davidson_io_write

  subroutine davidson_io_read(vec, filename, vec_num)
    ! inputs
    real(kind=real8), dimension(:) :: vec
    integer, intent(in) :: vec_num
    character(len=*), intent(in) :: filename 
    call read_vector_with_c(scratch_directory // filename, vec, size(vec), vec_num)
  end subroutine davidson_io_read      

end module davidson_io_mod
