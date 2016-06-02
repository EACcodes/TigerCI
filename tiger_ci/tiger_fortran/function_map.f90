! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
! Maps number of basis functions to angular momentum, and assorted related tasks
! Francis Ricci - July 2015
!

module function_map
  implicit none

contains
  

  ! Convert the number of basis functions to angular momentum
  subroutine funcs_to_am(nFuncs, am_i)
    implicit none

    integer(kind=8), intent(in) :: nFuncs
    integer(kind=8), intent(out) :: am_i

    if (nFuncs == 1) then
       am_i = 0
    else if (nFuncs == 3) then
       am_i = 1
    else if (nFuncs == 6 .or. nFuncs == 5) then
       am_i = 2
    else if (nFuncs == 10 .or. nFuncs == 7) then
       am_i = 3
    else
       write(*,*) "Error in funcs to am conversion!"
       flush(6)
       stop
    end if

  end subroutine funcs_to_am
              

  ! Find the angular momentum of a given basis function, i.
  ! Return value through am_i
  subroutine get_am(i,nFuncs,nShell,am_i)
    implicit none

    integer(kind=8), intent(in) :: i
    integer(kind=8), intent(in) :: nFuncs(:)
    integer(kind=8), intent(in) :: nShell
    integer(kind=8), intent(out) :: am_i

    integer(kind=8) :: num_funcs, index

    am_i = -1
    num_funcs = 0

    do index = 1, nShell
       num_funcs = num_funcs + nFuncs(index)
       if (num_funcs >= i) then
          call funcs_to_am(nFuncs(index),am_i)
          exit
       end if
    end do
      
    if (am_i == -1) then
       write(*,*) "Error in am conversion!"
       flush(6)
       stop
    end if

  end subroutine get_am

  function diag_address(i, j)
    implicit none

    integer(kind=8), intent(in) :: i,j
    integer(kind=8) :: diag_address
    integer(kind=8) :: tmp, i_loc, j_loc

    i_loc = i
    j_loc = j

    if (i < j) then
       tmp = i_loc
       i_loc = j_loc
       j_loc = tmp
    end if

    diag_address = i_loc * (i_loc - 1) / 2 + j_loc

  end function diag_address

end module function_map
