! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief Provides basic timing functionality using F90 intrinsics

module fortran_timing
  use global_var_mod

  type clock
     real(kind=real8) :: cpuStart
     integer :: wallStart, clock_rate, clock_max
  end type clock

contains


  !> \brief Makes a timer
  !> \param timer A clock data structure 
  subroutine start_clock(timer)
    implicit none
    type(clock), intent(out) :: timer

    ! First the cpu_time
    call cpu_time(timer%cpuStart) 

    ! Next the wall time
    call system_clock(timer%wallStart, timer%clock_rate, timer%clock_max)
  end subroutine start_clock

  !> \brief Prints out the amount of time that has passed in Wall/CPU since the timer was started
  !> \param timer A clock data structure
  !> \param title A string describing what was being measured
  subroutine print_clock(timer, title)
    implicit none
    type(clock) :: timer
    character(len=*), intent(in) :: title
    real :: cpu, wallTotal
    integer :: wall

    ! First the cpu time
    call cpu_time(cpu)
    cpu = cpu - timer%cpuStart

    ! Next the wall time
    call system_clock(wall, timer%clock_rate, timer%clock_max)
    wallTotal = real(wall - timer%wallStart) / real(timer%clock_rate)

    ! Print the results
800 format( A, "    wall time(s): ", E10.3, " cpu time(s): ", E10.3)
    write(*,800) title, wallTotal, cpu
    !write(*,'(A,A,A,E22.15,A,E22.15)') ' ', title, " wall time(s)=", wallTotal, " cpu time(s) = " , cpu
  end subroutine print_clock

  real(kind=real8) function get_clock_wall_time(timer)
    implicit none
    type(clock) :: timer
    real :: cpu, wallTotal
    integer :: wall

    ! First the cpu time
    call cpu_time(cpu)
    cpu = cpu - timer%cpuStart

    ! Next the wall time
    call system_clock(wall, timer%clock_rate, timer%clock_max)
    wallTotal = real(wall - timer%wallStart) / real(timer%clock_rate)

    get_clock_wall_time = wallTotal
  end function get_clock_wall_time


end module fortran_timing
