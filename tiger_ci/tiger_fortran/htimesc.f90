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
module SGGA

  logical :: first_time = .true.
  
  contains

!*****************************************************************
!//  HTIMESC - THIS IS THE DRIVER ROUTINE FOR THE MULTIPLICATION OF THE
!//  HAMILTONIAN TIMES THE CI VECTOR IN DIRECT MODE.
!// 
!//  WRITTEN BY DEREK WALTER, 1999
!//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS.
!*****************************************************************
  
  subroutine htimesc(civec, sigmavec, data)

    use global_var_mod
    use two_seg_mod_3
    use three_and_four_seg_mod_1
    use three_and_four_seg_mod_2
    use time_var_mod
    use make_all_integrals
    use diagtimesc_mod
    use blocked_locks_mod
    use sgga_struct_mod
    use IOBuffer
    use fortran_timing

    implicit none

    type(blockedLockVectorType) :: civec, sigmavec
    real(real8)::zero
    real(real8)::CPU1, CPU2, tim1, tim2, tim3
    integer::iteration,status
    type(clock) :: timer
    type(sgga_struct) :: data
    
    zero = real(0.0, real8)
    cpu2 = 0.0
    cpu1 = 0.0
    tim2 = 0.0
    tim3 = 0.0
    tim1 = 0.0

    ! Its now the responsibility of the SGGA code to figure out how many times it has been run
    if ( first_time ) then 
       first_time = .false.
       iteration = 1 
    else
       iteration = 2
    end if

    !// FIRST, LETS TAKE CARE OF THE PURELY DIAGONAL ELEMENTS AND THE THE
    !// OFFDIAGONAL ELEMENTS ASSOCIATED WITH DIFFERENT SPIN FUNCTIONS COMING FROM
    !// THE SAME SPATIAL FUNCTION.
    write(ioOutput,*)
    write(ioOutput,*) "* Entering diag"
    call flush(ioOutput)
    call start_clock(timer)
    call diagtimesc(civec, sigmavec, data%loc_scr)
    call print_clock(timer, "  + diag time")
    diagtimesc_time = diagtimesc_time + get_clock_wall_time(timer)
    write(ioOutput,*) " => Did diagonal elements"
    call flush(ioOutput)

    !// SECOND, THE TWO SEGMENT LOOPS.  THESE ARE THE TOUGHEST FROM 
    !// A LOGICAL PERSPECTIVE, SO THEY TAKE AN INORDINATE AMOUNT
    !// OF TIME
    write(ioOutput,*) "* Entering two seg"
    call flush(ioOutput)
    call start_clock(timer)
    call two_seg_driver(civec, sigmavec, iteration, data%cho_data, data%loc_scr, data%twoVars)
    call print_clock(timer, "  + two seg ")
    two_seg_time = two_seg_time + get_clock_wall_time(timer)
    !write(ioOutput,*) "   + two seg time: ",two_seg_time
    write(ioOutput,*) " => Did two seg loops"
    call flush(ioOutput)

    !// THIRD, THE THREE AND FOUR SEGMENT LOOPS
    write(ioOutput,*) "* Entering driver 1"
    call flush(ioOutput)
    call start_clock(timer)
    call three_and_four_seg_driver_1(civec, sigmavec, iteration,  data%cho_data,  data%loc_scr,  data%mod1vars)
    write(ioOutput,*) " => Did first part of 3,4 seg"


    write(ioOutput,*) "* Entering driver 2"
    call flush(ioOutput)
    call three_and_four_seg_driver_2(civec, sigmavec, iteration,  data%cho_data,  data%loc_scr,  data%mod2vars)
    write(IoOutput,*) " => Did second part of 3,4 seg"
    three_and_four_seg_time = three_and_four_seg_time + get_clock_wall_time(timer)
    call flush(ioOutput)
    
    if(iteration .eq. 1 .and. reference_ci_flag .ne. 1 .and. valence_ci_flag .ne. 1) then
       ! close our cd buffer after all integrals have been prepared
       write(ioOutput,*) " Trying to close buffer for CD vectors." 
       flush(ioOutput)
       if(.not. (integralDirect .or. fullyIntegralDirect .or. directFourInternal)) then
         if(for_buf_storeCDVecs) then
           call for_double_buf_closepool(for_buf_cd_poolID)
         else
           call for_double_buf_removepool(for_buf_cd_poolID)
         endif
         if(cdVecsInMemory) then
           deallocate(data%cho_data%cho_vectors,data%cho_data%cho_norms,stat=status)
           call deallocatecheck(status,"gtimesc CD data")
         endif
         ! resize the integral buffer to take the space from the CD vectors
         write(ioOutput,*) " Resizing integral buffer to maximum size ",(for_buf_maxmemInts+for_buf_maxmemCD)
         flush(ioOutput)
         call for_double_buf_changebuffersize(for_buf_int_poolID,(for_buf_maxmemInts+for_buf_maxmemCD)/(for_buf_blocksizeInts*real8))
       endif
    endif


  end subroutine htimesc
  
  !> \brief Resets the SGGA to think that it is at iteration 1 again
  !>
  !> This IS necessary after running a reference CI before going on to another ACPF/SDCI calculation
  subroutine SGGA_reboot()
#ifdef TIGER_USE_SLOW_OMP
      use arraylist_integer
      use three_and_four_seg_mod_2,only:omp_offsets_four_icount,omp_offsets_four_nodecount
      use global_var_mod
#endif
      
      implicit none
      first_time = .true.
      
#ifdef TIGER_USE_SLOW_OMP
      if(.not. fullyIntegralDirect .and. .not. integralDirect .and. .not. directFourInternal) then
        call destroy(omp_offsets_four_icount)
        call destroy(omp_offsets_four_nodecount)
      endif
#endif
      
  end subroutine SGGA_reboot

end module SGGA
