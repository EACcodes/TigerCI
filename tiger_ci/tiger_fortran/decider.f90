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
!*****************************************************************
  !//  DECIDER - THIS ROUTINE JUST DECIDES WHERE YOU GO ONCE YOU HAVE
  !//  FINISHED A CYCLE OF THE TREE SEARCH
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS.
!*****************************************************************    

recursive subroutine decider(next_task, cho_data, the_path, loc_scr, twoVars, mod1vars, mod2vars, a, b)
  use decider_symbols
  use utilities_mod
  use tree_search_mod, only: add_constrained_level
  use diag_element_mod, only: build_head_element
  use two_seg_var_mod, only: twoModVars
  use new_tree_search_structs, only: orbital_path
  use two_seg_mod_3
  use three_and_four_seg_mod_1
  use three_and_four_seg_mod_2  
  use global_var_mod
  use cholesky_structs
  use locist_var_mod,only:locist_scratch
  use three_four_seg_var_mod
  use tree_search_mod, only: write_path, test_debug
  use natural_mod, only: one_part_purely_internal_compl, one_part_internal_ext_compl
  use blocked_locks_mod

  implicit none

  ! inputs
  type(orbital_path)::the_path
  integer,intent(in)::next_task

  ! optional inputs (extra data structures that the next_task might need
  type(cholesky_data), optional     :: cho_data
  type(locist_scratch), optional    :: loc_scr
  type(threefourmod1vars), optional :: mod1vars
  type(threefourmod2vars), optional :: mod2vars
  type(twoModVars), optional        :: twoVars
  type(blockedLockVectorType), optional :: a, b
  
  !// ALL WE DO IN THIS SUBROUTINE IS DECIDE WHERE TO GO NEXT
 
  ! It's a fancy function pointer
  
  select case(next_task)
      
      case(BLD_HEAD_ELEM)
          if (.not. present(loc_scr)) call missing_optional_var("build_head_element", "loc_scr")
          call build_head_element(the_path,loc_scr)
      
      case(UNSUPPORTED)
          write(*,*) "This method of using TIGER w/ 2 electron integrals is no longer supported."
          stop
          
      case(ADD_CNSTR_LVL)
          ! note that we don't check the optional arguments here 
          ! add_constrained_level takes all the data structures as optional arguments
         call add_constrained_level(the_path,cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b)        
          
      case(PATH_WRITER)
          call write_path(the_path,"decider   ")  
  
      case(PURELY_INT_COMP)
          call one_part_purely_internal_compl(a, the_path)

      case(INT_EXT_COMP)
          call one_part_internal_ext_compl(a, the_path)

      case(DEBUG)
          call test_debug(the_path)    

      case(LMO_INT_2_COMP)
          if (.not. present(twoVars)) call missing_optional_var("intern_two_seg_compl_vec_lmo", "twoVars")
          call intern_two_seg_compl_vec_lmo(the_path,twoVars)

      case(LMO_1_INT_COMP)
          if (.not. present(loc_scr)) call missing_optional_var("one_intern_seg_compl_vec_lmo", "loc_scr")
          if (.not. present(a)) call missing_optional_var("one_intern_seg_compl_vec_lmo", "a")
          if (.not. present(b)) call missing_optional_var("one_intern_seg_compl_vec_lmo", "b")
          call one_intern_seg_compl_vec_lmo(a%v, b%v, the_path,loc_scr)

      case(LMO_1_EXT_COMP)
          if (.not. present(loc_scr)) call missing_optional_var("extern_one_seg_compl_vec_lmo", "loc_scr")
          if (.not. present(mod1Vars)) call missing_optional_var("extern_one_seg_compl_vec_lmo", "mod1Vars")
          if (.not. present(a)) call missing_optional_var("extern_one_seg_compl_vec_lmo", "a")
          if (.not. present(b)) call missing_optional_var("extern_one_seg_compl_vec_lmo", "b")
          call extern_one_seg_compl_vec_lmo(a,b,the_path,loc_scr,mod1vars)

      case(LMO_3_INT_COMP)
          if (.not. present(mod2Vars)) call missing_optional_var("three_intern_seg_compl_vec_lmo", "mod2Vars")
          call three_intern_seg_compl_vec_lmo(the_path,mod2vars)

      case(LMO_2_INT_COMP)
          if (.not. present(mod2Vars)) call missing_optional_var("two_intern_seg_compl_vec_lmo", "mod2Vars")
          call two_intern_seg_compl_vec_lmo(the_path,mod2vars)

      case(CHO_IABC)
          if (.not. present(mod1Vars)) call missing_optional_var("iabc_sig_cho_big", "mod1Vars")
          call iabc_sig_cho_big(the_path,mod1vars)
          
      case default
          write(*,*) "Decider: Unknown value in decider ",next_task
          write(*,*) "So long and thanks for all the fish!"
          flush(6)
          stop

  end select
  
  contains
  
      !> \brief Produces an error message if a simple message about a missing (not present) variable
    !> Note that it doesn't do the checking ... it just produces the error message
    subroutine missing_optional_var(routine_name, variable_name)
        ! inputs
        character(len=*) :: routine_name, variable_name
        write(*,*) "FATAL DECIDER ERROR ... NEEDED OPTIONAL VARIABLE NOT PRESENT"
        write(*,*) "routine = " , routine_name
        write(*,*) "missing variable = " , variable_name
        call cause_floating_point_exception
        stop
end subroutine missing_optional_var
  
end subroutine decider  
