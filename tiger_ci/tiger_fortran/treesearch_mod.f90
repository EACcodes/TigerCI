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
  !//  TREESEARCH_MOD.F90 - THIS MODULE CONTAINS THE GENERIC TREE SEARCH 
  !//  ALGORITHM.  WE ARE GOING TO DEFINE A DERIVED DATA TYPE TO STORE ALL 
  !//  THE INFORMATION WE NEED, AND THEN WE WILL BUILD THE APPROPRIATE 
  !//  COMBINATIONS OF INTEGRALS IN SMALLER COMPANION ROUTINES SPECIFIC
  !//  TO THE INTERACTING CONFIGURATIONS.  
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS.
!*****************************************************************
  
module tree_search_mod 
  
  use global_var_mod
  use graph_var_mod
  use molecule_var_mod
  use utilities_mod
  use ci_utilities_mod
  use time_var_mod
  use decider_symbols
  use new_tree_search_structs                        
  use cholesky_structs
  use locist_var_mod,only:locist_scratch
  use three_four_seg_var_mod
  use two_seg_var_mod
  use blocked_locks_mod
  implicit none
  
  ! Explicit interface for the decider. Required because sticking the decider in a module
  ! creates (unsurprisingly) recursive module dependencies
  INTERFACE
    RECURSIVE SUBROUTINE DECIDER(NEXT_TASK, CHO_DATA,THE_PATH,LOC_SCR,TWOVARS,MOD1VARS,MOD2VARS,A,B)
        USE THREE_FOUR_SEG_VAR_MOD
        USE GLOBAL_VAR_MOD
        USE TWO_SEG_VAR_MOD
        USE LOCIST_VAR_MOD, ONLY : LOCIST_SCRATCH
        USE CHOLESKY_STRUCTS
        USE NEW_TREE_SEARCH_STRUCTS
        USE BLOCKED_LOCKS_MOD
        IMPLICIT NONE
        ! Required input
        INTEGER(KIND=8) :: NEXT_TASK
        ! Optional inputs
        TYPE (CHOLESKY_DATA), OPTIONAL :: CHO_DATA
        TYPE (ORBITAL_PATH), OPTIONAL  :: THE_PATH
        TYPE (LOCIST_SCRATCH), OPTIONAL  :: LOC_SCR
        TYPE (TWOMODVARS), OPTIONAL  :: TWOVARS
        TYPE (THREEFOURMOD1VARS), OPTIONAL  :: MOD1VARS
        TYPE (THREEFOURMOD2VARS), OPTIONAL  :: MOD2VARS
        TYPE(BLOCKEDLOCKVECTORTYPE) :: A,B
    END SUBROUTINE DECIDER
   END INTERFACE

      
contains
  
  
!*****************************************************************
recursive subroutine tree_search(the_path, top_level, top_occ, bottom_level, bottom_occ, next_task, &
                                 cho_data, loc_scr, twoVars, mod1vars, mod2vars, a, b )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
                         
  ! inputs
  type(orbital_path)::the_path
  integer::next_task                             !// INTEGER VARIABLE WHICH DETERMINES WHAT TO DO NEXT
  integer::top_level                             !// TOP LEVEL OF THE PATH
  integer::bottom_level                          !// BOTTOM LEVEL OF PATH
  integer::top_occ                               !// NUMBER OF ELECTRONS AT TOP LEVEL
  integer::bottom_occ                            !// NUMBER OF ELECTRONS AT BOTTOM LEVEL
  
  ! optional inputs
  type(cholesky_data), optional     :: cho_data
  type(locist_scratch), optional    :: loc_scr
  type(threefourmod1vars), optional :: mod1vars
  type(threefourmod2vars), optional :: mod2vars
  type(twoModVars), optional        :: twoVars
  type(blockedLockVectorType), optional :: a,b
  
  ! local variables
  real(real8)::zero
  integer::level                                 !// THE LEVEL WE ARE ON
  integer::current_vertex                        !// THE VERTEX WE ARE ON
  integer::ending_vertex                         !// THE VERTEX WE END ON
  integer::step_type                             !// THE TYPE OF STEP WE ARE TAKING
  integer::path_elecs                            !// THE NUMBER OF ELECTRONS IN THE PATH

  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  zero = real(0.0,real8)
  
  level = top_level
  step_type = 2
  path_elecs = top_occ
  current_vertex = vertex(top_level, top_occ)
  ending_vertex = vertex(bottom_level, bottom_occ)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// NOW, WE HAVE THE TREE SEARCH ALGORITHM.  WE ARE GOING TO GO FROM THE TOP 
  !// VERTEX TO THE BOTTOM VERTEX.  ALL OF THESE VARIABLES SHOULD BE INITIALIZED IN
  !// THE CALLING ROUTINE.    
  
  !// LETS GET OUT OF HERE IF THE TOP_LEVEL EQUALS THE BOTTOM_LEVEL
  if (top_level == bottom_level) then
  
      if (current_vertex /= ending_vertex) return
  
      call decider(next_task,cho_data,the_path,loc_scr,twoVars,mod1vars,mod2vars,a,b)
      return
  endif
  
  levels: do 
      
      if (level >= 0) current_vertex = vertex(level, path_elecs)
      
      !// FIRST CHECK TO SEE IF WE SHOULD EXIT
      if ((level <= top_level).and.(step_type == -1)) then
          
          exit levels
      
      !// NOW CHECK TO SEE IF WE HAVE A PATH TO PROCESS
      elseif (current_vertex==ending_vertex) then
          
          call decider(next_task, cho_data, the_path, loc_scr,twoVars,mod1vars,mod2vars,a,b)
  
          !// NOW STEP BACK A LEVEL.  WE HAVE TO DO THIS EXPLICITLY HERE
          !// OR WE JUST GET KICKED BACK HERE.
  
          step_type = the_path%occupations(level) - 1
          path_elecs = path_elecs - the_path%occupations(level)
  
          the_path%arc_weights(level) = 0
          if (the_path%occupations(level) == 1) then
          
              the_path%singles(the_path%num_singles) = 0
              the_path%num_singles = the_path%num_singles - 1
          
          endif
          level = level - 1
          cycle levels
          
      !// NOW TRY A TWO STEP
      elseif (step_type == 2) then
      
          !// TRY TO SEE IF WE CAN MAKE A 2 STEP
          if (add2(level,path_elecs).and.&
              ((path_elecs+2) <= bottom_occ).and.&
              (level < bottom_level)) then
  
              level = level + 1
  
              the_path%occupations(level) = step_type
              the_path%arc_weights(level) = abs(y2(current_vertex))
              path_elecs = path_elecs + 2
              step_type = 2
          else
              step_type = 1
              cycle levels
          endif
      
      !// NOW TRY A 1 STEP    
      elseif (step_type == 1) then
      
          !// TRY TO SEE IF WE CAN MAKE A 1 STEP
          if (add1(level,path_elecs).and.&
              ((path_elecs+1) <= bottom_occ).and.&
              (level < bottom_level).and.&
              (the_path%num_singles+1 <= open_shells)) then
  
              level = level + 1
              
              the_path%occupations(level) = step_type
              the_path%arc_weights(level) = abs(y1(current_vertex))
              path_elecs = path_elecs + 1
              the_path%num_singles = the_path%num_singles + 1
              the_path%singles(the_path%num_singles) = level
              step_type = 2
          else
              step_type = 0
              cycle levels
          endif
  
      !// NOW TRY A 0 STEP
      elseif (step_type == 0) then
          
          !// TRY TO SEE IF WE CAN MAKE A 0 STEP
          if (add0(level,path_elecs).and.&
             (level < bottom_level)) then
              
              level = level + 1
              the_path%occupations(level) = step_type
              step_type = 2
          else
              step_type = -1
              cycle levels
          endif
          
      !// UH OH, NOW WE NEED TO STEP BACK.  REMEMBER, THE ARCS ARE INDEXED
      !// ACCORDING TO THE LEVELS THEY POINT TO, AND "LEVEL" IS LABELLING
      !// THE ARC WE JUST TRIED TO JUMP OFF FROM.  BASICALLY, ALL I'M TRYING TO
      !// SAY HERE IS PAY ATTENTION AND KEEP THE INDEXING STRAIGHT.   
      elseif (step_type == -1) then
  
          step_type = the_path%occupations(level) - 1
          path_elecs = path_elecs - the_path%occupations(level)
          
          the_path%arc_weights(level) = 0
  
          if (the_path%occupations(level) == 1) then
              the_path%singles(the_path%num_singles) = 0
              the_path%num_singles = the_path%num_singles - 1
          endif
  
          level = level - 1
          cycle levels
      endif        
  
  enddo levels
  
  
end subroutine tree_search
  
  
!*****************************************************************
recursive subroutine broken_constrained_search(the_path, starting_level,starting_occ, &
                                               cho_data, loc_scr,twoVars,mod1vars,mod2vars,a,b)

  !// IN THIS ROUTINE WE BUILD PATHS GOING FROM THE TOP OF THE GRAPH TO
  !// THE END OF THE INTERNAL SPACE.  CONSTRAINTS CAN BE INCLUDED IN THE 
  !// SEARCH.  THIS IS DONE BY BREAKING THE SEARCH INTO PARTS.  AT EVERY
  !// CONSTRAINED LEVEL WE CALL AN AUXILIARY ROUTINE WHICH ADDS THE LEVEL.
  !// THEN WE CALL THE SEARCH ALGORITHM AGAIN.  
  
  ! inputs                      
  type(orbital_path)::the_path            !// THE DERIVED TYPE
  integer::starting_level, starting_occ   !// WHERE WE ARE STARTING TO SEARCH
  
  ! optional inputs
  type(cholesky_data), optional     :: cho_data
  type(locist_scratch), optional    :: loc_scr
  type(threefourmod1vars), optional :: mod1vars
  type(threefourmod2vars), optional :: mod2vars
  type(twoModVars), optional        :: twoVars
  type(blockedLockVectorType), optional :: a,b
  
  ! local variable
  integer::ending_level, ending_occ       !// WHERE WE END THE SEARCH
  integer::constraint_count               !// LABELS THE CONSTRAINT WE ARE ON
  integer::left_elec,right_elec           !// LEFTMOST AND RIGHTMOST ELECTRONS
  integer::next_task                      !// INTEGER TO KEEP TRACK OF THE NEXT THING TO DO
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  constraint_count = the_path%constraints(2,6)
  
  ending_level = the_path%constraints(2,constraint_count)-1
  next_task = the_path%constraints(1,6)
  
  left_elec = left_border(ending_level)
  
  right_elec = right_border(ending_level)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// EXIT IF WE ARE ALREADY AT THE INTERNAL/EXTERNAL SPACE BORDER.  TAKES CARE OF 
  !// THE SPECIAL CASE THAT ONE OF THE SEGMENTS ENDS AT THE INTERNAL/EXTERNAL SPACE
  !// BORDER 
  if (starting_level == num_internal) then
       call decider(next_task, cho_data, the_path,loc_scr,twoVars,mod1vars,mod2vars,a,b)
       return
  endif
    
  
  if (ending_level == num_internal) then
  
      !// ONCE WE COMPLETE THIS SEARCH WE ARE DONE
      do ending_occ = left_elec, right_elec
          call tree_search(the_path,starting_level,starting_occ,ending_level,ending_occ,next_task, &
                           cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b)
      enddo
  
  else
      
      !// HERE WE WILL NEED TO SEARCH AGAIN.  THUS, AFTER COMPLETING THE
      !// SEARCH WE CALL THE ROUTINE TO ADD THE CONSTRAINED SEGMENT
      do ending_occ = left_elec, right_elec
          call tree_search(the_path,starting_level,starting_occ,ending_level,ending_occ,ADD_CNSTR_LVL, &
                           cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b)
      enddo
  
  endif    
  
  
end subroutine broken_constrained_search
  
  
    
!*****************************************************************
recursive subroutine add_constrained_level(the_path, cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b)
  !// IN THIS ROUTINE WE TRY TO ADD A CONSTRAINED LEVEL.  THE LEVEL
  !// IS STORED IN THE CONSTRAINTS ARRAY.  ONCE WE HAVE ADDED THE 
  !// LEVEL, WE CALL THE TREE SEARCH ROUTINE AGAIN.  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! inputs
  type(orbital_path)::the_path
  
  ! optional inputs
  type(cholesky_data), optional     :: cho_data
  type(locist_scratch), optional    :: loc_scr
  type(threefourmod1vars), optional :: mod1vars
  type(threefourmod2vars), optional :: mod2vars
  type(twoModVars), optional        :: twoVars
  type(blockedLockVectorType), optional :: a,b
  
  ! local variables
  integer::constrained_level      !// THE LEVEL WHERE THE CONSTRAINT IS
  integer::lambda_step_type       !// THE OCCUPATION OF THE CONSTRAINED ORBITAL IN LAMBDA
  integer::mu_step_type           !// OCCUPATION OF CONSTRAINED ORBITAL IN MU
  integer::constraint_count       !// LABLES THE CONSTRAINT
  integer::distance               !// SEPARATION BETWEEN MU AND LAMBDA PATHS
  integer::path_elecs             !// OCCUPATION IN LAMBDA PATH UP TO CONSTRAINED ORBITAL
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  constraint_count = the_path%constraints(2,6)
  
  constrained_level = the_path%constraints(2,constraint_count)
  lambda_step_type = the_path%constraints(1,constraint_count)
  mu_step_type = the_path%constraints(3,constraint_count)
  distance = sum(the_path%constraints(3,1:constraint_count-1))-&
             sum(the_path%constraints(1,1:constraint_count-1))
             
  path_elecs = sum(the_path%occupations(0:constrained_level-1))
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// NOW LETS TRY TO ADD THE STEP.  AFTER WE HAVE SUCCESFULLY ADDED THE STEP
  !// WE NEED TO CHANGE EVERYTHING BACK TO THE WAY IT WAS BEFORE EXITING THIS ROUTINE
  
  !// FIRST CHECK THE MU STEP.  WE DON'T DO ANYTHING HERE EXCEPT CHECK TO SEE IF
  !// WE CAN ADD IT.  IT IS BELIEVED THIS WILL SAVE WORK ON AVERAGE.
  
  if (mu_step_type == 2) then
      
      !// EXIT IF WE CAN'T ADD THE LEVEL
      if (.not.add2(constrained_level-1,path_elecs+distance)) return
      
  elseif (mu_step_type == 1) then
  
      !// EXIT IF WE CAN'T ADD THE LEVEL
      if (.not.add1(constrained_level-1,path_elecs+distance)) return    
      
  else    
      
      !// EXIT IF WE CAN'T ADD THE LEVEL
      if (.not.add0(constrained_level-1,path_elecs+distance)) return    
      
  endif    
  
  !// NOW DO THE LAMBDA STEP.  
  if (lambda_step_type == 2) then
  
      !// EXIT IF WE CAN'T ADD THE LEVEL
      if (.not.add2(constrained_level-1,path_elecs)) return
      
      !// ADD THE STEP
      the_path%occupations(constrained_level)= 2
      the_path%arc_weights(constrained_level)=abs(y2(vertex(constrained_level-1,path_elecs)))
      path_elecs = path_elecs + 2
      
      !// UPDATE THE COSNTRAINTS
      the_path%constraints(2,6) = constraint_count+1    
      
      if (constrained_level+1 == the_path%constraints(2,constraint_count+1).and.&
          constrained_level /= num_internal) then
          call add_constrained_level(the_path,cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b) 
      else
          call broken_constrained_search(the_path,constrained_level,path_elecs, &
                                         cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b)
      endif        
      
      !// CHANGE STUFF BACK
      the_path%constraints(2,6) = constraint_count
      
  elseif (lambda_step_type == 1)  then
  
      !// EXIT IF WE CAN'T ADD THE LEVEL
      if (.not.add1(constrained_level-1,path_elecs).or.&
          the_path%num_singles + 1 > open_shells) return
      
      !// ADD THE STEP
      the_path%occupations(constrained_level)=1
      the_path%arc_weights(constrained_level)=abs(y1(vertex(constrained_level-1,path_elecs)))
      the_path%num_singles = the_path%num_singles + 1
      the_path%singles(the_path%num_singles) = constrained_level
      path_elecs = path_elecs + 1
      
      !// UPDATE THE CONSTRAINTS
      the_path%constraints(2,6) = constraint_count+1    
      
      if (constrained_level+1 == the_path%constraints(2,constraint_count+1).and.&
          constrained_level /= num_internal) then
          call add_constrained_level(the_path,cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b) 
      else
          call broken_constrained_search(the_path,constrained_level,path_elecs, &
                                         cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b)
      endif        
      
      !// CHANGE STUFF BACK
      the_path%constraints(2,6) = constraint_count
      !the_path%num_singles = the_path%num_singles -1
      the_path%singles(the_path%num_singles) = 0
      the_path%num_singles = the_path%num_singles -1
  
  elseif (lambda_step_type == 0) then
  
      !// EXIT IF WE CAN'T ADD THE LEVEL
      if (.not.add0(constrained_level-1,path_elecs)) return
      
      !// ADD THE STEP
      the_path%occupations(constrained_level)=0
      the_path%arc_weights(constrained_level)=0
      
      !// UPDATE THE COSNTRAINTS
      the_path%constraints(2,6) = constraint_count+1    
      
      
      if (constrained_level+1 == the_path%constraints(2,constraint_count+1).and.&
          constrained_level /= num_internal) then
          call add_constrained_level(the_path,cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b) 
      else
          call broken_constrained_search(the_path,constrained_level,path_elecs, &
                                         cho_data,loc_scr,twoVars,mod1vars,mod2vars,a,b)
      endif            
  
      !// CHANGE STUFF BACK
      the_path%constraints(2,6) = constraint_count             
      
  endif    
  
end subroutine add_constrained_level
  

!  
!*****************************************************************
subroutine test_debug(the_path)
 
  type(orbital_path)::the_path

  write(ioOutput,*) "Well, I found the path.  I don't know why you can't!// "
  call write_path(the_path,"test      ")
  stop
  
end subroutine test_debug

subroutine copyPath(this, that)
  type(orbital_path)::this,that

  ! first the scalar values
  that%num_singles = this%num_singles
  that%weight = this%weight
  that%level1 = this%level1
  that%level2 = this%level2
  that%level3 = this%level3
  that%level4 = this%level4
  
  that%i_segment = this%i_segment
  that%j_segment = this%j_segment
  that%k_segment = this%k_segment
  that%l_segment = this%l_segment
  
  that%loop_type = this%loop_type
  that%rt_loop_weight = this%rt_loop_weight
  that%start_elec = this%start_elec
  
  that%element1 = this%element1
  that%element2 = this%element2
  that%element3 = this%element3
 
  ! now all array types
  if(allocated(this%occupations)) then
    call copyArrayInt(this%occupations,that%occupations)
  endif
  if(allocated(this%arc_weights)) then
    call copyArrayInt(this%arc_weights,that%arc_weights)
  endif
  if(allocated(this%singles)) then
    call copyArrayInt(this%singles,that%singles)
  endif
  if(allocated(this%constraints)) then
    call copyMatrixInt(this%constraints,that%constraints)
  endif
  if(allocated(this%encountered)) then
    call copyArrayBool(this%encountered,that%encountered)
  endif
  if(allocated(this%integrals)) then
    call copyArray(this%integrals,that%integrals)
  endif

end subroutine copyPath

end module tree_search_mod
   
