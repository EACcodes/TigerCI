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
  !//  TIME_VAR_MOD - THIS MODULE IS USED TO STORE THE TIMING VARIABLES
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  ! WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************
  
module time_var_mod
  
  use global_var_mod 
  
  implicit none
  
  
  real(real8)::sort_time               !// TIME SPENT SORTING INTEGRALS
  real(real8)::graph_time              !// TIME SPENT MAKING THE SPATIAL GRAPH AND INDEX VECTOR
  real(real8)::spin_time               !// TIME SPENT MAKING CYCLE BRANCHING DIAGRAM, SFT GRAPH, CYCLE MATRICES
  								     !// AND TRANSPOSITION MATRICES
  real(real8)::search_time             !// TIME SPENT SEARCHING								     
  real(real8)::diag_tot_time           !// TOTAL TIME SPENT IN ALL DIAGONALIZATION ROUTINES
  real(real8)::diag_element_time       !// TIME SPENT MAKING DIAGONAL ELEMENTS OF H.  
  real(real8)::diagtimesc_time         !// TIME SPENT MUTLIPLYING DIAGONAL ELEMENTS BY CI VECTOR
  real(real8)::lp26_time               !// TIME SPENT TAKING CARE OF CONTRIBUTIONS FROM {26} LOOPS
  real(real8)::two_seg_time            !// TIME SPENT ON TWO SEGMENT LOOPS
  real(real8)::three_and_four_seg_time !// TIME SPENT ON THREE AND FOUR SEGEMNT LOOPS
  
  
  real(real8)::ijkl_time               !// TIME SPENT TREATING VARIOUS INTEGRALS
  real(real8)::aijk_time
  real(real8)::abij_time               !// THIS ALSO INCLUDES THE (IA|JA)
  real(real8)::ikjk_time
  real(real8)::abcd_time               !// THIS ALSO INCLUDES (AC|BC)
  real(real8)::iabc_time          
  real(real8)::ijaj_time      
  
  real(real8)::rc_time
  real(real8)::coulomb_time
  real(real8)::exchange_time
  real(real8)::qbrc_time
  real(real8)::aqbrc_time
  
end module time_var_mod
