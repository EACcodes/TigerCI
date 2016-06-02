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
  !//  SPIN_VAR_MOD - A MODULE WHICH WE ARE GOING TO USE TO PASS AROUND
  !//  VARIABLES RELATED TO THE SPIN INTEGRALS.  
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  ! WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************
  
module spin_var_mod
  use utilities_mod
  use global_var_mod
  implicit none
  
  integer, dimension(:,:), allocatable::branch             !// THE BRANCHING DIAGRAM
  integer, dimension(:,:), allocatable::sft                !// THE SFT GRAPH
  integer, dimension(:,:), allocatable::spin_matrix        !// MATRIX WHICH RECORDS SPIN GENEOLGOY OF ALL SPIN FUNCTIONS
  
  real(real8), dimension(:,:,:), allocatable::cycles        !// STORES ALL THE CYCLES
  real(real8), dimension(:,:), allocatable::transpositions  !// STORES THE TRANSPOSITIONS
  
  
  contains
  
  subroutine clean_spin_var_mod()
      implicit none
      call try_deallocate_int_2D(branch, "branch")
      call try_deallocate_int_2D(sft, "SFT")
      call try_deallocate_int_2D(spin_matrix, "spin_matrix")
      call try_deallocate_real_2D(transpositions, "transpositions")
      call try_deallocate_real_3D(cycles, "cycles")
  end subroutine
      
  
end module spin_var_mod
  
  
