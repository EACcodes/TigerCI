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
  !//  GRAPH_VAR_MOD - THIS MODULE STORE INFORMATION ABOUT THE GRAPH.
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************
  
module graph_var_mod

  use global_var_mod
  use utilities_mod
  
  implicit none
  
  integer, dimension(:), allocatable::left_border        !// LEFT BORDER OF GRAPH
  integer, dimension(:), allocatable::right_border       !// RIGHT BORDER OF GRAPH  
  integer, dimension(:), allocatable::auxiliary          !// AUXILIARY ARRAY OF DUCH AND KARWOWSKI.  AUXILIARY(I) + J RETURNS THE 
                                                         !// NUMBER OF THE JTH VERTEX ON ROW I.
  integer, dimension(:), allocatable::y0                 !// THE VERTEX WEIGHT ARRAY
  integer, dimension(:), allocatable::y1                 !// THE FIRST ARC WEIGHT ARRAY
  integer, dimension(:), allocatable::y2                 !// THE SECOND ARC WEIGHT ARRAY
  integer, dimension(:), allocatable::fsn                !// NUMBER OF PSIN FUNCTIONS FOR NUMBER OF OPEN SHELLS
                                                         !// AND THE GIVEN SPIN MULTIPLICITY
  integer, dimension(:), allocatable::internal_index_vector0   !// FOR THROWING OUT INTERNAL CONFIGURATIONS; VALENCE
  integer, dimension(:), allocatable::internal_index_vector1   !// FOR THROWING OUT INTERNAL CONFIGURATIONS; N-1
  integer, dimension(:), allocatable::internal_index_vector2   !// FOR THROWING OUT INTERNAL CONFIGURATIONS; N-2, TWO VIRTUALS
  integer, dimension(:), allocatable::internal_index_vector3   !// FOR THROWING OUT INTERNAL CONFIGURATIONS; N-2, ONE VIRTUAL
  !integer, dimension(:), allocatable::index_vector_pao3        !// PAO INDEX VECTOR FOR N-2 ONE DOUBLY OCC. VIRTUAL 
  !integer, dimension(:), allocatable::index_vector_pao2        !// PAO INDEX VECTOR FOR N-2 TWO SINGLY OCC. VIRTUALS 
  !integer, dimension(:), allocatable::index_vector_pao1        !// PAO INDEX VECTOR FOR N-1 ONE SINGLY OCC. VIRTUAL
  
  
  integer, dimension(:), allocatable::v_singles      !// NUMBER OF INTERNAL SPACE SINGLES FOR VALENCE STATES
  integer, dimension(:), allocatable::nm1_singles    !// NUMBER OF INTERNAL SPACE SINGLES FOR N-1 STATES
  integer, dimension(:), allocatable::nm2_singles    !// NUMBER OF INTERNAL SPACE SINGLES FOR N-2 STATES
  
  integer, dimension(:), allocatable::fsn3               !//  
  integer,parameter::exilevel = 2                        !// EXCITATION LEVEL
  integer::open_shells                                   !// MAXIMUM NUMBER OF OPEN SHELLS
  integer::num_vert                                      !// NUMBER OF VERTICES
  integer::opt_level                                     !// OPTIMAL LEVEL FOR CUTTING GRAPH INTO TWO PARTS
  integer::total_csfs                                    !// TOTAL NUMBER OF CSFS 
  
  contains
  
  subroutine delete_graph_variables()
      implicit none
      integer :: ierr
      
      deallocate(left_border, right_border, auxiliary, y0, y1, y2, fsn, fsn3, internal_index_vector0, &
      internal_index_vector1, internal_index_vector2, internal_index_vector3, v_singles, nm1_singles,&
      nm2_singles, stat=ierr)
      call deallocatecheck(ierr, "delete_graph_variables")
      
      open_shells = 0
      num_vert = 0 
      opt_level = 0 
      total_csfs = 0 
  end subroutine
    
end module graph_var_mod 
