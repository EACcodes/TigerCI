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
module sgga_helpers

  use new_tree_search_structs
  use molecule_var_mod

  implicit none

  private
  
  public :: get_numbers_of_singles,get_numbers_of_singles2
  
  contains
       
  subroutine get_numbers_of_singles(s1,s2,s3,s4,lambda_path)

    !// IN THIS ROUTINE WE COMPUTE THE NUMBERS OF SINGLES IN BETWEEN LOOP SEGMENTS.
    !// FOR EXAMPLE, S1 IS THE NUMBER OF SINGLES FROM THE HEAD OF THE GRAPH UP TO
    !// AND INCLUDING THE FIRST LOOP SEGMENT.  S2 IS THE NUMBER OF SINGLES FROM THE
    !// FIRST LOOP SEGMENT UP TO AND INCLUDING THE SECOND LOOP SEGMENT.  WE COUNT UP 
    !// SINGLES UP UNTIL THE INTERNAL-EXTERNAL SPACE BORDER AT WHICH POINT WE JUST
    !// STORE THE LEFT OVER NUMBER OF SINGLES IN THE LOWEST REMAINING VARIABLE OF
    !// S1,S2,S3,S4.
    
    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(orbital_path),intent(in)::lambda_path

    integer,intent(out)::s1,s2,s3,s4            !// NUMBER OF SINGLES IN EACH SEGMENT OF THE LOOP.

    call get_numbers_of_singles2(s1,s2,s3,s4,lambda_path%occupations,lambda_path%constraints)

  end subroutine get_numbers_of_singles
  
  subroutine get_numbers_of_singles2(s1,s2,s3,s4,occupations,constraints)
  
  !// IN THIS ROUTINE WE COMPUTE THE NUMBERS OF SINGLES IN BETWEEN LOOP SEGMENTS.
  !// FOR EXAMPLE, S1 IS THE NUMBER OF SINGLES FROM THE HEAD OF THE GRAPH UP TO
  !// AND INCLUDING THE FIRST LOOP SEGMENT.  S2 IS THE NUMBER OF SINGLES FROM THE
  !// FIRST LOOP SEGMENT UP TO AND INCLUDING THE SECOND LOOP SEGMENT.  WE COUNT UP 
  !// SINGLES UP UNTIL THE INTERNAL-EXTERNAL SPACE BORDER AT WHICH POINT WE JUST
  !// STORE THE LEFT OVER NUMBER OF SINGLES IN THE LOWEST REMAINING VARIABLE OF
  !// S1,S2,S3,S4.  
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! USING ALLOCATABLE HERE IS A FIX FOR INTEL COMPOSER V2013, REMOVE ONCE FIXED!!!
  integer,dimension(:),allocatable,intent(in)::occupations
  integer,dimension(:,:),allocatable,intent(in)::constraints
  integer,intent(out)::s1,s2,s3,s4            !// NUMBER OF SINGLES IN EACH SEGMENT OF THE LOOP.
                                  !//   I.E. S1 IS THE NUMBER OF SINGLES FROM THE HEAD
                                  !//   OF THE GRAPH UP TO AND INCLUDING THE FIRST LOOP
                                  !//   SEGMENT
  integer::singles_count          !// TEMPORARY STORAGE OF THE NUMBER OF SINGLES                                 
  integer::constraint_count       !// KEEPS TRACK OF THE LOOP LEVELS
  integer::internal_orbital       !// LOOP VARIABLE FOR LOOPING OVER INTERNAL ORBITALS
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  singles_count = 0
  constraint_count = 1
  s1=0
  s2=0
  s3=0
  s4=0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// LOOP OVER THE ORBITALS AND COUNT UP THE NUMBER OF SINGLES
  do internal_orbital = 1, num_internal
  
      if (occupations(internal_orbital)==1) then
          singles_count = singles_count + 1
      endif
      
      if (constraints(2,constraint_count)==internal_orbital) then
          if (constraint_count == 1) then
              s1 = singles_count
              singles_count = 0
              constraint_count = constraint_count + 1
          elseif (constraint_count == 2) then
              s2 = singles_count
              singles_count = 0
              constraint_count = constraint_count + 1
          elseif (constraint_count == 3) then
              s3 = singles_count
              singles_count = 0
              constraint_count = constraint_count + 1
          elseif (constraint_count == 4) then
              s4 = singles_count
              singles_count = 0
              constraint_count = constraint_count + 1
          endif
      endif
  
  enddo
  
  !// PUT THE NUMBER OF SINGLES GOING FROM THE LAST REQUESTED SEGMENT
  !// TO THE INTERAL-EXTERNAL SPACE BORDER INTO THE APPROPRIATE SPOT
  if (constraint_count == 1) then
      s1 = singles_count
  elseif (constraint_count == 2) then
      s2 = singles_count
  elseif (constraint_count == 3) then
      s3 = singles_count
  elseif (constraint_count == 4) then
      s4 = singles_count
  endif
  
end subroutine get_numbers_of_singles2

end module sgga_helpers
