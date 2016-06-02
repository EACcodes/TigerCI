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
  !//  CI_UTILITIES_MOD - THIS MODULE CONTAINS A NUMBER OF UTILITY
  !//  FUNCTIONS WHICH ARE USED ONLY FOR CI CALCULATIONS
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************
  
module ci_utilities_mod
  
  use graph_var_mod
  use molecule_var_mod
  use utilities_mod
  use global_var_mod
  
  contains
  
function vertex(i,j)
  
      !// THIS FUNCTIONS RETURNS THE VERTEX NUMBER OF THE VERTEX WITH J ELECTRONS
      !// IN THE ITH ROW
  
      integer::vertex
      integer, intent(in)::i,j
  
      vertex = auxiliary(i) + j
  
end function vertex
  
function belong(i,j)
  
      !// THIS FUNCTION CHECKS TO SEE IF THE VERTEX HAVING J ELECTRONS IN
      !// ROW I IS IN THE GRAPH
  
      logical::belong
      integer, intent(in)::i,j
  
      if ((i<0).or.(j<0).or.(i > num_orbitals)) then
          belong =.false.
          return
      endif
      
      belong = ((j >= left_border(i)).and.(j <= right_border(i)).and.&
                (i <= num_orbitals).and.(i >= 0))
  
end function belong
  
function add0(i,j) 
  
      !// THIS FUNCTION CHECKS TO SEE IF VERTEX(I,J) HAS A 0 ARC HANGING OFF IT
  
      logical::add0
      integer, intent(in)::i,j
      
      if ((i<0).or.(j<0)) then
      add0 =.false.
      return
      endif
  
      add0 = (belong(i,j).and.(y0(vertex(i,j)) > 0))
  
end function add0
  
function add1(i,j) 
  
      !// THIS FUNCTION CHECKS TO SEE IF VERTEX(I,J) HAS A 1 ARC HANGING OFF IT
  
      logical::add1
      integer, intent(in)::i,j
  
      if ((i<0).or.(j<0)) then
      add1 =.false.
      return
      endif
      add1 = (belong(i,j).and.(y1(vertex(i,j)) >= 0))
  
end function add1
  
function add2(i,j) 
  
      !// THIS FUNCTION CHECKS TO SEE IF VERTEX(I,J) HAS A 2 ARC HANGING OFF IT
  
      logical::add2
      integer, intent(in)::i,j
  
      if ((i<0).or.(j<0)) then
      add2 =.false.
      return
      endif
      add2 = (belong(i,j).and.(y2(vertex(i,j)) >= 0))
  
end function add2
  
function num_paths(vertex)
  
      implicit none
  
      integer::allocatestatus      !// FOR DYNAMIC MEMORY
      integer::deallocatestatus    !// ALSO FOR DYNAMIC MEMORY
      integer::num_paths           !// THE PATHS WE ARE FINDING
      integer::level               !// THE CURRENT LEVEL
      integer::step_type           !// THE TYPE OF STEP WE TAKE
      integer::path_elecs          !// NUMBER OF ELECTRONS IN PATH
      integer, intent(in)::vertex  !// THE CURRENT TOP VERTEX
      integer, dimension(:), allocatable::step_array   !// KEEPS TRACK OF PATH WE MAKE

      allocate(step_array(opt_level:num_orbitals),stat=allocatestatus)
      call allocatecheck(allocatestatus,"step_arra")
  
      level = opt_level
      step_type = 2
      num_paths = 0
      path_elecs = vertex - auxiliary(opt_level)   
      step_array = 0
      
      repeater: do
  
      !// FIRST TEST TO EXIT
      if ((level==opt_level).and.(step_type < 0)) exit repeater
  
      !// TEST IF WE HAVE A NEW PATH
      if (level == num_orbitals) then
          num_paths = num_paths + 1
          level = level - 1
          path_elecs = path_elecs - step_array(level)
          step_type = step_array(level) - 1
          cycle repeater
          
      !// NOW TEST TO SEE IF WE HAVE RUN OUT OF STEPS
      elseif (step_type < 0) then
          level = level - 1
          path_elecs = path_elecs - step_array(level)
          step_type = step_array(level) - 1
          cycle repeater
  
      !// NOW, GO THROUGH THE STEP TYPES.  FIRST 2, THEN 1, THEN 0
      elseif (step_type == 2) then
          if (add2(level,path_elecs)) then  !// THE STEP IS ALLOWED
          step_array(level) = 2
          level = level + 1
          path_elecs = path_elecs + 2
          step_type = 2
          cycle repeater
          else   !// THE STEP IS NOT ALLOWED, TRY A 1 STEP FROM THE SAME LEVEL
          step_type = 1
          cycle repeater
          endif
      
      !// NOW A 1 STEP
      elseif (step_type == 1) then
          if (add1(level,path_elecs)) then  !// THE STEP IS ALLOWED
          step_array(level) = 1
          level = level + 1
          path_elecs = path_elecs + 1
          step_type = 2
          cycle repeater
          else   !// THE STEP IS NOT ALLOWED, TRY A 0 STEP FROM THE SAME LEVEL
          step_type = 0 
          cycle repeater
  
          endif
  
      !// NOW A 0 STEP
      elseif (step_type == 0) then
          if (add0(level,path_elecs)) then  !// THE STEP IS ALLOWED
          step_array(level) =  0
          level = level + 1
          step_type = 2
          path_elecs = path_elecs
          cycle repeater
          else   !// THE STEP IS NOT ALLOWED; STEP BY A LEVEL
          step_type = - 1
          cycle repeater
          endif
  
      endif
  
      enddo repeater
      
      deallocate(step_array, stat=deallocatestatus)
      call deallocatecheck(deallocatestatus,"step_arr")
      
end function num_paths 
  
function num_paths2(begin_lev,begin_elec,end_lev,end_elec)
      
      !// IN THIS FUNCTION WE COUNT THE NUMBER OF PATHS BETWEEN
      !// TWO VERTICES    
  
      implicit none
  
      integer::allocatestatus         !// FOR DYNAMIC MEMORY
      integer::deallocatestatus       !// ALSO FOR DYNAMIC MEMORY
      integer::num_paths2             !// THE PATHS WE ARE FINDING
      integer::level                  !// THE CURRENT LEVEL
      integer::step_type              !// THE TYPE OF STEP WE TAKE
      integer::path_elecs             !// NUMBER OF ELECTRONS IN PATH
      integer::num_needed             !// NUMBER OF ELECTRONS WE NEED
      integer::num_we_can_add         !// NUMBER OF ELECTRONS WE HAVE
      integer, intent(in)::begin_lev  !// BEGINNING LEVEL
      integer, intent(in)::end_lev    !// ENDING LEVEL
      integer, intent(in)::begin_elec !// BEGINNING NUMBER OF VERTICES
      integer, intent(in)::end_elec   !// ENDING NUMBER OF VERTICES
      integer, dimension(:), allocatable::step_array   !// KEEPS TRACK OF PATH WE MAKE

      allocate(step_array(begin_lev:end_lev),stat=allocatestatus)
      call allocatecheck(allocatestatus,"step_arra")
  
      level = begin_lev
      step_type = 2
      num_paths2 = 0
      path_elecs = begin_elec   
      step_array = 0
      
      if (end_elec < begin_elec) return
      if (level == end_lev.and.end_elec /= begin_elec) return
      
      repeater: do
  
      !// FIRST TEST TO EXIT
      if ((level==begin_lev).and.(step_type < 0)) exit repeater
  
      !// TEST IF WE HAVE A NEW PATH
      if (level == end_lev) then
          
          if (path_elecs==end_elec) num_paths2 = num_paths2 + 1
          level = level - 1
          path_elecs = path_elecs - step_array(level)
          step_type = step_array(level) - 1
          cycle repeater
          
      !// NOW TEST TO SEE IF WE HAVE RUN OUT OF STEPS
      elseif (step_type < 0) then
          level = level - 1
          path_elecs = path_elecs - step_array(level)
          step_type = step_array(level) - 1
          cycle repeater
  
      !// NOW, GO THROUGH THE STEP TYPES.  FIRST 2, THEN 1, THEN 0
      elseif (step_type == 2) then
          if (add2(level,path_elecs).and.&
              ((path_elecs+2) <= end_elec)) then  !// THE STEP IS ALLOWED
          step_array(level) = 2
          level = level + 1
          path_elecs = path_elecs + 2
          step_type = 2
          cycle repeater
          else   !// THE STEP IS NOT ALLOWED, TRY A 1 STEP FROM THE SAME LEVEL
          step_type = 1
          cycle repeater
          endif
      
      !// NOW A 1 STEP
      elseif (step_type == 1) then
          num_needed = end_elec-path_elecs-1
          num_we_can_add = 2*(end_lev - level - 1)
          if (add1(level,path_elecs).and.&
              ((path_elecs+1) <= end_elec).and.&
              (num_needed <= num_we_can_add)) then  !// THE STEP IS ALLOWED
          step_array(level) = 1
          level = level + 1
          path_elecs = path_elecs + 1
          step_type = 2
          cycle repeater
          else   !// THE STEP IS NOT ALLOWED, TRY A 0 STEP FROM THE SAME LEVEL
          step_type = 0 
          cycle repeater
  
          endif
  
      !// NOW A 0 STEP
      elseif (step_type == 0) then
          num_needed = end_elec-path_elecs
          num_we_can_add = 2*(end_lev - level - 1)
          if (add0(level,path_elecs).and.&
              (num_needed <= num_we_can_add)) then  !// THE STEP IS ALLOWED
          step_array(level) =  0
          level = level + 1
          step_type = 2
          path_elecs = path_elecs
          cycle repeater
          else   !// THE STEP IS NOT ALLOWED; STEP BY A LEVEL
          step_type = - 1
          cycle repeater
          endif
  
      endif
  
      enddo repeater
      
      deallocate(step_array, stat=deallocatestatus)
      call deallocatecheck(deallocatestatus,"step_arr")
      
end function num_paths2 
  
  
  
function maximum_singles(levels, number_electrons)
  
      !// THIS FUNCTION TELLS YOU THE MAXIMUM NUMBER A SINGLES
      !// THAT A PATH WITH A CERTAIN NUMBER OF LEVELS AND A 
      !// A CERTAIN NUMBER OF ELECTRONS CAN HAVE
      integer::maximum_singles
      integer, intent(in)::levels, number_electrons
  
      if (levels >= number_electrons) then
          maximum_singles = number_electrons 
      elseif (number_electrons <= 2*levels) then
          maximum_singles = levels - (number_electrons - levels)
      else
          write(ioError,450) number_electrons,levels 
          450 format (1x,"You have too many electrons and they won't fit!!",/,&
              1x,"Number electrons =",i4,/,&
              1x,"Number of levels =",i4)
                  stop
      endif
  
end function maximum_singles
  
  
function E_matrix(size,zopen_shells)
  
      integer,intent(in)::size
      integer,intent(in)::zopen_shells
      real(real8), dimension(size,size)::E_matrix
      
      integer::j  !// LOOP CONTROL VARIABLE
      integer::f2 !//  FSN(OPEN_SHELLS) - FSN(OPEN_SHELLS-2)
      
      f2 = fsn(zopen_shells-2)
      
      E_matrix = real(0.0,real8)
      
      do j = 1, f2
          E_matrix(j,j) = real(1.0,real8)
      enddo
      
      do j = f2+1,fsn(zopen_shells)
          E_matrix(j,j) = -real(1.0,real8)
      enddo
      
end function E_matrix
  
function E_matrix_as_vector(zopen_shells)
  
      integer,intent(in)::zopen_shells
      real(real8), dimension(fsn(zopen_shells))::E_matrix_as_vector
      
      integer::j
      integer::f2
      
      f2 = fsn(zopen_shells-2)
      
      do j = 1, f2
          E_matrix_as_vector(j) = real(1.0,real8)
      enddo
      
      do j = f2+1,fsn(zopen_shells)
          E_matrix_as_vector(j) = -real(1.0,real8)
      enddo
      
end function E_matrix_as_vector
  
  
!> \brief Given a set of orbital occupations, test if they match a reference
!> \param occupations   The internal orbital occupations for this configuration
!>
!> NOTE! If you pass part of a path (path%occupations) make sure you make the argument path%occupations(1:num_internal). A 0 element
!> will screw this up. 
logical function is_reference(occupations)
implicit none
integer, dimension(:) :: occupations
integer :: iRef
is_reference = .false.
do iRef = 1, num_ref
    if ( arraysAreEqual(occupations(1:num_internal), references(1:num_internal,iRef))) then 
        is_reference = .true.
        return
    end if
end do
end function is_reference

  
  
end module ci_utilities_mod
  
