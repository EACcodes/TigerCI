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
  !//  SCEPMAKER - THIS FILE CONTAINS ROUTINES FOR SCEPPING AND
  !//  UNSCEPPING COEFICIENT MATRICES.  BASICALLY, THESE MATRICES
  !//  ARE REDUNDANT MATRICES WHERE A COMBINED INDEX "[AB]" IS 
  !//  DUPLICATED AS THE "AB" AND "BA" ELEMENT OF A MATRIX.
  !// 
  !//  WRITTEN BY DEREK WALTER, 2001
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS.
  !//  CLEANED UP BY JOHANNES M DIETERICH, 2013
  !//
  !//  MODIFIED BY CAROLINE M KRAUTER, 2014
  !//  OLD CODE ACCESSED ADDRESSES WHICH DID NOT EXIST
  !//  SPLITTING SUBROUTINES INTO SEPARATE ONES FOR 
  !//  EACH VALUE OF MODE

!*****************************************************************
  
module scepmaker_mod
  
  use global_var_mod
  use blocked_locks_mod
  !use sigmavec_var_mod
  implicit none
  
  contains
  
  
!*****************************************************************
subroutine scepper(scep_matrix, straight_vector, offset, dimension, mode)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
  
  !// Does not have to be split because always used correctly.
  !// CMK

  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
  integer,intent(in)::mode               !// FOR MAKING SYMMETRIC OR SKEW MATRICES
  integer,intent(in)::offset             !// OFFSET IN THE VECTOR
  real(real8),intent(inout)::scep_matrix(:,:)
  real(real8),intent(in)::straight_vector(:)
  
  integer::i,j,count                     !// FOR LOOPING OVER ELEMENTS
  real(real8)::matrix_element
  real(real8),parameter::zero= real(0.0, real8)
  
  count = offset-1
  
  if (mode == 1) then
  
      do i = 1, dimension
          do j = 1,i-1
      
              count = count + 1
              matrix_element = straight_vector(count)
          
              scep_matrix(j,i) = matrix_element
              scep_matrix(i,j) = matrix_element
          
          enddo
      
          scep_matrix(i,i) = zero
      
      enddo    
      
  elseif (mode == -1) then
  
      do i = 1, dimension
          do j = 1,i-1
      
              count = count + 1
              matrix_element = straight_vector(count)
          
              scep_matrix(j,i) =  matrix_element
              scep_matrix(i,j) = -matrix_element
          
          enddo
      
          scep_matrix(i,i) = zero
      
      enddo    
      
  endif    
  
end subroutine scepper
  
!*****************************************************************

subroutine scepper_diag_p1(scep_matrix, straight_vector, off_straight, dimension, diag_vector, off_diag)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
  !// THIS ROUTINE IS IDENTICAL TO SCEPPER WITH THE EXCEPTION THAT 
  !// THE DIAGONAL ELEMENTS ARE ALSO INCLUDED FOR SINGLET COUPLED
  !// SPIN FUNCTIONS.
  
  !// Uses diag_vector. Only call if aa_address valid!
  !// mode=+1 (singlet case?)
  !// CMK

  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR MAKING SYMMETRIC OR SKEW MATRICES
  integer,intent(in)::off_straight,off_diag
  real(real8),intent(inout)::scep_matrix(:,:)
  real(real8),intent(in)::diag_vector(:)
  real(real8),intent(in)::straight_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  real(real8)::matrix_element
  real(real8),parameter::zero=real(0.0, real8), sqrt2=real(sqrt(real(2.0,real8)), real8)
  
  count = off_straight-1

  
  do i = 1, dimension
  
      do j = 1,i-1
  
          count = count + 1
          matrix_element = straight_vector(count)
      
          scep_matrix(j,i) = matrix_element
          scep_matrix(i,j) = matrix_element
      
      enddo
  
      scep_matrix(i,i) = diag_vector(i+off_diag-1) * sqrt2
  
  enddo    
      
end subroutine scepper_diag_p1

!*****************************************************************

subroutine scepper_diag_m1(scep_matrix, straight_vector, off_straight, dimension)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
  !// THIS ROUTINE IS IDENTICAL TO SCEPPER WITH THE EXCEPTION THAT 
  !// THE DIAGONAL ELEMENTS ARE ALSO INCLUDED FOR SINGLET COUPLED
  !// SPIN FUNCTIONS.
  
  !// Does not use diag_vector!
  !// mode=-1 (triplet case?)
  !// Is the same as scepper for mode=-1. Replace?
  !// CMK

  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR MAKING SYMMETRIC OR SKEW MATRICES
  integer,intent(in)::off_straight !,off_diag
  real(real8),intent(inout)::scep_matrix(:,:)
!  real(real8),intent(in)::diag_vector(:)
  real(real8),intent(in)::straight_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  real(real8)::matrix_element
  real(real8),parameter::zero=real(0.0, real8), sqrt2=real(sqrt(real(2.0,real8)), real8)
  
  count = off_straight-1

  do i = 1, dimension
      do j = 1,i-1
  
          count = count + 1
          matrix_element = straight_vector(count)
      
          scep_matrix(j,i) =  matrix_element
          scep_matrix(i,j) = -matrix_element
      
      enddo
  
      scep_matrix(i,i) = zero
  
  enddo
      
end subroutine scepper_diag_m1

!*****************************************************************

subroutine scepper_diag_m2(scep_matrix, straight_vector, off_straight, dimension)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
  !// THIS ROUTINE IS IDENTICAL TO SCEPPER WITH THE EXCEPTION THAT 
  !// THE DIAGONAL ELEMENTS ARE ALSO INCLUDED FOR SINGLET COUPLED
  !// SPIN FUNCTIONS.

  !// Does not use diag_vector!
  !// mode=-2 
  !// CMK
  
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR MAKING SYMMETRIC OR SKEW MATRICES
  integer,intent(in)::off_straight !,off_diag
  real(real8),intent(inout)::scep_matrix(:,:)
!  real(real8),intent(in)::diag_vector(:)
  real(real8),intent(in)::straight_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  real(real8)::matrix_element
  real(real8),parameter::zero=real(0.0, real8), sqrt2=real(sqrt(real(2.0,real8)), real8)
  
  count = off_straight-1

  do i = 1, dimension
      do j = 1,i-1
  
          count = count + 1
          matrix_element = straight_vector(count)
      
          scep_matrix(j,i) = -matrix_element
          scep_matrix(i,j) =  matrix_element
      
      enddo
  
      scep_matrix(i,i) = zero
  
  enddo
      
end subroutine scepper_diag_m2

!*****************************************************************

subroutine scepper_diag(scep_matrix, straight_vector, off_straight, dimension, mode, diag_vector, off_diag)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
  !// THIS ROUTINE IS IDENTICAL TO SCEPPER WITH THE EXCEPTION THAT 
  !// THE DIAGONAL ELEMENTS ARE ALSO INCLUDED FOR SINGLET COUPLED
  !// SPIN FUNCTIONS.
  
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
  integer,intent(in)::mode               !// FOR MAKING SYMMETRIC OR SKEW MATRICES
  integer,intent(in)::off_straight,off_diag
  real(real8),intent(inout)::scep_matrix(:,:)
  real(real8),intent(in)::diag_vector(:)
  real(real8),intent(in)::straight_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  real(real8)::matrix_element
  real(real8),parameter::zero=real(0.0, real8), sqrt2=real(sqrt(real(2.0,real8)), real8)
  
  count = off_straight-1

  
  if (mode == 1) then
  
      do i = 1, dimension
  
          do j = 1,i-1
      
              count = count + 1
              matrix_element = straight_vector(count)
          
              scep_matrix(j,i) = matrix_element
              scep_matrix(i,j) = matrix_element
          
          enddo
      
          scep_matrix(i,i) = diag_vector(i+off_diag-1) * sqrt2
      
      enddo    
      
  elseif (mode == -1) then
  
      do i = 1, dimension
          do j = 1,i-1
      
              count = count + 1
              matrix_element = straight_vector(count)
          
              scep_matrix(j,i) =  matrix_element
              scep_matrix(i,j) = -matrix_element
          
          enddo
      
          scep_matrix(i,i) = zero
      
      enddo
      
  elseif (mode == -2) then
  
      do i = 1, dimension
          do j = 1,i-1
      
              count = count + 1
              matrix_element = straight_vector(count)
          
              scep_matrix(j,i) = -matrix_element
              scep_matrix(i,j) =  matrix_element
          
          enddo
      
          scep_matrix(i,i) = zero
      
      enddo
      
  endif    
  
end subroutine scepper_diag

!*****************************************************************

subroutine unscepper(scep_matrix, straight_vector, offset, dimension, mode)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
 
  !// Does not have to be split because always used correctly.
  !// CMK
 
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
  integer,intent(in)::offset
  real(real8),intent(in)::scep_matrix(:,:)
  real(real8),intent(inout)::straight_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  
  
  count = offset
  if (mode == 1) then
  
      do i = 1, dimension
          do j = 1,i-1
      
              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)+&
                                                                scep_matrix(i,j)
              count = count + 1
          
          enddo
      enddo    
      
  elseif (mode == -1) then
  
      do i = 1, dimension
          do j = 1,i-1
      
              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)-&
                                                                scep_matrix(i,j)
              count = count + 1
          
          enddo
      enddo    
  
  endif    
  
end subroutine unscepper

!*****************************************************************

!subroutine unscepper_sigma_atomic(scep_matrix, straight_vector, offset, dimension, mode)
!  
!  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
!  
!  implicit none
!  
!  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
!  integer,intent(in)::offset
!  real(real8),intent(in)::scep_matrix(:,:)
!  real(real8),intent(inout)::straight_vector(:)
!
!  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
!  integer::startp,endp
!
!  startp = offset
!  endp = startp + (dimension*dimension-dimension)/2
!
!#ifdef TIGER_USE_OMP
!  call setSigmaLocks(startp,endp)
!#endif
!  
!  if (mode == 1) then
!  
!      count = offset-1
!      do i = 1, dimension
!          do j = 1,i-1
!      
!              count = count + 1
!#ifdef DEBUG_TIGER
!              if(count > endp) then
!                write(*,*) "scepper: ",count,endp
!                stop
!              endif
!#endif
!              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)+scep_matrix(i,j)
!          
!          enddo
!      enddo    
!      
!  elseif (mode == -1) then
!  
!      count = offset-1
!      do i = 1, dimension
!          do j = 1,i-1
!      
!              count = count + 1
!#ifdef DEBUG_TIGER
!              if(count > endp) then
!                write(*,*) "scepper: ",count,endp
!                stop
!              endif
!#endif
!              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)-scep_matrix(i,j)
!          
!          enddo
!      enddo    
!  
!  endif    
!#ifdef TIGER_USE_OMP
!  call unsetSigmaLocks(startp,endp)
!#endif
!  
!end subroutine unscepper_sigma_atomic
  
!*****************************************************************

subroutine scep_half(scep_matrix, straight_vector, dimension, dc2, mode)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP

  !// Never used. 
  !// CMK
  
  implicit none
  
  integer::dimension          !// SIZE OF THE SCEP MATRIX
  integer::dc2                !// LENGTH OF STRAIGHT VECTOR = DIMENSION CHOSE 2
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  integer::mode               !// FOR MAKING SYMMETRIC OR SKEW MATRICES
  
  
  real(real8)::scep_matrix(:,:)
  real(real8)::straight_vector(dc2)
  real(real8)::matrix_element
  real(real8),parameter::zero = real(0.0, real8)

  count = 0
  
  
  if (mode == 1) then
  
      !// UPPER HALF
  
      do i = 1, dimension
          do j = 1,i-1
      
              count = count + 1
              matrix_element = straight_vector(count)
              scep_matrix(j,i) = matrix_element
          
          enddo
      
          scep_matrix(i,i) = zero
      
      enddo    
      
  elseif (mode == -1) then
  
      !// LOWER HALF
  
      do i = 1, dimension
          do j = 1,i-1
      
              count = count + 1
              matrix_element = straight_vector(count)
              scep_matrix(i,j) = matrix_element
          
          enddo
      
          scep_matrix(i,i) = zero
      
      enddo    
      
  endif    
  
end subroutine scep_half
  
!*****************************************************************

subroutine unscepper_diag_p1(scep_matrix, straight_vector, off_straight, dimension, diag_vector, off_diag)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP

  !// Uses diag_vector. Only call if aa_address valid!!!!!!!
  !// mode=+1 (singlet case?)
  !// CMK
 
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
  integer,intent(in)::off_diag,off_straight
  real(real8),intent(in)::scep_matrix(:,:)
  real(real8),intent(inout)::straight_vector(:)
  real(real8),intent(inout)::diag_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  real(real8),parameter::sqrt2=real(sqrt(real(2.0,real8)), real8)
  
  count = off_straight-1
  do i = 1, dimension
      do j = 1,i-1
                  
          count = count + 1
          straight_vector(count) = straight_vector(count) + scep_matrix(j,i)+&
                                                            scep_matrix(i,j)
      enddo
      diag_vector(i+off_diag-1)  = diag_vector(i+off_diag-1) + scep_matrix(i,i)*sqrt2
  enddo

end subroutine unscepper_diag_p1

!*****************************************************************

subroutine unscepper_diag_m12(scep_matrix, straight_vector, off_straight, dimension)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP

  !// Does not use diag_vector!
  !// mode=-1 or mode=-2 
  !// CMK
 
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
  integer,intent(in)::off_straight   !,off_diag
  real(real8),intent(in)::scep_matrix(:,:)
  real(real8),intent(inout)::straight_vector(:)
!  real(real8),intent(inout)::diag_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  real(real8),parameter::sqrt2=real(sqrt(real(2.0,real8)), real8)
  
  count = off_straight-1
  do i = 1, dimension
      do j = 1,i-1
          count = count + 1
          straight_vector(count) = straight_vector(count) + scep_matrix(j,i)-&
                                                            scep_matrix(i,j)
      enddo
  enddo    
      
end subroutine unscepper_diag_m12

!*****************************************************************

subroutine unscepper_diag(scep_matrix, straight_vector, off_straight, dimension, mode, diag_vector, off_diag)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
  
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
  integer,intent(in)::off_diag,off_straight
  real(real8),intent(in)::scep_matrix(:,:)
  real(real8),intent(inout)::straight_vector(:)
  real(real8),intent(inout)::diag_vector(:)
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  real(real8),parameter::sqrt2=real(sqrt(real(2.0,real8)), real8)
  
  if (mode == 1) then
  
      count = off_straight-1
      do i = 1, dimension
          do j = 1,i-1
                      
              count = count + 1
              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)+&
                                                                scep_matrix(i,j)
          enddo
          diag_vector(i+off_diag-1)  = diag_vector(i+off_diag-1) + scep_matrix(i,i)*sqrt2
      enddo

  elseif (mode == -1 .or. mode == -2) then
      
      count = off_straight-1
      do i = 1, dimension
          do j = 1,i-1
              count = count + 1
              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)-&
                                                                scep_matrix(i,j)
          enddo
      enddo    
      
  endif
  
end subroutine unscepper_diag

!*****************************************************************

subroutine unscepper_diag_sigma_atomic_p1(scep_matrix, straight_vector, off_straight, dimension, off_diag)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
 
  !// Uses diag_vector. Only call if aa_address valid!!!!!!!
  !// mode=+1 (singlet case?)
  !// CMK

  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
  integer,intent(in)::off_diag,off_straight
  real(real8),intent(in)::scep_matrix(:,:)
  type(blockedLockVectorType),intent(inout):: straight_vector
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  integer::startp,endp
  real(real8),parameter::sqrt2 = real(sqrt(real(2.0,real8)), real8)
  real(real8)::tmp
    
  startp = off_straight
  endp = startp+(dimension*dimension-dimension)/2-1
  
#ifdef TIGER_USE_OMP
  call setLocks(straight_vector%l,startp,endp)
#endif

  count = off_straight-1
  do i = 1, dimension
      do j = 1,i-1
                      
          count = count + 1
#ifdef DEBUG_TIGER
          if(count > endp) then
            write(*,*) "scepper(sigma) : ",count,endp
            stop
          endif
#endif
          straight_vector%v(count) = straight_vector%v(count) + scep_matrix(j,i)+&
                                                            scep_matrix(i,j)
      enddo
      tmp = scep_matrix(i,i)*sqrt2
      !$omp atomic
      straight_vector%v(i+off_diag-1)  = straight_vector%v(i+off_diag-1) + tmp
  enddo
#ifdef TIGER_USE_OMP
  call unsetLocks(straight_vector%l,startp,endp)
#endif
  
end subroutine unscepper_diag_sigma_atomic_p1

!*****************************************************************

subroutine unscepper_diag_sigma_atomic_m12(scep_matrix, straight_vector, off_straight, dimension)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP

  !// Does not use diag_vector!
  !// mode=-1 or mode=-2 
  !// CMK
  
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
  integer,intent(in)::off_straight !,off_diag
  real(real8),intent(in)::scep_matrix(:,:)
  type(blockedLockVectorType),intent(inout):: straight_vector
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  integer::startp,endp
    
  startp = off_straight
  endp = startp+(dimension*dimension-dimension)/2-1
  
#ifdef TIGER_USE_OMP
  call setLocks(straight_vector%l,startp,endp)
#endif

  count = off_straight-1
      
  do i = 1, dimension
      do j = 1,i-1
          count = count + 1
#ifdef DEBUG_TIGER
          if(count > endp) then
            write(*,*) "scepper(sigma) : ",count,endp
            stop
          endif
#endif
          straight_vector%v(count) = straight_vector%v(count) + scep_matrix(j,i)-&
                                                            scep_matrix(i,j)
      enddo
  enddo    
#ifdef TIGER_USE_OMP
  call unsetLocks(straight_vector%l,startp,endp)
#endif
  
end subroutine unscepper_diag_sigma_atomic_m12

!*****************************************************************

subroutine unscepper_diag_sigma_atomic(scep_matrix, straight_vector, off_straight, dimension, mode, off_diag)
  
  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
  
  implicit none
  
  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
  integer,intent(in)::off_diag,off_straight
  real(real8),intent(in)::scep_matrix(:,:)
  type(blockedLockVectorType),intent(inout):: straight_vector
  
  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
  integer::startp,endp
  real(real8),parameter::sqrt2 = real(sqrt(real(2.0,real8)), real8)
  real(real8)::tmp
    
  startp = off_straight
  endp = startp+(dimension*dimension-dimension)/2-1
  
#ifdef TIGER_USE_OMP
  call setLocks(straight_vector%l,startp,endp)
#endif

  count = off_straight-1
  if (mode == 1) then
      do i = 1, dimension
          do j = 1,i-1
                      
              count = count + 1
#ifdef DEBUG_TIGER
              if(count > endp) then
                write(*,*) "scepper(sigma) : ",count,endp
                stop
              endif
#endif
              straight_vector%v(count) = straight_vector%v(count) + scep_matrix(j,i)+&
                                                                scep_matrix(i,j)
          enddo
          tmp = scep_matrix(i,i)*sqrt2
          !$omp atomic
          straight_vector%v(i+off_diag-1)  = straight_vector%v(i+off_diag-1) + tmp
      enddo
  elseif (mode == -1 .or. mode == -2) then
      do i = 1, dimension
          do j = 1,i-1
              count = count + 1
#ifdef DEBUG_TIGER
              if(count > endp) then
                write(*,*) "scepper(sigma) : ",count,endp
                stop
              endif
#endif
              straight_vector%v(count) = straight_vector%v(count) + scep_matrix(j,i)-&
                                                                scep_matrix(i,j)
          enddo
      enddo    
  endif
#ifdef TIGER_USE_OMP
  call unsetLocks(straight_vector%l,startp,endp)
#endif
  
end subroutine unscepper_diag_sigma_atomic

!*****************************************************************

!subroutine unscepper_diag_ci_atomic(scep_matrix, straight_vector, off_straight, dimension, mode, diag_vector, off_diag)
!  
!  !// VERY STRAIGHTFORWARD.  LOOPOVER ELEMENTS AND DIVY THEM UP
!  
!  implicit none
!  
!  integer,intent(in)::dimension          !// SIZE OF THE SCEP MATRIX
!  integer,intent(in)::mode               !// FOR SINGLET OR TRIPLET SPIN COUPLING
!  integer,intent(in)::off_diag,off_straight
!  real(real8),intent(in)::scep_matrix(:,:)
!  real(real8),intent(inout)::straight_vector(:)
!  real(real8),intent(inout)::diag_vector(:)
!  
!  integer::i,j,count          !// FOR LOOPING OVER ELEMENTS
!  integer::startp,endp
!  real(real8),parameter::sqrt2 = real(sqrt(real(2.0,real8)), real8)
!  real(real8)::tmp
!  
!  startp = off_straight
!  endp = startp+(dimension*dimension-dimension)/2-1
!  
!#ifdef TIGER_USE_OMP
!  call setCILocks(startp,endp)
!#endif
!
!  count = off_straight-1
!  if (mode == 1) then
!      do i = 1, dimension
!          do j = 1,i-1
!                      
!              count = count + 1
!#ifdef DEBUG_TIGER
!              if(count > endp) then
!                write(*,*) "scepper(ci) : ",count,endp
!                stop
!              endif
!#endif
!              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)+&
!                                                                scep_matrix(i,j)
!          enddo
!          tmp = scep_matrix(i,i)*sqrt2
!          !$omp atomic
!          diag_vector(i+off_diag-1)  = diag_vector(i+off_diag-1) + tmp
!      enddo
!  elseif (mode == -1) then
!      do i = 1, dimension
!          do j = 1,i-1
!              count = count + 1
!#ifdef DEBUG_TIGER
!              if(count > endp) then
!                write(*,*) "scepper(ci) : ",count,endp
!                stop
!              endif
!#endif
!              straight_vector(count) = straight_vector(count) + scep_matrix(j,i)-&
!                                                                scep_matrix(i,j)
!          enddo
!      enddo    
!  endif
!#ifdef TIGER_USE_OMP
!  call unsetCILocks(startp,endp)
!#endif
!  
!end subroutine unscepper_diag_ci_atomic

end module scepmaker_mod
  
  
  
  
  
