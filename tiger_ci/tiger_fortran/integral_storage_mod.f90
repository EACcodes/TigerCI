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
  !//  INTEGRAL_STORAGE_MOD - THIS MODULE JUST DEFINES ARRAYS WHICH ARE USED
  !//  TO PASS AROUND BLOCKS OF INTEGRALS.  
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
  ! JMD: HOLY... this is one nasty mod. next target: remove this stuff..
!**************************************************************
  
module integral_storage_mod
  
  use global_var_mod
  
  implicit none
  
  real(real8), dimension(:), allocatable::pp             !// (PP) INTEGRALS
  real(real8), dimension(:), allocatable::pppp           !// (PPPP) INTEGRALS
  real(real8), dimension(:), allocatable::pprr           !// (PP|RR) INTEGRALS
  real(real8), dimension(:), allocatable::prpr           !// (PR|PR) INTEGRALS
  real(real8), dimension(:), allocatable::pqrr           !// (PQ|RR) AND (PQ|QQ) INTEGRALS
  real(real8), dimension(:), allocatable::prqr           !// (PR|QR) INTEGRALS
  real(real8), dimension(:), allocatable::pqrs           !// (PQ|RS) INTEGRALS
  real(real8), dimension(:), allocatable::iaib           !// (IB|AB) INTEGRALS
  real(real8), dimension(:), allocatable::iaja           !// (IA|JA) INTEGRALS
  real(real8), dimension(:), allocatable::ijaa           !// (IJ|AA) INTEGRALS
  real(real8), dimension(:), allocatable::ijij           !// (IJ|IJ) INTEGRALS
  real(real8), dimension(:), allocatable::iaia           !// (IA|IA) INTEGRALS 
 
  real(real8), dimension(:,:), allocatable::integral_buffer     !// (AC|BC),(AB|CD), (IA|BC)
  !real(real4), dimension(:,:), allocatable::integral_buffer_r4
  
  !// FOR USE IN THE VECTORIZED TREATMENT OF THE (AB|CD) INTEGRALS
  !real(real4), dimension(:,:,:), allocatable::abcd_buffer4
  !real(real8), dimension(:,:,:), allocatable::abcd_bufferE
  
  !// THE PSEUDOSPECTRAL QUANTITIES
  real(real8), dimension(:,:),  allocatable::q_mat       !// THE Q MATRIX
  real(real8), dimension(:,:),  allocatable::r_mat       !// THE R MATRIX
  real(real8), dimension(:,:,:),allocatable::Agkl        !// THE Akl(g) INTEGRALS
  real(real8), dimension(:,:),  allocatable::rc          !// R TIMES CI VECTOR
  real(real8), dimension(:,:),  allocatable::arc         !// Agkl TIMES R TIMES CI VECTOR
  real(real8), dimension(:,:,:),  allocatable::rc_gbt    !// R TIMES CI VECTOR FOR EACH T
  real(real8), dimension(:,:,:),  allocatable::arc_tgb   !// Agkl TIMES R TIMES CI VECTOR
  real(real8), dimension(:,:,:),  allocatable::qbrc      !// QB TIMES R TIMES CI VECTOR
  real(real8), dimension(:),  allocatable::Bkq_mat       !// B TIMES Q MATRIX
  
  
                                                         
  !// THIS NEXT VARIABLE WILL NEED TO BE CHANGED 
  
  real(real8), dimension(:,:,:), allocatable::iabc_buffer
 
end module integral_storage_mod
  
  
