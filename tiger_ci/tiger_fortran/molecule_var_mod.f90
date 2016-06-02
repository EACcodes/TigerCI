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
  !//  MOLECULE_VAR_MOD - THIS MODULE IS USED TO PASS AROUND INFORMATION
  !//  REGARDING THE MOLECULE.  IT INCLUDES THINGS LIKE THE SPIN, THE
  !//  NUMBER OF ATOMS, THE NUMBER OF ELECTRONS, ETC. WE ALSO PASS AROUND
  !//  ANY RELEVANT BASIS SET INFORMATION THAT WE MAY NEED.
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  ! WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************
  
module molecule_var_mod
  use global_var_mod, ONLY: real8
  implicit none
  
  integer::num_elec        !// NUMBER OF ELECTRONS
  integer::num_atoms       !// NUMBER OF ATOMS
  integer::spinM           !// SPIN MULTIPLICITY
  integer::num_internal    !// NUMBER OF INTERNAL ORBITALS (OCCUPIED IN A REFERENCE)
  integer::num_external    !// NUMBER OF EXTERNAL ORBITALS (UNOCCUPIED IN REFERENCES)
  integer::num_active      !// NUMBER OF ACTIVE ORBITALS (SINGLY OR VARIABLY OCCUPIED IN REFERENCES)
  integer::num_inactive    !// NUMBER OF INACTIVE ORBITALS (DOUBLE OCCUPIED IN REFERENCES; NUM_ACTIVE + NUM_INACTIVE = 
  			             !// NUM_INTERNAL)
  integer::num_frozen      !// NUMBER OF FROZEN ORBITALS
  integer::nm2csfs         !// TOTAL NUMBER OF N-2 CSFS (DOUBLY EXCITED)			             
  integer::num_extC2       !// NUMBER OF EXTERNALS CHOOSE 2			             
  integer::num_orbitals    !// TOTAL NUMBER OF ORBITALS (NUM_ORBITALS = NUM_ACTIVE + NUM_INACTVE + NUM_EXTERNAL
  			             !//                                        = NUM_INTERNAL + NUM_EXTERNAL
  			             !//                                        = NUMBER_BAS)

  integer::num_orbitalsC2  !// TOTAL NUMBER OF ORBITALS CHOOSE TWO

  integer::number_bas      !// NUMBER OF BASIS FUNCTIONS. SHOULD EQUAL NUM_ORBITALS, OR WE HAVE A PROBLEM.  
  integer::number_basC2    !// NUMBER OF BASIS FUNCTIONS CHOOSE TWO
  integer::num_ref         !// NUMBER OF REFERENCES
  
  integer::numcho          !// MAXIMUM RANK OF CHOLESKY VECTORS
  
  integer, dimension(:,:), allocatable::references  !// THIS ARRAY SHOULD BE DIMENSIONED TO HAVE NUM_REF COLUMNS AND
  						  !// NUM_ORBITALS ROWS.  THE ELEMENTS OF EACH COLUMN WILL THEN
  						  !// CONTAIN THE OCCUPATION PATTERN FOR EACH REFERENCE. 
  Real(kind=real8),dimension(:,:),allocatable::ci_guess     !// USER INPUT WEIGHT AND CI COEFFICIENTS
  
end module molecule_var_mod
