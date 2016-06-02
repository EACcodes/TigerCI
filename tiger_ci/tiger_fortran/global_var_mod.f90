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
  !//  GLOBAL_VAR_MOD - THE GLOBAL VARIABLES MODULE.  STUFF THAT ALMOST ALL
  !//  ROUTINES SHOULD NEED.  HOUSEKEEPING STUFF.
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************
  
module global_var_mod
  use iso_c_binding  
  implicit none
  
  !// FILENAMES, ETC.
  character(:), allocatable::scratch_directory
    
  !// VARIABLE TYPES
  integer, parameter :: real8 = c_double
  integer, parameter :: real4 = c_float
  integer, parameter :: int8  = c_long
  integer, parameter :: int4  = c_int

  !// CONTROL FLAGS
  logical::restart_flag
  integer::valence_ci_flag
  integer::reference_ci_flag
  integer::nat_orb_flag
  integer::nonlocal_flag    
  logical::density_fitting
  logical::CPP_DECOMPOSED_INTS      !// integral decomposed by either DF or CD in C++
  logical::CPP_TRANSFORMED_INTS     !// decomposed integral transformed by either DF or CD in C++
  
  ! Numerical cutoff flags
  real(kind=real8)::integral_threshold
  real(kind=real8)::ao_integral_threshold
  real(kind=real8)::energy_tol         !// CI CONVERGENCE TOLERANCE ON ENERGY
  real(kind=real8)::norm_tol     !// CI CONVERENGCE TOLERANCE ON |HC-EC|

  !// OTHER GLOBAL VARIABLES THAT EVERYONE SHOULD HAVE ACCESS TOO, OR SHOULD
  !// BE PUT IN ONE PLACE SO THAT THEY CAN BE EASILY CHANGED BY THE USER
  integer, parameter::fsn_size=150
  integer, parameter::kill_number = -1000000
  real(real8),parameter::magic_number = real(42,real8)
  
  !> Parameter for different excited state solvers
  integer::num_roots

  !> WP parameters
  real(kind=real8) :: wp_default_radius, wp_multiplier, internal_threshold
  
  !> TOV parameters for occupied orbitals
  real(kind=real8) :: tov_occupied_default_radius, tov_occupied_multiplier, tov_occupied_threshold
  
  !> TOV parameters for virtual orbitals
  real(kind=real8) :: tov_virtual_default_radius, tov_virtual_multiplier, virtual_threshold
  
  !> Cylinder parameters
  real(kind=real8) :: wp_cylinder_radius
  real(kind=real8) :: tov_cylinder_radius
  
  logical :: sphere_based_integral_truncations
  
  !> options for the LOVOs
  logical :: localize_inactive
     
  !// PRINTIING VARIABLES
  integer::diagprint=0        !// PRINT DIAGONAL ELEMENTS
  integer::transprint=0       !// PRINT THE TRANSPOSITIONS
  integer::cycleprint=0       !// PRINT THE CYCLES
  integer::borderprint=0      !// PRINT THE ORBITAL GRAPH BORDER
  integer::drtprint=0         !// PRINT THE DRT
  integer::auxprint=0         !// PRINT THE AUXILIARY VECTOR
  integer::indexprint=0       !// PRINT THE INDEX VECTOR
  integer::sftprint=0         !// PRINT THE SPIN FUNCTION TRANSFORMATION GRAPH
  integer::branchprint=0      !// PRINT THE SPIN FUNCTION BRANCHING DIAGRAM
  integer::graphprint=0       !// PRINT GRAPHICAL REPRESENTATION OF DRT
  integer::pqrr_prqrprint=0   !// PRINT (PQ|RR) AND (PR|QR) INTEGRALS
  integer::ppprint=0          !// PRINT (PP) INTEGRALS
  integer::contractedprint=0  !// PRINT CONTRACTED CYCLES
  logical::sphereprint      !// PRINT INFO ON LOCAL CI ORBITAL SPHERES

  !// ACPF VARIABLE 
  integer::acpf_flag
  logical :: ACPF_ROOT_FOLLOW, ACPF_ROOT_FOLLOW_HELP
  logical:: use_input_ref_energy
  real(kind=real8), dimension(:), allocatable :: user_ref_energy
  real(kind=real8) :: custom_g_val

  !// CHOLESKY DECOMPOSITION VARIABLES
  real(kind=real8)::cd_thresh
  integer::max_mem
  logical :: restart_from_old_CD

  !// INTEGRAL ASSEMBLY MEMORY SIZE
  logical::integralDirect
  logical::fullyIntegralDirect
  logical::directFourInternal
  logical::directLowMem
  logical::directSuperLowMem
  logical::cdVecsInMemory
  integer::max_mem_ints = 500 * 1024 * 1024 / 8
  integer::integral_buffer_size
  
  !// BUFFERED IO VARIABLES
  !// INITIAL POOLIDs ONLY FOR DEBUGGING
  integer::for_buf_blocksizeInts
  integer::for_buf_maxmemInts
  integer::for_buf_maxmemAIJK = 10 * 1024 * 1024
  logical::for_buf_storeIntegrals
  integer::for_buf_int_poolID = -20
  integer::for_buf_integer_int_poolID = -20
  integer::for_buf_aijk_poolID = -55

  integer::for_buf_maxmemCD
  logical::for_buf_storeCDVecs
  integer::for_buf_cd_poolID = -22
  
  integer::for_buf_maxmemSegs
  integer::for_buf_blocksizeSegs
  integer::for_buf_twostwoSeg_poolID = -23
  integer::for_buf_twostwoSegReal_poolID = -24
  integer::for_buf_abcdSeg_poolID = -25
  integer::for_buf_iabcSeg_poolID = -26
  integer::for_buf_twoSeg_poolID = -27
  
  !//NUMBER OF THREADS
#ifdef TIGER_USE_OMP
  integer::numberOfThreads
  integer,parameter::OMP_LOOP_ITER=10
  integer::for_buf_threeiseg_poolID = -56
  integer::for_buf_threeiseg_maxmem = 10*1024*1024 ! keep 10 MB in memory
  integer::for_buf_threeiseg_block = 5*1024 ! really, 5K entries (so 40K blocks) should be fine
#else
  integer,parameter::numberOfThreads=1
#endif
  integer::blockSizeSigma
  integer::blockSizeCI
  
  ! hard drive space for the Davidson subspace vectors (user controlled but defaults to 30)
  ! in GB. Note that this only effects the block Davidson solvers. Doesn't do anything for
  ! the more complicated excited state solvers
  real(kind=real8) :: ASSUMED_DAV_DISK_SPACE = 30 
  
end module global_var_mod
