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
!> \brief THIS MODULE CONTAINS A POST CI CALCULATION WHICH
!> ATTEMPTS TO EXTRACT THE NATURAL ORBITALS.
!> THE PROCEDURE IS AS FOLLOWS:
!>         1) OBTAIN THE ONE PARTICLE REDUCED DENSITY MATRIX IN
!>            THE MO BASIS (THE BULK OF THE WORK IS HERE SO WE
!>            BREAK THIS INTO MANAGABLE PARTS)
!>               A) COMPUTE DIAGONAL ELEMENTS
!>               B) COMPUTE PURELY EXTERNAL ELEMENTS
!>               C) COMPUTE PURELY INTERNAL ELEMENTS
!>               D) COMPUTE INTERNAL / EXTERNAL ELEMENTS
!>         2) CONVERT THE REDUCED DENSITY MATRIX TO THE AO BASIS
!>         3) DIAGONALIZE THIS MATRIX TO GET THE NATURAL ORBITALS
!>            IN THE AO BASIS
!>         4) USE THESE ORBITALS TO MAP OUT A DENSITY ON A REAL
!>            SPACE GRID
!>
!> THE OUTPUT IS CURENTLY BEING DUMPED TO THE FILE, "N_ORBS.OUT"
!*****************************************************************

! To reduce the size of symbols a few abbreviations have been added
! David Krisiloff 12/22/2010
!
! particle   --> part
! complement --> compl
! external   --> ext

module natural_mod
  
  use global_var_mod
  use molecule_var_mod
  use integral_storage_mod
  use ci_utilities_mod
  use new_tree_search_mod
  use two_seg_mod_3
  use spin_var_mod
  use spin_mod
  use math_utils, ONLY: diagonalize_real_symmetric
  use cholesky_structs
  use locist_mod
  use decider_symbols
  use tree_search_mod
  use sgga_helpers
  use blocked_locks_mod
  
  implicit none
  
  integer, parameter, dimension(5,2,2)::full_loops_natural=&
         reshape((/0,1,&    !// {13} LOOP
                   1,0,&    
                           
                   1,1,&    !// {35} LOOP
                   0,2,&    
                         
                   1,1,&    !// {53} LOOP
                   2,0,&    
                         
                   2,1,&    !// {75} LOOP
                   1,2,&   
                         
                   0,2,&    !// {26} LOOP
                   2,0/),&    
                   shape = (/5,2,2/),&
                   order = (/3,2,1/))
                   
  real(real8), dimension(:,:), allocatable::permut_mat
  real(real8), dimension(:,:), allocatable::one_density_matrix
  real(real8), dimension(:,:), allocatable::one_density_matrix_ao
  real(real8), dimension(:,:), allocatable::coefficient_mat
  
  
  contains
  
  
!*****************************************************************
subroutine natural_orbital_driver(civec)
  
  !// THIS SUBROUTINE IS THE MAIN DRIVER FOR THE GENERATION OF NATURAL
  !// ORBITALS.  THE MAIN CHALLENGE IS TO PRODUCE THE ONE PARTICLE REDUCED 
  !// DENSITY MATRIX IN THE MO BASIS. THE OTHER STEPS (AS YOU WILL SEE 
  !// BELOW) ARE SOMEWHAT TRIVIAL. IN THE END WE WILL OBTAIN THE ONE MATRIX
  !// IN THE AO BASIS AND THE NATURAL ORBITALS IN THE AO BASIS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use locist_var_mod,only:locist_scratch
  use three_four_seg_var_mod
  use two_seg_var_mod
  
  implicit none
  
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec
  
  integer::allocatestatus             !// FOR DYNAMIC MEMORY ALLOCATION
  integer::deallocatestatus
  integer::orbital_a, orbital_b       !// ORBITAL COUNTERS USED FOR PRINTING
  integer::max_spin_dim               !// MAXIMUM NUMBER OF SPIN COUPLINGS
  real(real8)::zero                   !// THE NUMBER ZERO


  !// THESE VARIABLES ARE BEING DEFINED TO GET THE NATURAL ORBITALS.

  real(real8), dimension(:), allocatable::eigen_values
  real(real8), dimension(:,:), allocatable::eigen_vectors
    
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(lambda_path%arc_weights(0:num_orbitals), stat=allocatestatus)
! AUTOGENERATED INITALIZATION
lambda_path%arc_weights = 0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"%2sg_wei")
  
  allocate(lambda_path%occupations(0:num_orbitals), stat=allocatestatus)
! AUTOGENERATED INITALIZATION
lambda_path%occupations = 0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"%lam_occ")
  
  allocate(lambda_path%singles(0:num_orbitals), stat=allocatestatus)
! AUTOGENERATED INITALIZATION
lambda_path%singles = 0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"%lam_sin")
  
  allocate(lambda_path%constraints(3,6), stat=allocatestatus)
! AUTOGENERATED INITALIZATION
lambda_path%constraints = 0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"%lam_con")
  
  allocate(one_density_matrix(num_orbitals,num_orbitals),&
       one_density_matrix_ao(num_orbitals,num_orbitals),&
       coefficient_mat(num_orbitals,num_orbitals),&
       stat = allocatestatus)
! AUTOGENERATED INITALIZATION
one_density_matrix = 0
one_density_matrix_ao = 0
coefficient_mat = 0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"one_density")             

  allocate(eigen_vectors(num_orbitals,num_orbitals),&
       eigen_values(num_orbitals),&
       stat = allocatestatus)
! AUTOGENERATED INITALIZATION
eigen_vectors = 0.0
eigen_values = 0.0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"eigen   ")             

  max_spin_dim = fsn(open_shells)
  allocate(permut_mat(max_spin_dim, max_spin_dim), &
       stat=allocatestatus)
! AUTOGENERATED INITALIZATION
permut_mat = 0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"permut     ")         


  !// INITIALIZE THE REDUCED DENSITY MATRIX
  zero = real(0.0, real8)
  one_density_matrix = zero

  !// FORMATS FOR OUTPUTTING OUT RESULTS IN AN ORDERLY MANNER 
!1650 format(4f19.15)
1200 format(i5) 
1211 format(3f13.9) 
!1222 format(13i3) 
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !// WE START BY BUILDING THE REDUCED DENSITY MATRIX PIECE BY PIECE. HOPEFULLY
  !// THE LABELS ARE SELF EXPLANATORY.


  !!!!!!!!!!!!!!!!!
  !//
  !// DIAGONAL ELEMENTS
  !//
  !!!!!!!!!!!!!!!!!

  call one_particle_diagonal(civec, lambda_path)


  !!!!!!!!!!!!!!!!!
  !// 
  !// PURELY EXTERNAL PART
  !//
  !!!!!!!!!!!!!!!!!

  call one_particle_purely_external(civec, lambda_path)                        


  !!!!!!!!!!!!!!!!!
  !// 
  !// PURELY INTERNAL PART
  !//
  !!!!!!!!!!!!!!!!!

  call one_particle_purely_internal(civec, lambda_path)                   
  

  !!!!!!!!!!!!!!!!!
  !// 
  !// INTERNAL / EXTERNAL PART
  !//
  !!!!!!!!!!!!!!!!!  

  call one_particle_internal_external(civec,lambda_path)


  !// WE NOW HAVE THE ENTIRE ONE PARTICLE REDUCED DENSITY MATRIX 
  !// IN THE MOLECULAR ORBITAL BASIS. OUR NEXT ORDER OF BUSINESS IS
  !// TO CONVERT THIS MATRIX TO THE AO BASIS.  START BY READING IN 
  !// THE MO COEFFICIENTS.  THEN PERFORM THE TRANSFORMATION.


  open(unit=15,file=scratch_directory // "coefficients", access="sequential" )
  read(15,*)

  do orbital_a=1, num_orbitals
      read(15,*)
      read(15,*) (coefficient_mat(orbital_a,orbital_b), orbital_b=1, num_orbitals)
  enddo



  !// HERE'S THE TRANSFORM (YEP, IT IS THAT SIMPLE)
  
  one_density_matrix_ao = matmul(one_density_matrix, coefficient_mat)
  coefficient_mat = transpose(coefficient_mat)
  one_density_matrix_ao = matmul(coefficient_mat, one_density_matrix_ao)



  !// DIAGONALIZE THE MOLECULAR ORBITAL REDUCED DENSITY MATRIX TO GET
  !// THE NATURAL ORBITALS IN THE MO BASIS.  
  !// HERE I SIMPLY CALL A SUBROUTINE THAT SHOULD DO THIS FOR ME.  I 
  !// NEED TO TEST THIS GARBAE SO LETS SEE HOW THINGS GO AFTERWARDS.

  call diagonalize_real_symmetric(one_density_matrix, eigen_values, &
       eigen_vectors, num_orbitals)


  !// NOW WE CAN MAP OUR DENSITY ONTO A REAL SPACE GRID.  HOWEVER, 
  !// THIS PROGRAM IS NOT YET CAPABLE OF DOING THIS.  PLEASE STAY
  !// ON THE LINE AND WE WILL HELP YOU AS SOON AS POSSIBLE...
  !// YOUR BUSINESS IS IMPORTANT TO US... PLEASE STOP BY OFTEN TO 
  !// VIEW POSSIBLE UPDATES. WE APOLOGIZE FOR ANY INCONVIENCE...

  print*, "natural orbitals complete"

  !//  PRINT THE RESULTS TO THE OUTPUT FILE


  open(unit=14,file=scratch_directory // "n_orbs.out", access="sequential" )

  write(14,*) "MO one electron reduced density matrix"
  write(14,*) 
  do orbital_a=1, num_orbitals
      write(14,1200) orbital_a 
      write(14,1211) (one_density_matrix(orbital_a,orbital_b),orbital_b=1,num_orbitals) 
      write(14,*)
  enddo

  write(14,*) "AO one electron reduced density matrix"
  write(14,*) 
  do orbital_a=1, num_orbitals
      write(14,1200) orbital_a 
      write(14,1211) (one_density_matrix_ao(orbital_a,orbital_b),orbital_b=1,num_orbitals) 
      write(14,*) 
  enddo


!  write(14,*) "MO natural orbitals"
!  write(14,*) 
!  do orbital_a=1, num_orbitals
!      write(14,1200) orbital_a 
!      write(14,1211) ((eigen_vectors(orbital_a,orbital_b)),orbital_b=1,num_orbitals) 
!      write(14,*) 
!  enddo


  !// CONVERT TO AO BASIS

  
!  coefficient_mat = transpose(coefficient_mat)
!  eigen_vectors = matmul(coefficient_mat,eigen_vectors)
!
!  write(14,*) "AO natural orbitals"
!  write(14,*) 
!  do orbital_a=1, num_orbitals
!      write(14,1200) orbital_a 
!      write(14,1211) ((eigen_vectors(orbital_a,orbital_b)),orbital_b=1,num_orbitals) 
!      write(14,*) 
!  enddo


  !// RE-ORDERED PRINT OF NATURAL ORBITALS IN AO BASIS


  eigen_vectors = matmul(coefficient_mat,eigen_vectors)

  write(14,*) "AO natural orbitals after re-ordering"
  write(14,*) 
  do orbital_a=1, num_orbitals
      write(14,1200) orbital_a 
      write(14,1211) (eigen_vectors(orbital_b, num_orbitals - orbital_a + 1), &
           orbital_b=1,num_orbitals) 
      write(14,*) 
  enddo



  !// CLEAN UP BEFORE WE EXIT THE ROUTINE 
  close(14) 

  !deallocate(civec, stat = deallocatestatus)
  !call deallocatecheck(deallocatestatus,"civec   ")

  deallocate(one_density_matrix, one_density_matrix_ao, coefficient_mat,&
           stat = allocatestatus)
  call deallocatecheck(allocatestatus,"one_density")             

  deallocate(eigen_vectors,eigen_values, &
       stat = deallocatestatus)
  call deallocatecheck(allocatestatus,"eigen   ")             

  deallocate(lambda_path%arc_weights, stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"%2sg_wei")
  
  deallocate(lambda_path%occupations, stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"%lam_occ")
  
  deallocate(lambda_path%singles, stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"%lam_sin")
  
  deallocate(lambda_path%constraints, stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"%lam_con")
  
  deallocate(permut_mat,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"permut     ")               


end subroutine natural_orbital_driver
  
  
!*****************************************************************
subroutine one_particle_diagonal(civec, lambda_path)
  
  !// HERE WE OBTAIN THE DIAGONAL ELEMENTS OF THE REDUCED DENSITY MATRIX.    
  !// WE KNOW THE ELEMENTS OF THE REDUCED DENSITY MATRIX CAN BE RE-WRITTEN
  !// IN TERMS OF CREATION AND ANNIHILATION OPERATORS.  FOR THE DIAGONAL 
  !// ELEMENTS WE ANNIHILATE AN ELECTRON IN ORBITAL "I" AND RE-CREATE IT IN 
  !// ORBITAL "I".  THEREFORE, OUR LINE UP PERMUTATIONS (FROM THE SPIN PARTS
  !// OF OUR WAVEFUNCTION) ARE EQUAL TO THE IDENTITY THEREBY SIMPLIFYING OUR
  !// CALCULATION (BY MAKING THE DIAGONAL ELEMENTS DEPENDENT UPON THE CI 
  !// COEFFICIENTS ONLY).
  
  use locist_var_mod,only:locist_scratch
  use three_four_seg_var_mod

  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec
  type(graph_search_state) :: graph
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// INITIALIZE DERIVED TYPE ARRAYS
  lambda_path%occupations = 0
  lambda_path%singles = 0
  lambda_path%arc_weights = 0
  lambda_path%constraints = -1
  lambda_path%num_singles = 0

 
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !// AGAIN WE DO NOT HAVE TO WORRY ABOUT SPIN SO WE CAN LOOP THROUGH ALL POSSIBLE 
  !// CSFS TO GET THE RESPECTIVE CI COEFFICIENTS.  TO FOLLOW ALONG IN A LOGICAL 
  !// MANNER WE DIVIDE UP THE CSFS INTO THREE GROUPS: NO EXCITATIONS, ONE EXCITATION, 
  !// AND TWO EXCITATIONS. TREESEARCH WILL FIND THE INTERNAL PART OF OUR CONFIGURATION 
  !// PATH BEFORE CALLING THE NEXT TASK NAMELY: "one_part_diagonal_compl".

  !// NO EXCITATIONS: SEARCH THE INTERNAL SPACE
  call init_tree_search(graph, 0, 0, num_internal, num_elec)
  do while ( get_internal_path(lambda_path, graph))
      call one_part_diagonal_compl(civec, lambda_path)
  end do


  !// ONE EXCITATION: SEARCH THE INTERNAL SPACE
  
  !// REINITIALIZE
  lambda_path%occupations = 0
  lambda_path%singles = 0
  lambda_path%arc_weights = 0
  lambda_path%constraints = -1
  lambda_path%num_singles = 0
  
  call init_tree_search(graph, 0, 0, num_internal, num_elec-1)
  do while ( get_internal_path(lambda_path, graph))
      call one_part_diagonal_compl(civec, lambda_path)
  end do
  
  
  !// TWO EXCITATIONS: SEARCH THE INTERNAL SPACE
  
  !// REINITIALIZE
  lambda_path%occupations = 0
  lambda_path%singles = 0
  lambda_path%arc_weights = 0
  lambda_path%constraints = -1
  lambda_path%num_singles = 0
  
  call init_tree_search(graph, 0, 0, num_internal, num_elec-2)
  do while ( get_internal_path(lambda_path, graph))
      call one_part_diagonal_compl(civec, lambda_path)
  end do

  
end subroutine one_particle_diagonal
  
!*****************************************************************
subroutine one_part_diagonal_compl(civec, lambda_path)

  !// IN THIS SUBROUTINE WE ALREADY HAVE THE INTERNAL CONFIGURATION
  !// PATH.  WE FORM THE EXTERNAL PARTS AND POINT TO THE CORRESPONDING
  !// CI COEFFICIENTS.  WITH THE CI COEFFICIENTS IN HAND WE COMPUTE THE 
  !// CONFIGURATION'S CONTRIBUTION TO THE DIAGONAL ELEMENTS OF THE ONE-
  !// PARTICLE REDUCED DENSITY MATRIX.
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec
  
  
  integer::start_elec             !// NUMBER OF ELECTRONS IN THE INTERNAL SPACE
  integer::b,c                    !// LABELS THE LOOP LEVELS
  integer::n                      !// LOOP COUNTER
  integer::internal_weight        !// WEIGHT OF INTERNAL PART OF PATH
  integer::internal_singles       !// NUMBER OF SINGLES IN INTERNAL PART OF PATH
  integer::internal_address       !// ADDRESS OF INTERNAL SPACE CSF IN INTERNAL INDEX
  integer::lambda_singles         !// NUMBER OF SINGLES IN PATHS
  integer::lambda_dim             !// DIMENSIONS OF MU PATH
  integer::lambda_address         !// PATH ADDRESSES
  integer::num_virtual            !// NUMBER OF VIRTUAL ORBITALS
  integer::ab_count               !// FOR LOOPING OVER A,B
  integer::spin_lambda            !// FOR LOOPING OVER SPIN FUNCTIONS
  real(real8)::two                !// THE NUMBER TWO
! START AUTOGENERATED INITIALIZATION 
ab_count = 0
two = 0.0
internal_weight = 0
lambda_dim = 0
spin_lambda = 0
lambda_singles = 0
internal_singles = 0
start_elec = 0
num_virtual = 0
lambda_address = 0
b = 0
c = 0
n = 0
internal_address = 0
! END AUTOGENERATED INITIALIZATION 
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  start_elec = sum(lambda_path%occupations(0:num_internal))
  internal_weight = sum(lambda_path%arc_weights(0:num_internal))
  internal_singles = lambda_path%num_singles
  num_virtual = num_orbitals-num_internal
  two = real(2.0, real8)
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!
  !//  
  !// FIRST TREAT THE NUM_ELEC-0 CASE
  !//
  !!!!!!!!!!!!!!!!!!!!!!!!


  !// SINCE THESE ARE REFERENCE STATES, THERE IS NO NEED TO LOOP OVER THE EXTERNAL
  !// ELEMENTS. THE CALCULATION OF THE MATRIX ELEMENTS SHOULD BE STRAIGHT FORWARD
  !// FROM THE FOLLOWING EQUATION.

  startelec: if (start_elec == num_elec - 0) then

      !// SET UP INTERNAL ADDRESS FOR OUR PATH
      internal_address = internal_index_vector0(internal_weight+1)

      if (internal_address > 0) then

          !// SET UP SINGLES AND DIMENSIONS
          lambda_singles     = internal_singles 
          lambda_dim         = fsn(lambda_singles)

          if (lambda_dim > 0) then
              do spin_lambda = 1, lambda_dim
                  lambda_address = internal_address + (spin_lambda-1) 
                
                  do n=1, num_internal
                      one_density_matrix(n,n) = one_density_matrix(n,n) + &
                           civec%v(lambda_address) * &
                           civec%v(lambda_address) * lambda_path%occupations(n)
                  enddo
              enddo
          endif
      endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!
  !//  
  !// NEXT TREAT THE NUM_ELEC-1 CASE
  !//
  !!!!!!!!!!!!!!!!!!!!!!!!
 
  elseif (start_elec == num_elec - 1) then

      
      !// SET UP INTERNAL ADDRESS FOR OUR PATH

      internal_address = internal_index_vector1(internal_weight+1)
      if (internal_address > 0) then
      
          !// SET UP SINGLES AND DIMENSIONS
          lambda_singles     = internal_singles + 1
          lambda_dim         = fsn(lambda_singles)
      
          !// COMPUTE MATRIX ELEMENTS AND DO MULTIPLICATIONS
          if (lambda_dim > 0) then
              do spin_lambda = 1, lambda_dim

                  lambda_address = internal_address + (spin_lambda-1)*num_virtual
                  do c = 1, num_external
                      
                      !//  COMPUTE DIAGONAL ELEMENTS OF THE INTERNAL SPACE
                      do n=1, num_internal
                          one_density_matrix(n,n) = one_density_matrix(n,n) + &
                               civec%v(lambda_address - 1 + c) * &
                               civec%v(lambda_address - 1 + c) * & 
                               lambda_path%occupations(n)
                      enddo

                      !//  COMPUTE DIAGONAL ELEMENTS OF VIRTUAL SPACE
                      n = num_internal + c
                      
                      one_density_matrix(n,n) = one_density_matrix(n,n) +  &
                           civec%v(lambda_address - 1 + c) * &
                           civec%v(lambda_address - 1 + c) 


                  enddo
              enddo
          endif
      endif


  !!!!!!!!!!!!!!!!!!!!!!!!
  !//  
  !// NOW THE NUM_ELEC-2 CASE
  !//
  !!!!!!!!!!!!!!!!!!!!!!!!    

  elseif (start_elec == num_elec - 2) then startelec
 
      internal_address = internal_index_vector2(internal_weight+1)
      if (internal_address > 0) then

          !!!!!!!!!!!!!!!!!!!
          !//
          !// TWO VIRTUALS SINGLY OCCUPIED
          !//
          !!!!!!!!!!!!!!!!!!!

     
          !// SET UP SINGLES AND DIMENSIONS
          lambda_singles     = internal_singles + 2
          lambda_dim         = fsn(lambda_singles)
          
          
          !// COMPUTE MATRIX ELEMENTS AND DO MULTIPLICATIONS
          if (lambda_dim > 0) then

              
              !// COMPUTE MATRIX ELEMENTS AND DO MULTIPLICATIONS
              do spin_lambda = 1, lambda_dim
                  lambda_address = internal_address + (spin_lambda-1)*num_extC2
                  
                  ab_count = 0
                  do c = 1, num_external
                      do b = 1, c-1
                          ab_count = ab_count+1
                          
                          !//  COMPUTE DIAGONAL ELEMENTS OF THE INTERNAL SPACE
                          do n=1, num_internal
                              one_density_matrix(n,n) = &
                                   one_density_matrix(n,n) + &
                                   civec%v(lambda_address - 1 + ab_count) * &
                                   civec%v(lambda_address - 1 + ab_count) * &
                                   lambda_path%occupations(n)                          
                          enddo
                          
                          !//  COMPUTE DIAGONAL ELEMENTS OF VIRTUAL SPACE
                          n = num_internal + b
                          one_density_matrix(n,n) = &
                               one_density_matrix(n,n) + &
                               civec%v(lambda_address - 1 + ab_count) * &
                               civec%v(lambda_address - 1 + ab_count)                           

                          n = num_internal + c
                          one_density_matrix(n,n) = &
                               one_density_matrix(n,n) + &
                               civec%v(lambda_address - 1 + ab_count) * &
                               civec%v(lambda_address - 1 + ab_count)                           

                      enddo
                  enddo
              enddo
          endif

              
          !!!!!!!!!!!!!!!!!!!
          !//
          !// ONE VIRTUAL DOUBLY OCCUPIED
          !//
          !!!!!!!!!!!!!!!!!!!
  

          !// SET UP SINGLES AND DIMENSIONS
          lambda_singles     = internal_singles 
          lambda_dim         = fsn(lambda_singles)

                  
          !// COMPUTE MATRIX ELEMENTS AND DO MULTIPLICATIONS
          if (lambda_dim > 0) then
                  
              do spin_lambda = 1, lambda_dim
                  
                  lambda_address = internal_index_vector3(internal_weight+1) + &
                       (spin_lambda - 1) * num_external
                  
                  do c = 1, num_external
                          
                      !//  COMPUTE DIAGONAL ELEMENTS OF THE INTERNAL SPACE
                      do n=1, num_internal
                          one_density_matrix(n,n) = &
                               one_density_matrix(n,n) + &
                               civec%v(lambda_address + c - 1) * &
                               civec%v(lambda_address - 1 + c) * &
                               lambda_path%occupations(n)  
                      enddo

                      
                      !// COMPUTE EXTERNAL DIAGONAL ELEMENTS.  HERE I MULTIPLY BY 
                      !// 2.00 BECAUSE I KNOW THESE ARE DOUBLY OCCUPIED 
                      n = num_internal + c

                      one_density_matrix(n,n) = &
                           one_density_matrix(n,n) + &
                           civec%v(lambda_address + c - 1) * &
                           civec%v(lambda_address - 1 + c) * two
                      
                  enddo
              enddo
          endif
          
      endif
      
  endif startelec       

end subroutine one_part_diagonal_compl

!*****************************************************************
subroutine one_particle_purely_external(civec, lambda_path)
  
  !// THE CONTRIBUTIONS TO THE EXTERNAL PART OF THE 1-PARTICLE
  !// REDUCED DENSITY MATRIX COME FROM THE TWO SEGMENT LOOPS LYING
  !// WITHIN THE EXTERNAL SPACE.  AFTER FINDING  A GIVEN INTERNAL 
  !// PATH WE FIND ALL SUCH EXTERNAL LOOPS, MATCH UP THE CI 
  !// COEFFICIENTS, AND FIND THE CORRESPONDING LINE-UP PERMUTATION 
  !// RESULTING FROM THE SPIN.
  
  use locist_var_mod,only:locist_scratch
  use three_four_seg_var_mod
  use two_seg_var_mod
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::lambda_path
  type(graph_search_state) :: graph
  type(blockedLockVectorType) :: civec

  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// INITIALIZE DERIVED TYPE ARRAYS
  lambda_path%occupations = 0
  lambda_path%singles = 0
  lambda_path%arc_weights = 0
  lambda_path%constraints = -1
  lambda_path%num_singles = 0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !// FIRST WE BUILD THE INTERNAL PATHES.  THEN, WE DIRECT THE ROUTINE TO 
  !// THE NEXT TASK, "one_part_purely_ext_compl".
  
  
  !// FIRST SEARCH TO NUM_ELEC-2
  call init_tree_search(graph, 0, 0, num_internal, num_elec-2)
  do while ( get_internal_path(lambda_path, graph))
      call one_part_purely_ext_compl(civec, lambda_path)
  end do
  
  !// REINITIALIZE DERIVED TYPE
  lambda_path%occupations = 0
  lambda_path%singles = 0
  lambda_path%arc_weights = 0
  lambda_path%constraints = -1
  lambda_path%num_singles = 0

  
  !// NOW SEARCH TO NUM_ELEC-1
  call init_tree_search(graph, 0, 0, num_internal, num_elec-1)
  do while ( get_internal_path(lambda_path, graph))
      call one_part_purely_ext_compl(civec, lambda_path)
  end do
  

end subroutine one_particle_purely_external
    
!*****************************************************************
subroutine one_part_purely_ext_compl(civec, lambda_path)

  !// IN THIS SUBROUTINE WE EXTRACT THE PURELY EXTERNAL PART OF THE 
  !// 1-PARTICLE REDUCED DENSITY MATRIX BY BUILDING THE TWO SEGMENT
  !// LOOPS WHICH RESIDE ENTIRELY IN THE EXTERNAL SPACE AND JOINING
  !// THEM TO THE (PASSED IN) INTERNAL PATH  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec
  
  
  integer::start_elec             !// NUMBER OF ELECTRONS IN THE INTERNAL SPACE
  integer::a,b,c                  !// LABELS THE LOOP LEVELS
  integer::internal_weight        !// WEIGHT OF INTERNAL PART OF PATH
  integer::internal_singles       !// NUMBER OF SINGLES IN INTERNAL PART OF PATH
  integer::internal_address       !// ADDRESS OF INTERNAL SPACE CSF IN INTERNAL INDEX
  integer::lambda_singles         !// NUMBER OF SINGLES IN PATHS
  integer::mu_singles
  integer::lambda_dim             !// DIMENSIONS OF MU PATH
  integer::mu_dim
  integer::lambda_address         !// PATH ADDRESSES
  integer::mu_address
  integer::allocatestatus         !// FOR DYNAMIC MEMORY
  integer::deallocatestatus       
  integer::ibar, jbar             !// THE POSITIONS OF THE UNPAIRED SINGLES
  integer::num_virtual            !// NUMBER OF VIRTUAL ORBITALS
  integer::ab_count               !// FOR LOOPING OVER A,B
  integer::spin_lambda,spin_mu    !// FOR LOOPING OVER SPIN FUNCTIONS
  integer::spin_max_dim           !// MAXIMUM NUMBER OF SPIN COUPLINGS

  real(real8)::zero               !// THE NUMBER ZERO

  real(real8), dimension(:,:), allocatable::identity
! START AUTOGENERATED INITIALIZATION 
mu_address = 0
mu_dim = 0
mu_singles = 0
spin_lambda = 0
a = 0
b = 0
c = 0
spin_mu = 0
lambda_dim = 0
internal_address = 0
lambda_address = 0
jbar = 0
lambda_singles = 0
ab_count = 0
! END AUTOGENERATED INITIALIZATION 


    
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  start_elec = sum(lambda_path%occupations(0:num_internal))
  internal_weight = sum(lambda_path%arc_weights(0:num_internal))
  internal_singles = lambda_path%num_singles
  num_virtual = num_orbitals-num_internal
  spin_max_dim = fsn(open_shells)
  zero = real(0.0, real8)
  
  if (skip_this_internal(internal_weight, start_elec)) return
  
  !// GET THE IDENTITY MATRIX IN THE REQUIRED DIMENSIONS
  allocate(identity(spin_max_dim, spin_max_dim),&
       stat = allocatestatus)
  call allocatecheck(allocatestatus,"identity")             


  identity = zero
  identity = identity_matrix(spin_max_dim)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !// THE INTERNAL PATH CAN END IN TWO LOCATIONS NAMELY: NUM_ELEC-1 AND NUM_ELEC-2.
  !// FOR THE NUM_ELEC-1 CASE ONLY THE {13} LOOP IS POSSIBLE.  FOR THE THE OTHER CASE 
  !// (INTERNAL PATH ENDS AT NUM_ELEC-2) WE CAN BUILD THE {13}, {35}, AND {53} 
  !// LOOPS.  WE SEPARATE THE TWO CASES FOR THE MAXIMUM CLARITY.


 
  !!!!!!!!!!!!!!!!!!!!!!!!
  !//  
  !// FIRST TREAT THE NUM_ELEC-1 CASE
  !//
  !!!!!!!!!!!!!!!!!!!!!!!!
  startelec: if (start_elec == num_elec - 1) then

      
      !!!!!!!!!!!!!!!!!!!
      !//
      !// JUST THE {13} LOOP IS POSSIBLE
      !//
      !!!!!!!!!!!!!!!!!!!

      !// DEFINE THE INTERNAL ADDRESS FOR THE PATH BEING USED.
      internal_address = internal_index_vector1(internal_weight+1)

      if (internal_address > 0) then
      
          !// SET UP SINGLES AND DIMENSIONS
          lambda_singles     = internal_singles + 1
          mu_singles         = lambda_singles
          lambda_dim         = fsn(lambda_singles)
          mu_dim             = lambda_dim
      
         
          if (lambda_dim > 0) then


              !// LOOP OVER ALL POSSIBLE COMBINATIONS OF SPIN.  THEN LOOP OVER THE POSSIBLE 
              !// EXTERNAL CONFIGURATIONS THAT MAKE 13 LOOPS.
              do spin_lambda = 1, lambda_dim

                  lambda_address = internal_address + (spin_lambda-1)*num_virtual
                  do spin_mu = 1, mu_dim

                      mu_address = internal_address + (spin_mu-1)*num_virtual

                      do a = 1, num_external
                          do b = 1, a-1
                              
                              !// ADD IN THE CONTRIBUTION TO THE REDUCED DENSITY MATRIX 
                              !// RESULTING FROM THE RESPECTIVE CSFs.  FOR THESE LOOPS IN
                              !// PARTICULAR, THE LINE-UP PERMUTATIONS ARE THE IDENTITY.  
                              !// CONSEQUENTLY, WE JUST MULTIPLY THROUGH BY THE IDENTITY 
                              !// RATHER THAN COMPUTING THE PERMUTATION MATRIX (WHICH WE
                              !// KNOW WOULD BE THE IDENTITY ANYWAYS).

                              one_density_matrix(a+num_internal, b+num_internal) = &
                                   one_density_matrix(a+num_internal, b+num_internal) + &
                                   civec%v(lambda_address - 1 + a) * &
                                   civec%v(mu_address - 1 + b) * &
                                   identity(spin_lambda,spin_mu)

                              one_density_matrix(b+num_internal, a+num_internal) = &
                                   one_density_matrix(b+num_internal, a+num_internal)  + &
                                   civec%v(lambda_address - 1 + b) * &
                                   civec%v(mu_address - 1 + a) * &
                                   identity(spin_lambda,spin_mu)


                          enddo
                      enddo
                  enddo
              enddo            
          endif
     
      endif


  !!!!!!!!!!!!!!!!!!!!!!!!
  !//  
  !// NOW THE NUM_ELEC-2 CASE
  !//
  !!!!!!!!!!!!!!!!!!!!!!!!    
  elseif (start_elec == num_elec - 2) then startelec
  
      
      !// ONCE AGAIN DEFINE THE INTERNAL ADDRESS FOR THE PATH BEING USED.
      internal_address = internal_index_vector2(internal_weight+1)
 
     if (internal_address > 0) then
          
          !// SET UP SINGLES AND DIMENSIONS
          lambda_singles     = internal_singles + 2
          mu_singles         = lambda_singles
          lambda_dim         = fsn(lambda_singles)
          mu_dim             = lambda_dim


          if (lambda_dim > 0) then


              !// GET THE LINE UP PERMUTATIONS IN A MATRIX REPRESENTATION.
              !// THIS MATRIX IS COMPUTED IN CONTRACTED_CYCLE BUT WE ARE
              !// REQUIRED TO COMPUTE AND PASS IN IBAR AND JBAR.

              jbar = internal_singles + 1
              ibar = internal_singles + 2

              !permut_mat = zero
              call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
              !permut_mat = contracted


              !!!!!!!!!!!!!!!!!!!
              !//
              !// {13} LOOPS
              !//
              !!!!!!!!!!!!!!!!!!!


              !// LOOP OVER THE POSSIBLE SPIN AND SPATIAL PARTS TO GET A NEW SET OF CSFs.
              do spin_lambda = 1, lambda_dim
                  lambda_address = internal_address + (spin_lambda-1)*num_extC2
                  
                  do spin_mu = 1, mu_dim
                      mu_address = internal_address + (spin_mu-1)*num_extC2

                      do a = 1, num_external

                          do b = a+1, num_external

                              do c = b+1, num_external

                                  one_density_matrix(a+num_internal, c+num_internal) = &
                                       one_density_matrix(num_internal+a, num_internal+c) + &
                                       civec%v(lambda_address - 1 + index2m1(a,b)) * &
                                       civec%v(mu_address - 1 + index2m1(b,c)) * &
                                       permut_mat(spin_lambda,spin_mu)

                                  one_density_matrix(c+num_internal, a+num_internal) = &
                                       one_density_matrix(num_internal+c, num_internal+a) + &
                                       civec%v(lambda_address - 1 + index2m1(b,c)) * &
                                       civec%v(mu_address - 1 + index2m1(a,b)) * &
                                       permut_mat(spin_lambda,spin_mu)

                                  

                                  !// AGAIN THE LINE UP PERMUTATIONS FOR THE FOLLOWING TERMS 
                                  !// HAPPEN TO JUST BE THE IDENTITY.  THEREFORE, WE PERFORM
                                  !// OUR MULTIPLICATIONS APPROPRIATELY.


                                  one_density_matrix(a+num_internal, b+num_internal) = &
                                       one_density_matrix(num_internal+a, num_internal+b) + &
                                       civec%v(lambda_address - 1 + index2m1(a,c)) * &
                                       civec%v(mu_address - 1 + index2m1(b,c)) * &
                                       identity(spin_lambda, spin_mu)
                                      
                                  one_density_matrix(b+num_internal, c+num_internal) = &
                                       one_density_matrix(num_internal+b, num_internal+c) + &
                                       civec%v(lambda_address - 1 + index2m1(a,b)) * &
                                       civec%v(mu_address - 1 + index2m1(a,c)) * &
                                       identity(spin_lambda, spin_mu)
                                  
                                  one_density_matrix(b+num_internal, a+num_internal) = &
                                       one_density_matrix(num_internal+b, num_internal+a) + &
                                       civec%v(lambda_address - 1 + index2m1(b,c)) * &
                                       civec%v(mu_address - 1 + index2m1(a,c)) * &
                                       identity(spin_lambda, spin_mu)
                                      
                                  one_density_matrix(c+num_internal, b+num_internal) = &
                                       one_density_matrix(num_internal+c, num_internal+b) + &
                                       civec%v(lambda_address - 1 + index2m1(a,c)) * &
                                       civec%v(mu_address - 1 + index2m1(a,b)) * &
                                       identity(spin_lambda, spin_mu)


                              enddo
                          enddo
                              
                      enddo
                  enddo
              enddo
          
              
              !!!!!!!!!!!!!!!!!!!
              !//
              !// {35}, {53} LOOPS
              !//
              !!!!!!!!!!!!!!!!!!!  

              !// SET UP SINGLES AND DIMENSIONS.  WE DEFINE LAMBDA TO BE THE 
              !// PATH WITH TWO SINGLES IN THE EXTERNAL SPACE.
              lambda_singles     = internal_singles + 2
              mu_singles         = internal_singles
              lambda_dim         = fsn(lambda_singles)
              mu_dim             = fsn(mu_singles)


              !// COMPUTE THE PERMUTATION MATRIX FOR THESE LOOPS.  THIS REQUIRES
              !// US TO RE-DEFINE THE SIZE OF OUR PERMUTATION MATRIX.
              jbar = internal_singles + 1
              ibar = internal_singles + 2

              !permut_mat = zero
              call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
              !permut_mat = contracted


              !// THE MU PATH HAS ONE EXTERNAL DOUBLY OCCUPIED.  THEREFORE WE HAVE TO DEFINE
              !// THE PATH ADDRESS ACCORDINGLY.
              mu_address = internal_index_vector3(internal_weight+1)

              if (mu_address > 0) then
                  
                  if (lambda_dim > 0 .and. mu_dim > 0) then
                      

                      !// NOW FIND ALL THE CSFS THAT FORM 35/53 LOOPS BY LOOPING OVER THE
                      !// SPIN PARTS AND THE EXTERNAL ORBITALS.
                      do spin_lambda = 1, lambda_dim
                          lambda_address = internal_address + (spin_lambda-1)*num_extc2
                          
                          do spin_mu = 1, mu_dim
                              mu_address = internal_index_vector3(internal_weight+1) + &
                                   (spin_mu - 1) * num_external

                              ab_count = 0
                              do a = 1, num_external
                                  
                                  do b=1, a-1
                                      ab_count = ab_count + 1


                                      !// ADD IN THE CONTRIBUTION OF THE GIVEN 35/53 LOOPS TO 
                                      !// THEIR RESPECTIVE MATRIX ELEMENTS.  A FACTOR OF SQUARE
                                      !// ROOT OF TWO IS ADDED IN FOR NORMALIZATION.
                                      
                                      one_density_matrix(a+num_internal, b+num_internal) = &
                                           one_density_matrix(a+num_internal, b+num_internal) + &
                                           civec%v(lambda_address + ab_count - 1) * &
                                           civec%v(mu_address - 1 + b) * &
                                           permut_mat(spin_lambda,spin_mu)
                                      
                                      one_density_matrix(b+num_internal, a+num_internal) = &
                                           one_density_matrix(b+num_internal, a+num_internal) + &
                                           civec%v(lambda_address - 1 + ab_count) * &
                                           civec%v(mu_address - 1 + b) * &
                                           permut_mat(spin_lambda,spin_mu)
                                      
                                      one_density_matrix(a+num_internal, b+num_internal) = &
                                           one_density_matrix(a+num_internal, b+num_internal) + &
                                           civec%v(lambda_address + ab_count - 1) * &
                                           civec%v(mu_address - 1 + a) * &
                                           permut_mat(spin_lambda,spin_mu)
                                  
                                      one_density_matrix(b+num_internal, a+num_internal) = &
                                           one_density_matrix(b+num_internal, a+num_internal) + &
                                           civec%v(lambda_address - 1 + ab_count) * &
                                           civec%v(mu_address - 1 + a) * &
                                           permut_mat(spin_lambda,spin_mu)
                                  enddo
                              enddo
                          enddo
                      enddo
                  endif
              endif
          endif
      endif

  endif startelec
      

  !// CLEAN UP
  deallocate(identity, stat = deallocatestatus)
  call deallocatecheck(allocatestatus,"identity")             

  
end subroutine one_part_purely_ext_compl
  
  
!********************************************************************************
subroutine one_particle_internal_external(civec, lambda_path)
  
  !// MAIN DRIVER ROUTINE FOR TWO SEGMENT LOOPS WITH ONE SEGMENT IN THE 
  !// INTERNAL SPACE AND THE OTHER IN THE EXTERNAL SPACE.
  
  use locist_var_mod,only:locist_scratch
  use three_four_seg_var_mod
  use two_seg_var_mod
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec

  
  integer::i                        !// LABELS LOOP LEVELS 
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// INITIALIZE DERIVED TYPE ARRAYS
  lambda_path%occupations = 0
  lambda_path%singles = 0
  lambda_path%arc_weights = 0
  lambda_path%constraints = -1
  lambda_path%num_singles = 0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// LOOP OVER THE INTERNAL WHERE THE LOOP IS GOING TO BE FORMED
  do i = 1, num_internal
  
  
      !!!!!!!!!!!!!!!!!
      !//
      !// FIRST THE {1} SEGMENT
      !//
      !!!!!!!!!!!!!!!!!
      
      !// INITIALIZE DERIVED TYPE ARRAYS
      lambda_path%occupations = 0
      lambda_path%singles = 0
      lambda_path%arc_weights = 0
      lambda_path%constraints = -1
      lambda_path%num_singles = 0
      
      !// STORE THE LOOP TYPE
      lambda_path%loop_type = 1
                  
      !// NOW LOAD CONSTRAINTS UP AND TRY TO BUILD THE LOOPS
      lambda_path%constraints(2,6) = 1  !// CONSTRAINT COUNT
      lambda_path%constraints(1,6) = INT_EXT_COMP !// NEXT TASK
      lambda_path%constraints(2,2) = num_internal+1
                          
      !// FOR BUILDING THE MU PATH
      lambda_path%level1 = i
  
      !// LOOP LEVELS
      lambda_path%constraints(2,1) = i
                          
      !// LAMBDA, MU PATHS
      lambda_path%constraints(1,1) = 0
      lambda_path%constraints(3,1) = 1
                                  
      !// NOW DO THE SEARCH
      lambda_path%num_singles = 0
      call broken_constrained_search(lambda_path,0,0,a=civec)      
      
  
      !!!!!!!!!!!!!!!!!
      !//
      !// NOW THE {2} SEGMENT
      !//
      !!!!!!!!!!!!!!!!!
      
      !// INITIALIZE DERIVED TYPE ARRAYS
      lambda_path%occupations = 0
      lambda_path%singles = 0
      lambda_path%arc_weights = 0
      lambda_path%constraints = -1
      lambda_path%num_singles = 0
      
      !// STORE THE LOOP TYPE
      lambda_path%loop_type = 2
                  
      !// NOW LOAD CONSTRAINTS UP AND TRY TO BUILD THE LOOPS
      lambda_path%constraints(2,6) = 1  !// CONSTRAINT COUNT
      lambda_path%constraints(1,6) = INT_EXT_COMP !// NEXT TASK
      lambda_path%constraints(2,2) = num_internal+1
                          
      !// FOR BUILDING THE MU PATH
      lambda_path%level1 = i
  
      !// LOOP LEVELS
      lambda_path%constraints(2,1) = i
                          
      !// LAMBDA, MU PATHS
      lambda_path%constraints(1,1) = 1
      lambda_path%constraints(3,1) = 2
                                  
      !// NOW DO THE SEARCH
      lambda_path%num_singles = 0
      call broken_constrained_search(lambda_path,0,0,a=civec)
      
  enddo    
  
end subroutine one_particle_internal_external
  
      
!*****************************************************************
subroutine one_part_internal_ext_compl(civec,lambda_path)
  
  !// IN THIS ROUTINE WE BUILD THE COMPLEMENT TO THE TWO SEGMENT LOOPS
  !// WHERE ONE SEGMENT RESIDES IN THE INTERNAL SPACE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec
  
  integer::ibar,jbar              !// ORBITAL POSITION INDICES
  integer::internal_singles       !// NUMBER OF SINGLES IN INTERNAL PATH OF PATH
  integer::lambda_singles         !// SINGLES IN LAMBDA AND MU PATH
  integer::mu_singles
  integer::lambda_dim,mu_dim      !// DIMENSIONS FOR LAMBDA AND MU CONFIGURATIONS
  integer::lambda_address         !// ADDRESS OF MU AND LAMBDA PATHS
  integer::mu_address
  integer::current_vertex         !// USED IN BUILDING MU PATH
  integer::step_type              !// KEEPS TRACK OF STEP IN MU PATH CONSTRUCTION
  integer::path_elecs             !// CUMULATIVE OCCUPATION IN MU PATH
  integer::top_level              !// LEVEL OF FIRST LOOP SEGMENT
  integer::mu_levels              !// LABELS LEVELS IN MU PATH
  integer::constraint_count       !// KEEP TRACK OF CONSTRAINTS IN BUILDING MU PATH
  integer::start_elec             !// OCCUPATION WHERE INTERNAL CSF MEETS UP WITH EXTERNAL SPACE
  integer::a,b                    !// USED TO INDEX EXTERNAL ORBITALS
  integer::ab_count
  integer::s1,s2,s3,s4            !// FOR THE CALL TO GET NUMBERS OF SINGLES
  integer::loop_type              !// LABELS THE LOOP
  integer::i,j                    !// LABELS THE INTERNAL LOOP LEVELS
  integer::k                      !// LOOP INDEX
  
  !// IN THE FOLLOWING v => VALENCE, ab => TWO VIRTS SINGLY OCCPIED,
  !// aa => ONE VIRT DOUBLY OCCUPIED, a => ONE VIRT SINGLY OCUPIED
  integer::v_address,a_address,ab_address,aa_address
  integer::ab_spin 
  integer::a_spin
  integer::v_spin
  integer::ab_start,ab_end 
  integer::aa_start,aa_end 
  integer::a_start,a_end 
  integer::v_start,v_end
  integer::lambda_spin
  integer::mu_spin
  integer::lambda_start,lambda_end 
  integer::mu_start,mu_end 
  
  real(real8)::zero                       !// THE NUMBER ZERO 
! START AUTOGENERATED INITIALIZATION 
ab_spin = 0
i = 0
ab_address = 0
ab_count = 0
a_address = 0
zero = 0.0
ab_end = 0
mu_address = 0
a_start = 0
v_start = 0
a_spin = 0
internal_singles = 0
mu_levels = 0
top_level = 0
aa_address = 0
lambda_end = 0
aa_end = 0
a = 0
j = 0
k = 0
lambda_spin = 0
v_spin = 0
loop_type = 0
jbar = 0
mu_spin = 0
ibar = 0
mu_end = 0
lambda_singles = 0
v_end = 0
aa_start = 0
lambda_address = 0
lambda_dim = 0
current_vertex = 0
v_address = 0
mu_singles = 0
mu_dim = 0
mu_start = 0
ab_start = 0
start_elec = 0
constraint_count = 0
step_type = 0
lambda_start = 0
s4 = 0
s2 = 0
s3 = 0
s1 = 0
b = 0
a_end = 0
path_elecs = 0
! END AUTOGENERATED INITIALIZATION 

  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  zero = real(0.0, real8)

  loop_type = lambda_path%loop_type
  internal_singles = lambda_path%num_singles
  i = lambda_path%constraints(2,1)
  j = lambda_path%level1
  k = lambda_path%constraints(1,1)
  start_elec = sum(lambda_path%occupations(0:num_internal))
  
  !// INITIALIZE VARIABLES FOR BUILDING MU PATH
  top_level = lambda_path%level1
  constraint_count = 1
  lambda_path%rt_loop_weight = 0
  path_elecs = sum(lambda_path%occupations(0:top_level-1))
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// 
  !//  BUILD THE MU PATH
  !// 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do mu_levels = top_level, num_internal
  
      current_vertex = vertex(mu_levels-1,path_elecs)
  
      !// CHECK TO SEE IF THIS IS A CONSTRAINED LEVEL AND GET THE STEP_TYPE
      if ((mu_levels)==lambda_path%constraints(2,constraint_count)) then
          step_type = lambda_path%constraints(3,constraint_count)
          constraint_count = constraint_count + 1
      else
          step_type = lambda_path%occupations(mu_levels)
      endif
      
      !// TRY TO ADD THE STEP
      if (step_type == 2) then
          
          if (add2(mu_levels-1,path_elecs)) then
              path_elecs = path_elecs + 2
              lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                                           abs(y2(current_vertex))

          else
              return
          endif
          
      elseif(step_type == 1) then
  
          if (add1(mu_levels-1,path_elecs)) then
              path_elecs = path_elecs + 1
              lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                                           abs(y1(current_vertex))

          else
              return
          endif
          
      elseif(step_type == 0) then
  

          if (.not.add0(mu_levels-1,path_elecs)) return
  
      endif
      
  enddo  
  
  !// STORE THE WEIGHT OF THE LAMBDA_PATH
  lambda_path%weight = sum(lambda_path%arc_weights(0:num_internal))
  
  !// FOR THE RT_LOOP_WEIGHT ADD IN THE WEIGHTS OF THE RELEVANT PARTS
  !// OF THE HEAD AND TAIL PATHS.  THIS WILL END UP BEING THE MU PATH WEIGHT
  lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                               sum(lambda_path%arc_weights(0:lambda_path%level1-1))
                               
  if (skip_these_internals(lambda_path%weight, lambda_path%rt_loop_weight,&
                           start_elec, path_elecs)) return
                           
  !!!!!!!!!!!!!!!
  !//
  !// GET JBAR... THIS WILL BE THE CORRECT JBAR REGARDLESS OF THE LOOP TYPE,
  !//             NUMBER OF EXCITATIONS, ETC.
  !//
  !!!!!!!!!!!!!!!                             

  call get_numbers_of_singles(s1,s2,s3,s4,lambda_path)
  jbar = s1+1
  if(loop_type == 2) jbar = s1


  !////////////////////////////////////////////////////////////////////////////////
  !//
  !//  FINISHED COMPUTING THE INTERNAL PATHES... BEGIN WITH NUM_ELEC-1 CASE
  !//
  !////////////////////////////////////////////////////////////////////////////////

  if (start_elec == num_elec-1) then    

  

      !// SETUP THE ADDRESSES FOR THE MU AND LAMBDA PATHES SO THAT WE CAN POINT
      !// TO THE CI COEFFICIENTS THAT CORRESPOND TO OUR FORMED CSFS

      call get_address(start_elec,v_address,a_address,ab_address,aa_address, &
                               lambda_path%weight+1, lambda_path%rt_loop_weight+1)


      !//////////////////
      !//               
      !// BUILD {13} and {53} LOOPS 
      !//               
      !//////////////////


      !// GET THE PERMUTATION MATRIX
      
      if (loop_type == 1) then
          ibar = internal_singles + 1
          lambda_singles = internal_singles + 1
          mu_singles = internal_singles           
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)

          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted

      elseif (loop_type == 2) then
          ibar = internal_singles + 1
          lambda_singles = internal_singles + 1
          mu_singles = internal_singles - 1
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)

          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      endif

      
      if (lambda_dim /=0 .and. mu_dim /=0 ) then 
          
          if (a_address > 0 .and. v_address > 0) then

              !// FORM CSFS AND COMPUTE THE LOOP'S CONTRIBUTION
              !// TO THE REDUCED DENSITY MATRIX.
              
              do a_spin =1,lambda_dim
                  a_start = a_address + num_external*(a_spin - 1) 
                  a_end   = a_start   + num_external - 1 
                  do v_spin =1,mu_dim 
                      v_start = v_address + v_spin -1
                      v_end   = v_start 

                      do a=1, num_external
                          
                          one_density_matrix(i,a+num_internal) = &
                               one_density_matrix(i,a+num_internal) + &
                               civec%v(a_start + a -1) * civec%v(v_start) * &
                               permut_mat(a_spin, v_spin) 

                          one_density_matrix(a+num_internal,i) = &
                               one_density_matrix(a+num_internal,i) + &
                               civec%v(a_start + a -1) * civec%v(v_start) * &
                               permut_mat(a_spin, v_spin) 
                        
                      enddo


                  enddo
              enddo
             
          endif
      
      endif
  

      !///////////////////////////////////////////////////////////////////////
      !//
      !//      NOW WE FIND THE CONTRIBUTION FROM THE VIRTUALS
      !//              THAT CONTAIN TWO ELECTRONS!
      !//
      !///////////////////////////////////////////////////////////////////////


  elseif (start_elec == num_elec-2) then    


      call get_address(start_elec,v_address,a_address,ab_address,aa_address, &
                               lambda_path%weight+1, lambda_path%rt_loop_weight+1)


      !////////////////////////////////////////////////////////////////////
      !//
      !// FIRST LET'S TREAT THE LOOPS WHERE ONE EXTERNAL IS DOUBLY OCCUPIED
      !//
      !////////////////////////////////////////////////////////////////////

      !// GET THE PERMUTATION MATRIX FOR THE CASE WHERE ONE VIRTUAL IS DOUBLY OCCUPIED
      
      if (loop_type == 1) then
          !// WE HAVE A {35} LOOP, BUT WE MUST SWAP OUR MU AND LAMBDA PATHES
          ibar = internal_singles + 2
          lambda_singles = internal_singles + 2  
          mu_singles = internal_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)

          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted

      elseif (loop_type == 2) then
          !// WE HAVE A {75} LOOP
          ibar = internal_singles 
          lambda_singles = internal_singles
          mu_singles = internal_singles
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix)
          permut_mat = -permut_mat
          !permut_mat = -contracted

      endif


      if (lambda_dim /= 0 .and. mu_dim /= 0) then 
         
          if ((ab_address > 0 .or. aa_address > 0) .and. a_address > 0) then
                  

              do lambda_spin = 1,lambda_dim
                   
                  !// SETUP THE LAMBDA AND MU ADDRESSES TO REFLECT THE LOOP TYPE
                  !// WE HAVE TO DO THIS HERE BECAUSE WE SWAP THE LAMBDA AND MU 
                  !// PATHES FOR THE {35} LOOPS. 
                  
                  if (loop_type == 1) then
                      lambda_address = a_address
                      mu_address = aa_address
                  else
                      lambda_address = aa_address
                      mu_address = a_address
                  endif

                  lambda_start = lambda_address + num_external*(lambda_spin - 1)
                  lambda_end   = lambda_start   + num_external - 1 


                  do mu_spin = 1,mu_dim 
                      
                      mu_start = mu_address + num_external*(mu_spin - 1)
                      mu_end   = mu_start   + num_external - 1 


                      do a=1, num_external


                          one_density_matrix(i,a+num_internal) = &
                               one_density_matrix(i,a+num_internal) + &
                               civec%v(lambda_start + a -1) * civec%v(mu_start + a -1) * &
                               permut_mat(lambda_spin, mu_spin)


                          one_density_matrix(a+num_internal,i) = &
                               one_density_matrix(a+num_internal,i) + &
                               civec%v(lambda_start + a -1) * civec%v(mu_start + a -1) * &
                               permut_mat(lambda_spin, mu_spin) 

                      enddo
                  enddo
              enddo

          endif
      endif


      !////////////////////////////////////////////////////
      !//
      !// NOW TREAT THE TWO VIRTUALS SINGLY OCCUPIED CASE 
      !//
      !////////////////////////////////////////////////////
      

      
      if ((ab_address > 0 .or. aa_address > 0) .and. a_address > 0) then
          
          
          !// GET PERMUTATION MATRIX (FOR THIS CASE IBAR EQUALS INTERNAL SINGLES + 1.
          !// LATER WE TREAT THE CASE WHERE IBAR EQUALS INTERNAL SINGLES + 2 BUT WE
          !// SEPARATE THE TWO CASES TO EASILY OBTAIN THE CORRECT PERMUTATION MATRIX.)

          if (loop_type == 1) then
              !// THIS IS A {13} LOOP
              ibar = internal_singles + 1
              lambda_singles = internal_singles + 2
              mu_singles = lambda_singles
              lambda_dim = fsn(lambda_singles)
              mu_dim = fsn(mu_singles)
              
              !permut_mat = zero
              call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
              !permut_mat = contracted
              
          elseif (loop_type == 2) then
              !// THIS IS A {53} LOOP              
              ibar = internal_singles + 1
              lambda_singles = internal_singles + 2
              mu_singles = internal_singles 
              lambda_dim = fsn(lambda_singles)
              mu_dim = fsn(mu_singles)
              
              !permut_mat = zero          
              call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
              !permut_mat = contracted
              
          endif
          
          if (lambda_dim /= 0 .and. mu_dim /= 0) then 

              !// LOOP OVER THE SPIN AND EXTERNAL ORBITALS TO FORM THE CSFs, AND 
              !// THEN CALCULATE THEIR CONTRIBUTION TO THE ONE DENSITY MATRIX.
              do ab_spin = 1,lambda_dim
                   
                  ab_start = ab_address + num_extC2*(ab_spin-1) 
                  ab_end   = ab_start   + num_extC2 - 1 
                  
                  aa_start = aa_address + num_external*(ab_spin-1) 
                  aa_end   = aa_start   + num_external - 1 
  
                  do a_spin = 1,mu_dim 
                                  
                      a_start = a_address + num_external*(a_spin - 1)
                      a_end   = a_start   + num_external - 1 


                      ab_count = 0
                      do a=1,num_external
                          do b=1,a-1 
                              ab_count = ab_count+1


                              one_density_matrix(i,b+num_internal) = & 
                                   one_density_matrix(i,b + num_internal) + &
                                   civec%v(ab_start + ab_count - 1) * civec%v(a_start + a -1) * &
                                   permut_mat(ab_spin, a_spin)
                              
                              one_density_matrix(b+num_internal,i) = &
                                   one_density_matrix(b + num_internal,i) + &
                                   civec%v(ab_start + ab_count - 1) * civec%v(a_start + a -1) * &
                                   permut_mat(ab_spin, a_spin)
                                                            
                          enddo
                      enddo
                  enddo
              enddo
          endif

          !// GET PERMUTATION MATRIX (IBAR = INTERNAL SINGLES +2)
          if (loop_type == 1) then
              !// THIS IS A {13} LOOP
              ibar = internal_singles + 2
              lambda_singles = internal_singles + 2
              mu_singles = lambda_singles
              lambda_dim = fsn(lambda_singles)
              mu_dim = fsn(mu_singles)
              
              !permut_mat = zero
              call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
              !permut_mat = contracted
              
          elseif (loop_type == 2) then
              !// THIS IS A {53} LOOP
              ibar = internal_singles + 2
              lambda_singles = internal_singles + 2
              mu_singles = internal_singles 
              lambda_dim = fsn(lambda_singles)
              mu_dim = fsn(mu_singles)
              
              !permut_mat = zero
              call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
              !permut_mat = contracted
              
          endif

          if (lambda_dim /= 0 .and. mu_dim /= 0) then 

              !// LOOP OVER THE SPIN AND EXTERNAL ORBITALS TO FORM THE CSFs, AND 
              !// THEN CALCULATE THEIR CONTRIBUTION TO THE ONE DENSITY MATRIX.
              do ab_spin = 1,lambda_dim
                   
                  ab_start = ab_address + num_extC2*(ab_spin-1) 
                  ab_end   = ab_start   + num_extC2 - 1 
                  
                  aa_start = aa_address + num_external*(ab_spin-1) 
                  aa_end   = aa_start   + num_external - 1 
  
                  do a_spin = 1,mu_dim 
                                  
                      a_start = a_address + num_external*(a_spin - 1)
                      a_end   = a_start   + num_external - 1 


                      ab_count = 0
                      do a=1,num_external
                          do b=1,a-1 
                              ab_count = ab_count+1


                              one_density_matrix(i,a+num_internal) = & 
                                   one_density_matrix(i,a + num_internal) + &
                                   civec%v(ab_start + ab_count - 1) * civec%v(a_start + b -1) * &
                                   permut_mat(ab_spin, a_spin)

                              one_density_matrix(a+num_internal,i) = &
                                   one_density_matrix(a + num_internal,i) + &
                                   civec%v(ab_start + ab_count - 1) * civec%v(a_start + b -1) * &
                                   permut_mat(ab_spin, a_spin)
                              
                          enddo
                      enddo
                  enddo
              enddo
          endif
      endif
  endif    
  
end subroutine one_part_internal_ext_compl

!*****************************************************************
subroutine get_address(start_elec,v_address,a_address,ab_address,aa_address, &
                       lambda_weight, mu_weight)
  
  implicit none
  
  !// THIS ROUTINE COMPUTES ADDRESSES FOR VALENCE, N-1 and N-2 STATES.
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::start_elec 
  integer::mu_weight, lambda_weight
  integer::v_address,a_address,ab_address,aa_address
! START AUTOGENERATED INITIALIZATION 
! END AUTOGENERATED INITIALIZATION 
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (start_elec == num_elec - 1) then 
  
      v_address      = internal_index_vector0(mu_weight)
      a_address      = internal_index_vector1(lambda_weight)
  
  elseif (start_elec == num_elec - 2) then 
  
      a_address      = internal_index_vector1(mu_weight)
      ab_address     = internal_index_vector2(lambda_weight)
      aa_address     = internal_index_vector3(lambda_weight)
  
  endif 
  
end subroutine get_address
  
  
!*****************************************************************
subroutine one_particle_purely_internal(civec, lambda_path)
  
  !// IN THIS ROUTINE WE TREAT THE PURELY INTERNAL TWO SEGMENT
  !// LOOPS.  THIS WILL NOT INCLUDE THE SPECIAL {26} LOOPS BECAUSE 
  !// THEY DIFFER BY MORE THAN ONE ORBITAL AND WILL GO TO ZERO IN
  !// THE REDUCED DENSITY MATRIX DUE TO ORTHOGONALITY
  
  use locist_var_mod,only:locist_scratch
  use three_four_seg_var_mod
  use two_seg_var_mod
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec

  
  integer::i,j                        !// LABELS LOOP LEVELS
  integer::loop                       !// LABELS THE LOOP
  
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// INITIALIZE DERIVED TYPE ARRAYS
  lambda_path%occupations = 0
  lambda_path%singles = 0
  lambda_path%arc_weights = 0
  lambda_path%constraints = -1
  lambda_path%num_singles = 0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// LOOP OVER PAIRS OF INTERNALS I > J

  do i = 1, num_internal
      
      do j = 1, i-1
      
          do loop = 1,4
              
              !// REINITIALIZE
              lambda_path%occupations = 0
              lambda_path%singles = 0
              lambda_path%arc_weights = 0
              lambda_path%constraints = -1
              lambda_path%num_singles = 0
                      
              !// STORE THE LOOP TYPE
              lambda_path%loop_type = loop
                  
              !// NOW LOAD CONSTRAINTS UP AND TRY TO BUILD THE LOOPS
              lambda_path%constraints(2,6) = 1  !// CONSTRAINT COUNT
              lambda_path%constraints(1,6) = PURELY_INT_COMP !// NEXT TASK
              lambda_path%constraints(2,3) = num_internal+1
                          
              !// FOR BUILDING THE MU PATH
              lambda_path%level1 = j
  
              !// LOOP LEVELS
              lambda_path%constraints(2,1:2) = (/j,i/)
                          
              !// LAMBDA, MU PATHS
              lambda_path%constraints(1,1:2) = &
                  full_loops_natural(loop,1,1:2)
                  
              lambda_path%constraints(3,1:2) = &
                  full_loops_natural(loop,2,1:2)
                                  
              !// NOW DO THE SEARCH
              lambda_path%num_singles = 0
              call broken_constrained_search(lambda_path,0,0,a=civec)
                          
          enddo
          
      enddo
  enddo        
    
end subroutine one_particle_purely_internal
  
  
!*****************************************************************
subroutine one_part_purely_internal_compl(civec, lambda_path)
  
  !// IN THIS ROUTINE WE BUILD THE COMPLEMENT TO THE TWO
  !// SEGMENT LOOPS WHICH RESIDE ENTIRELY IN THE INTERNAL SPACE.
  !// 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use spin_mod 
  
  implicit none
  
  type(orbital_path)::lambda_path
  type(blockedLockVectorType) :: civec
  
  integer::ibar,jbar               !// ORBITAL POSITION INDICES
  integer::internal_singles        !// NUMBER OF SINGLES IN INTERNAL PATH OF PATH
  integer::lambda_singles          !// SINGLES IN LAMBDA AND MU PATH
  integer::mu_singles
  integer::lambda_dim,mu_dim       !// DIMENSIONS FOR LAMBDA AND MU CONFIGURATIONS
  integer::current_vertex          !// USED IN BUILDING MU PATH
  integer::step_type               !// KEEPS TRACK OF STEP IN MU PATH CONSTRUCTION
  integer::path_elecs              !// CUMULATIVE OCCUPATION IN MU PATH
  integer::top_level               !// LEVEL OF FIRST LOOP SEGMENT
  integer::mu_levels               !// LABELS LEVELS IN MU PATH
  integer::constraint_count        !// KEEP TRACK OF CONSTRAINTS IN BUILDING MU PATH
  integer::start_elec              !// OCCUPATION WHERE INTERNAL CSF MEETS UP WITH EXTERNAL SPACE
  integer::a,b                     !// USED TO INDEX EXTERNAL ORBITALS
  integer::s1,s2,s3,s4             !// FOR THE CALL TO GET NUMBERS OF SINGLES
  integer::loop_type               !// LABELS THE LOOP
  integer::i,j                     !// LABELS THE INTERNAL LOOP LEVELS
  
  integer::v1_address,v2_address   !// ADDRESSES FOR VALENCE-VALENCE INTERACTIONS
  integer::a1_address,a2_address   !// ADDRESSES FOR N-1/N-1 INTERACTIONS
  integer::v1_spin,v2_spin         !// INDICES FOR VALENCE STATE SPINS
  integer::a1_spin,a2_spin         !// INDICES FOR N-1 STATE SPINS
  integer::v1_start,v2_start       !// STARTING INDEX FOR CI AND SIGMA VECTOR FOR VALENCE STATES
  integer::a1_start,a2_start       !// STARTING INDEX FOR CI AND SIGMA VECTOR FOR N-1  STATES 
  integer::v1_end,v2_end           !// ENDING INDEX FOR CI AND SIGMA VECTOR FOR VALENCE STATES 
  integer::a1_end,a2_end           !// ENDING INDEX FOR CI AND SIGMA VECTOR FOR N-1 STATES 
  
  integer::aa_lam_add,aa_mu_add    !// ADDRESS FOR LAMBDA AND MU PATHS FOR ONE DOUBLY OCC. VIRTUAL
  integer::ab_lam_add,ab_mu_add    !// ADDRESS FOR LAMBDA AND MU PATHS FOR TWO SINGLY OCC. VIRTUALS
  integer::aa1_spin,aa2_spin       !// SPIN INDEX FOR ONE DOUBLY OCC. VIRTUAL 
  integer::aa_lam_start,aa_lam_end !// STARTING AND ENDING PORTIONS OF CI AND SIGMA VECTOR FOR LAMBDA PATH 
                                   !// FOR ONE DOUBLY OCC. VIRTUAL
  integer::ab_lam_start,ab_lam_end !// STARTING AND ENDING PORTIONS OF CI AND SIGMA VECTOR FOR LAMBDA PATH 
                                   !// FOR TWO SINGLY OCC. VIRTUALS
  integer::aa_mu_start,aa_mu_end   !// STARTING AND ENDING PORTIONS OF CI AND SIGMA VECTOR FOR MU PATH 
                                   !// FOR ONE DOUBLY OCC. VIRTUAL
  integer::ab_mu_start,ab_mu_end   !// STARTING AND ENDING PORTIONS OF CI AND SIGMA VECTOR FOR MU PATH 
                                   !// FOR TWO SINGLY OCC. VIRTUALS
  integer::ab_count                !// COUNT FOR A,B 
  integer::index_ab                !// Nv CHOSE 2
  
  real(real8)::zero                !// THE NUMBER ZERO
! START AUTOGENERATED INITIALIZATION 
ab_mu_start = 0
internal_singles = 0
a2_end = 0
v1_start = 0
lambda_singles = 0
top_level = 0
v1_end = 0
a2_address = 0
i = 0
b = 0
a = 0
aa_mu_end = 0
constraint_count = 0
a1_spin = 0
path_elecs = 0
a1_address = 0
loop_type = 0
aa_mu_start = 0
ab_lam_end = 0
aa_lam_start = 0
aa_lam_end = 0
aa_lam_add = 0
step_type = 0
aa1_spin = 0
v2_spin = 0
ab_lam_start = 0
zero = 0.0
ab_count = 0
a2_spin = 0
current_vertex = 0
v2_end = 0
mu_singles = 0
ab_mu_end = 0
v1_spin = 0
a1_start = 0
s1 = 0
s3 = 0
s2 = 0
s4 = 0
mu_dim = 0
start_elec = 0
ibar = 0
v2_address = 0
a2_start = 0
jbar = 0
ab_lam_add = 0
v2_start = 0
a1_end = 0
lambda_dim = 0
mu_levels = 0
ab_mu_add = 0
index_ab = 0
aa_mu_add = 0
v1_address = 0
aa2_spin = 0
j = 0
! END AUTOGENERATED INITIALIZATION 
      
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  index_ab = num_external*(num_external-1)/2
  
  zero = real(0.0, real8)
  
  loop_type = lambda_path%loop_type
  internal_singles = lambda_path%num_singles
  j = lambda_path%constraints(2,1)
  i = lambda_path%constraints(2,2)

  start_elec = sum(lambda_path%occupations(0:num_internal))
  
  !// INITIALIZE VARIABLES FOR BUILDING MU PATH
  top_level = lambda_path%level1
  constraint_count = 1
  lambda_path%rt_loop_weight = 0
  path_elecs = sum(lambda_path%occupations(0:top_level-1))
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// 
  !//  BUILD THE MU PATH
  !// 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do mu_levels = top_level, num_internal
  
      current_vertex = vertex(mu_levels-1,path_elecs)
  
      !// CHECK TO SEE IF THIS IS A CONSTRAINED LEVEL AND GET THE STEP_TYPE
      if ((mu_levels)==lambda_path%constraints(2,constraint_count)) then
          step_type = lambda_path%constraints(3,constraint_count)
          constraint_count = constraint_count + 1
      else
          step_type = lambda_path%occupations(mu_levels)
      endif
      
      !// TRY TO ADD THE STEP
      if (step_type == 2) then
          
          if (add2(mu_levels-1,path_elecs)) then
              path_elecs = path_elecs + 2
              lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                   abs(y2(current_vertex))

          else
              return
          endif
          
      elseif(step_type == 1) then
  
          if (add1(mu_levels-1,path_elecs)) then
              path_elecs = path_elecs + 1
              lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                   abs(y1(current_vertex))
          else
              return
          endif
          
      elseif(step_type == 0) then

          if (.not.add0(mu_levels-1,path_elecs)) return
          
      endif
      
  enddo
  
  !// STORE THE WEIGHT OF THE LAMBDA_PATH
  lambda_path%weight = sum(lambda_path%arc_weights(0:num_internal))
  
  !// FOR THE RT_LOOP_WEIGHT ADD IN THE WEIGHTS OF THE RELEVANT PARTS
  !// OF THE HEAD AND TAIL PATHS.  THIS WILL END UP BEING THE MU PATH WEIGHT
  lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
       sum(lambda_path%arc_weights(0:lambda_path%level1-1))
                               
  if (skip_these_internals(lambda_path%weight, lambda_path%rt_loop_weight, &
       start_elec,path_elecs)) return
                               
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !//
  !// GET SINGLES AND CONVERT TO IBAR, JBAR
  !//
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call get_numbers_of_singles(s1,s2,s3,s4,lambda_path)
  
  if (loop_type == 1.or.loop_type == 4) then
      ibar = s1 + s2
      jbar = s1 + 1
  else
      ibar = s1 + s2
      jbar = s1 
  endif


  !// NOW THAT THE MU PATH IS BUILT (AND WE HAVE IBAR AND JBAR) LETS CALCULATE
  !// THE PURELY INTERNAL PART OF THE REDUCED DENSITY MATRIX RESULTING FROM OUR
  !// LAMBDA PATH.  START WITH THE DOUBLE EXCITATION CASES.

  if (start_elec == num_elec-2) then

  
      !!!!!!!!!!!!!!!!!
      !// 
      !// TWO VIRTUALS SINGLY OCCUPIED             
      !// 
      !!!!!!!!!!!!!!!!!
      
      !// GET DIMENSIONS AND PERMUTATIONS

      if (loop_type == 1) then
          !// {13} LOOP
          lambda_singles = internal_singles + 2
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)  

          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 2 .or. loop_type == 3) then
          !// {35} or {53} LOOP
          lambda_singles = internal_singles + 2
          mu_singles = internal_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 4) then
          !// {75} LOOP
          lambda_singles = internal_singles + 2
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix)
          permut_mat = -permut_mat
          !permut_mat = -contracted
      endif



      !// GET THE ADDRESSES THAT POINT TO THE APROPRIATE SECTION OF THE CI VECTOR
      ab_lam_add = internal_index_vector2(lambda_path%weight +1)
      ab_mu_add  = internal_index_vector2(lambda_path%rt_loop_weight +1)
      

      !// HERE ARE SOME SAFETY CHECKS TO MAKE SURE WE SHOULD CONTINUE ONWARD
      if (lambda_dim /= 0.and.mu_dim /= 0) then
      
          if (ab_lam_add > 0 .and. ab_mu_add > 0) then 
          
              
              do aa1_spin  = 1,lambda_dim 
                  
                  ab_lam_start = ab_lam_add   + num_extC2*(aa1_spin - 1)
                  ab_lam_end   = ab_lam_start + num_extC2 - 1 
                  do aa2_spin = 1, mu_dim
                      
                      ab_mu_start = ab_mu_add   + num_extC2*(aa2_spin - 1)
                      ab_mu_end   = ab_mu_start + num_extC2 - 1 
                      
                      !// COMPUTE THE REDUCED DENSITY MATRIX.  NOW I AM USING
                      !// THE PERMUTATION MATRIX THAT I COMPUTED EARLIER. 
                      
                      ab_count = 0
                                            
                      do a=1, num_external
                          do b=1, a-1
                              ab_count = ab_count + 1
                              
                              
                              !// HERE'S THE ACTUAL CALCULATION OF THIS INTERNAL'S CONTRIBUTION
                              !// TO THE ONE PARTICLE REDUCED DENSITY MATRIX.
                              one_density_matrix(i,j) = &
                                   one_density_matrix(i,j) + &
                                   permut_mat(aa1_spin,aa2_spin) * &
                                   civec%v(ab_lam_start + ab_count -1) * &
                                   civec%v(ab_mu_start + ab_count -1) 
                              
                              one_density_matrix(j,i) = &
                                   one_density_matrix(j,i) + &
                                   permut_mat(aa1_spin,aa2_spin) * &
                                   civec%v(ab_lam_start + ab_count -1) * &
                                   civec%v(ab_mu_start + ab_count -1) 
                          enddo
                          
                      enddo
                  enddo
              enddo
              
          endif
          
      endif
      
      

      !!!!!!!!!!!!!!!!!!!
      !//                
      !// ONE VIRTUAL DOUBLY OCCUPIED              
      !//                                          
      !!!!!!!!!!!!!!!!!!!
      
      !// GET DIMENSIONS AND PERMUTATIONS

      if (loop_type == 1) then
          !// {13} LOOP
          lambda_singles = internal_singles 
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)

          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 2 .or. loop_type == 3) then
          !// {35} or {53} LOOP
          lambda_singles = internal_singles 
          mu_singles = internal_singles - 2
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)

          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 4) then
          !// {75} LOOP
          lambda_singles = internal_singles 
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)

          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix)
          permut_mat = -permut_mat
          !permut_mat = -contracted
      endif



      aa_lam_add = internal_index_vector3(lambda_path%weight +1)
      aa_mu_add  = internal_index_vector3(lambda_path%rt_loop_weight +1)
      
      if (lambda_dim /= 0.and.mu_dim /= 0) then
          

          if (aa_lam_add > 0 .and. aa_mu_add > 0) then 

              do aa1_spin  = 1,lambda_dim 
                  
                  aa_lam_start = aa_lam_add   + num_external*(aa1_spin - 1)
                  aa_lam_end   = aa_lam_start + num_external - 1 
                  
                  do aa2_spin = 1, mu_dim
                      aa_mu_start  = aa_mu_add   + num_external*(aa2_spin - 1)
                      aa_mu_end    = aa_mu_start + num_external - 1 
                                            

                      do a=1, num_external

                          one_density_matrix(i,j) = &
                               one_density_matrix(i,j) + &
                               permut_mat(aa1_spin,aa2_spin) * &
                               civec%v(aa_lam_start + a - 1) * &
                               civec%v(aa_mu_start + a - 1) 
                          
                          one_density_matrix(j,i) = &
                               one_density_matrix(j,i) + &
                               permut_mat(aa1_spin,aa2_spin) * &
                               civec%v(aa_lam_start + a - 1) * &
                               civec%v(aa_mu_start + a - 1) 
                          
                      enddo
                      
                  enddo

              enddo
               
          endif
           
      endif

  elseif (start_elec == num_elec -1) then

      !!!!!!!!!!!!!!!!!!!!
      !//                 
      !// ONE VIRTUAL SINGLY OCCUPIED             
      !//                                         
      !!!!!!!!!!!!!!!!!!!!
      
      !// GET DIMENSIONS and PERMUTATIONS

      if (loop_type == 1) then
          lambda_singles = internal_singles + 1
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 2 .or. loop_type == 3) then
          lambda_singles = internal_singles + 1
          mu_singles = internal_singles - 1
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 4) then
          lambda_singles = internal_singles + 1
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix)
          permut_mat = -permut_mat
          !permut_mat = -contracted
      endif



      a1_address = internal_index_vector1(lambda_path%weight +1)
      a2_address = internal_index_vector1(lambda_path%rt_loop_weight +1)
      
      !// MAKE SURE DIMENSIONS ARE NON ZERO    
      if (lambda_dim /= 0.and.mu_dim /= 0) then
          
          if (a1_address > 0 .and. a2_address > 0 ) then 
              

              do a1_spin  = 1,lambda_dim 
                  
                  a1_start = a1_address + num_external*(a1_spin - 1)
                  a1_end   = a1_start   + num_external - 1 
                  
                  do a2_spin = 1, mu_dim
                      
                      a2_start = a2_address + num_external*(a2_spin - 1)
                      a2_end   = a2_start   + num_external - 1 
                      

                      do a=1, num_external
                          
                          one_density_matrix(i,j) = &
                               one_density_matrix(i,j) + &
                               permut_mat(a1_spin,a2_spin) * &
                               civec%v(a1_start + a - 1) * &
                               civec%v(a2_start + a - 1) 
                          
                          one_density_matrix(j,i) = &
                               one_density_matrix(j,i) + &
                               permut_mat(a1_spin,a2_spin) * &
                               civec%v(a1_start + a - 1) * &
                               civec%v(a2_start + a - 1) 
                      enddo
                  enddo
              enddo
          endif
           
      endif



      !// HERE'S OUR LAST CASE... NO VIRTUALS OCCUPIED.  THIS IS WHERE WE LOOK FOR LOOPS 
      !// WITHIN OUR REFERENCE STATES AND INTERNAL EXCITATIONS. 
       
       
  elseif (start_elec == num_elec) then
       
      !!!!!!!!!!!!!!!!!!
      !//               
      !// VALENCE STATE; NO EXTERNALS OCCUPIED        
      !//                                             
      !!!!!!!!!!!!!!!!!!
  
  
      !// GET DIMENSIONS AND PERMUTATIONS

      if (loop_type == 1) then
          lambda_singles = internal_singles 
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 2 .or. loop_type == 3) then
          lambda_singles = internal_singles 
          mu_singles = internal_singles - 2
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix) 
          !permut_mat = contracted
      elseif (loop_type == 4) then
          lambda_singles = internal_singles
          mu_singles = lambda_singles 
          lambda_dim = fsn(lambda_singles)
          mu_dim = fsn(mu_singles)                  
          
          !permut_mat = zero
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,permut_mat,cycles,sft,spin_matrix)
          permut_mat = -permut_mat
          !permut_mat = -contracted
      endif

      
      v1_address = internal_index_vector0(lambda_path%weight+1)
      v2_address = internal_index_vector0(lambda_path%rt_loop_weight+1)
      
      !// MAKE SURE DIMENSIONS ARE NON ZERO    
      if (lambda_dim /= 0 .and. mu_dim /= 0) then
          
          if (v1_address > 0 .and. v2_address > 0 ) then 
          
              do v1_spin  = 1,lambda_dim 
                  v1_start = v1_address + v1_spin - 1
                  v1_end   = v1_start
                  
                  do v2_spin = 1, mu_dim
                      v2_start = v2_address + v1_spin - 1
                      v2_end   = v2_start
                  

                      one_density_matrix(i,j) = &
                           one_density_matrix(i,j) + &
                           permut_mat(v1_spin,v2_spin) * &
                           civec%v(v1_start) * &
                           civec%v(v2_start) 
                      
                      one_density_matrix(j,i) = &
                           one_density_matrix(j,i) + &
                           permut_mat(v1_spin,v2_spin) * &
                           civec%v(v1_start) * &
                           civec%v(v2_start) 
                      
                  enddo
              enddo
          endif
      endif
  endif
  
end subroutine one_part_purely_internal_compl
!*****************************************************************
  
end module natural_mod
  



