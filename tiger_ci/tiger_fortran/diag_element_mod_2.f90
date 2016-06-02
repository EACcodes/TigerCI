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
  !>  \brief THIS MODULE IS RESPONSIBLE FOR COMPUTING THE DIAGONAL
  !>  ELEMENTS OF THE HAMILTONIAN.
  !>  WHEN I SAY DIAGONAL ELEMENTS, I MEAN
  !>  TRULY DIAGONAL ELEMENTS AND ALSO WE COMPUTE THE OFF DIAGONAL
  !>  ELEMENTS CORRESPONDING TO DIFFERENT SPIN FUNCTIONS FOR THE
  !>  SAME SPATIAL CONFIGURATION.  THESE ELEMENTS WILL BE USED FOR
  !>  DAVIDSON PRECONDITIONING AND FORMATION OF THE SIGMA VECTOR HC.
  !>
  !>  FIRST REVISION: 9/1/99 - REWRITTEN TO USE SINGLE TREESEARCH MODULE
  !>  AND DERIVED ORBITAL_PATH DATA TYPE.
  !>
  !>   \date 1999
  !>  \author Derek Walter
!**************************************************************

! Note to reduce the size of symbol names external is abbreviated ext
! add truncation is abbreviated trunc
! David Krisiloff 12/22/2010
!


module diag_element_mod
  
  use global_var_mod
  use molecule_var_mod
  use locist_mod
  use new_tree_search_mod
  use utilities_mod
  use two_electron_integrals
  implicit none


  integer::test_count=0
  integer, dimension(:), allocatable::head_occupations
  integer, dimension(:), allocatable::tail_occupations 
  integer, dimension(:), allocatable::E_vector
  
  real(real8),allocatable,dimension(:)::diag_head     !// MATRIX TO STORE PARTIAL HAMILTONIAN MATRIX ELEMENTS
  real(real8),allocatable,dimension(:)::diag_tail     !// MATRIX TO STORE PARTIAL HAMILTONIAN MATRIX ELEMENTS
  
  real(real8), dimension(:), allocatable::diag_ab,full_matrix_ab !// DIAGONAL_AB AND FULL MATRIX ELEMENTS 
  real(real8), dimension(:), allocatable::diag_aa,full_matrix_aa !// DIAGONAL_AA AND FULL MATRIX ELEMENTS 
  real(real8), dimension(:), allocatable::diag_a,full_matrix_a   !// DIAGONAL_A  AND FULL MATRIX ELEMENTS 
  real(real8), dimension(:),allocatable::Hab,Haa,Ha1,Valence
  
  real(real8), dimension(:), allocatable::diagonal_elements      !// DIAGONAL ELEMENTS OF HAMILTONIAN MATRIX

  
  contains
      
subroutine diag_element(cho_data, loc_scr)
  
  !// THIS SUBROUTINE IS THE MAIN DRIVER ROUTINE FOR COMPUTATION OF THE DIAGONAL ELEMENTS
  !// CALLED BY:                CALLS TO:
  !//    DIAGONALIZE               GETPP
  !//                              GETPPPP
  !//                              GETPPRR
  !//                              GETPRPR  
  !//                              OPTIMAL_LEVEL
  !//                              TREE_SEARCH 
  !// 
  !//   FLOW CHART FOR THIS MODULE:
  !//                      DIAG_ELEMENT
  !//       ___________________|_____________
  !//       |              |                |                                       
  !//    GETPP        OPTIMAL_LEVEL      TREE_SEARCH               
  !//    GETPPPP                            |                    
  !//    GETPPRR                      BUILD_HEAD_ELEMENT  
  !//    GETPRPR                            |
  !//                                   TAIL_PATH
  !//                                       |
  !//                                   TREE_SEARCH
  !//                                       |
  !//                                   BUILD_TAIL_ELEMENT
  !//                                       |
  !//                                    JOIN_EM_UP
  
  
  use integral_storage_mod                 !// MODULE TO STORE INTEGRALS
  use get_integrals_mod                    !// MODULE FOR OBTAINING INTEGRALS
  use graph_var_mod                        !// GRAPH VARIABLES
  use ci_utilities_mod     
  use cholesky_structs
  use locist_var_mod,only:locist_scratch

 
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::diag_path
  type(locist_scratch)::loc_scr
  
  integer::allocatestatus          !// FOR DYNAMIC MEMORY
  integer::deallocatestatus        !// ALSO FOR DYNAMIC MEMORY
  integer::i,j,l                 !// LOOP CONTROL VARIABLES
  integer::matrix_size             !// SIZE OF A MATRIX
  
  real(real8),parameter::zero = real(0.0,real8)                !// 0.0
  type(cholesky_data)::cho_data
  type(graph_search_state) :: G
  integer :: tmp 
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// WE NEED TO ALLOCATE THE ARRAYS WHICH MAKE UP THE DERIVED TYPE
  allocate(diag_path%occupations(0:num_orbitals), stat = allocatestatus)
  call allocatecheck(allocatestatus,"%occupat")
  
  allocate(diag_path%arc_weights(0:num_orbitals), stat = allocatestatus)
  call allocatecheck(allocatestatus,"%arc_wei")
  
  allocate(diag_path%singles(0:num_orbitals), stat = allocatestatus)
  call allocatecheck(allocatestatus,"%singles")
  
  tmp = num_orbitals
  
  !// ALLOCATE THE MODULE WIDE STORAGE ARRAYS
  allocate(head_occupations(num_orbitals + 2), stat = allocatestatus)
  call allocatecheck(allocatestatus,"head_occ")
  
  allocate(tail_occupations(num_orbitals + 2), stat = allocatestatus)
  call allocatecheck(allocatestatus,"tail_occ")
  
  !// NOW LETS ALLOCATE ARRAYS FOR STORING MATRIX ELEMENTS
  matrix_size = fsn(open_shells)
  matrix_size = matrix_size*(matrix_size + 1)/2
  allocate(diag_head(matrix_size), stat = allocatestatus)
  call allocatecheck(allocatestatus,"diag_hea")
  
  allocate(E_vector(fsn(open_shells)),stat = allocatestatus)
  call allocatecheck(allocatestatus,"E_vector")
  
  allocate(diag_tail(matrix_size), stat = allocatestatus)
  call allocatecheck(allocatestatus,"diag_tai")
  
  !// NOW ALLOCATE THE ARRAY TO STORE THE DIAGONAL ELEMENTS, WHICH WE WILL
  !// KEEP IN CORE
  if (.not.allocated(diagonal_elements)) then
      allocate(diagonal_elements(total_csfs), stat = allocatestatus)
      call allocatecheck(allocatestatus, "diagonal")
  endif    
 
  write(6,*)"THE MATRIX SIZE IS", matrix_size
  
  test_count = 0

!  intdebug = 0.0D0
  
  diagonal_elements = zero  !INITIALIZE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// NOW WE NEED TO GET THE INTEGRALS WE ARE GOING TO NEED.  WE ARE
  !// GOING TO NEED THE (PP), (PP|PP), (PR|PR), AND (PP|RR).  I HAVE SEPARATE
  !// ROUTINES WHICH ARE RESPONSIBLE FOR GETTING THE INTEGRALS.  THIS WILL MAKE
  !// IT EASIER TO PORT THE CODE FOR USE WITH DIFFERENT PACKAGES THAT HAVE
  !// DIFFERENT ROUTINES FOR GENERATING INTEGRALS

  
  call getpp
 
!*****************************************
! Diagonal integrals from Cholesky Vectors
  call getpppp_cho(cho_data)
  call getpprr_cho(cho_data)
  call getprpr_cho(cho_data)
!******************************************

  !!///ARUN .. 

  allocate(Haa(matrix_size), stat = allocatestatus)
  call allocatecheck(allocatestatus,"Haa     ")
  Haa = zero 

  allocate(Hab(matrix_size), stat = allocatestatus)
  call allocatecheck(allocatestatus,"Hab     ")
  Hab = zero 

  allocate(Ha1(matrix_size), stat = allocatestatus)
  call allocatecheck(allocatestatus,"Ha1     ")
  Ha1 = zero 

  allocate(Valence(matrix_size), stat = allocatestatus)
  call allocatecheck(allocatestatus,"Valence ")
  Valence = zero  

  write(6,*) "Allocated matrix"
  !// REWIND THE DATA FILE
  rewind iodiag1 

  !// LOOP OVER THE VERTICES AT THE OPTIMAL LEVEL COMPUTING THE OCCUPATION PATTERNS
  !// AND PARTIAL MATRIX ELEMENTS.  ONCE THIS IS DONE, JOIN THEM
  do i = num_elec-2,num_elec
  
     !// LETS ZERO OUT THE APPROPRIATE ARRAYS AND INITIALIZE VARIABLES
     diag_path%occupations = 0
     diag_path%arc_weights = 0
     diag_path%singles = 0
     diag_path%num_singles = 0
     diag_path%weight = 0
      
      !// COMPUTE ALL THE PARTIAL MATRIX ELEMENTS   
       call init_tree_search(G, 0, 0, num_internal, i)
       do while ( get_internal_path(diag_path, G)) 
            call build_head_element(diag_path,loc_scr)
       end do 
  enddo
  
  !// O.K., LETS COMPUTE THE AVERAGE DIAGONAL ELEMENT.  THIS IS A GOOD WAY TO COMPARE
  !// THE RESULTS OF THIS CODE WITH THE RESULTS OF OTHER CODES (LIKE GUGA CODES) BECAUSE
  !// THE TRACE OF THE HAMILTONIAN IS INVARIANT WITH RESPECT TO UNITARY TRANSFORMATIONS
  
  write(ioOutput,*) "total_csfs diag", total_csfs
  write(ioOutput,500) sum(diagonal_elements)/real(total_csfs, real8)
  write(ioOutput,600) minval(diagonal_elements)
  call flush(ioOutput)
  500 format (1x, "Average Diagonal Element................",f25.13)
  600 format (1x, "Lowest  Diagonal Element ...............",f25.13)
  
  !// WRITE OUT DIAGONAL ELEMENTS IF REQUESTED
  if (diagprint > 0) then
      write(ioOutput,*) "Diagonal elments printed below."
      j = 1
      do i = 1, total_csfs/10
         write(ioOutput, 510) (diagonal_elements(l), l = j, j+9)
         j = j + 10
      enddo
      write(ioOutput,510) (diagonal_elements(l), l = j,total_csfs)
      510 format(1x,10(f10.3))
      call flush(ioOutput)
  endif

!  write(6,*) "intdebug",intdebug
  
  !// CLEANUP
  deallocate(diag_path%occupations, stat = deallocatestatus)
  call deallocatecheck(deallocatestatus,"%occupat")
  
  deallocate(diag_path%arc_weights, stat = deallocatestatus)
  call deallocatecheck(deallocatestatus,"%arc_wei")
  
  deallocate(diag_path%singles, stat = deallocatestatus)
  call deallocatecheck(deallocatestatus,"%singles")
  
  deallocate(E_vector, stat = deallocatestatus)
  call deallocatecheck(deallocatestatus,"E_vector")
  
  deallocate(Haa, Hab, Ha1, Valence, stat = deallocatestatus)
  call deallocatecheck(deallocatestatus,"Hab_aa  ")

end subroutine diag_element
  
  
!**************************************************************
subroutine build_head_element(diag_path,loc_scr)
  
  !// WE NEED TO CALCULATE THE SUMS:
  !// NP(PP) + D(NP,2)(PP|PP) + SUM NP*NR*((PP|RR) - (1/2)(1-D(NP*NR,1))(PR|PR))
  !//                            R
  !// WE WILL STORE THE RESULT IN ELEMENT1 OF THE ORBITAL_PATH DERIVED TYPE
  
  use graph_var_mod
  use ci_utilities_mod
  use integral_storage_mod
  use utilities_mod
  use locist_var_mod,only:locist_scratch
 
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path)::diag_path
  type(locist_scratch)::loc_scr
  
  integer::n1,n2              !// KEEPS TRACK OF NUMBER OF SINGLES, DOUBLES
  integer::r                  !// LOOP VARIABLE
  integer::singles, doubles   !// MORE LOOP VARIABLES,MORE DESCRIPTIVE
  integer::p                  !// USED TO HOLD AN ORBITAL LABEL
  integer::number_singles     !// NUMBER OF SINGLES IN THE HEAD PATH
  integer::internal_elec      !// NUMBER OF ELECTRONS IN INTERNALS SPACE
  
  integer::ab_address,aa_address,a_address,v_address
  integer::i,j,count
  integer::dim0,dim1,dim2             !// SIZE OF A MATRIX
  integer::internal_singles
  integer::internal_weight
  integer::interaction_index
!  integer::allocatestatus
!  integer::deallocatestatus        !// ALSO FOR DYNAMIC MEMORY
  
  real(real8),parameter::zero = real(0.0,real8),two = real(2.0,real8)       !// 0.0 AND 2.0
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  n1 = 0
  n2 = 0
  number_singles = 0
  diag_path%element1 = zero 
  internal_elec = sum(diag_path%occupations(0:num_internal))
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// FIRST GET NUMBER OF SINGLES IN THE HEAD PATH
  do r = 1,num_internal  
      if (diag_path%occupations(r) == 1) number_singles = number_singles + 1
  enddo
  
  head_occupations(1) = number_singles 
  levels: do r = 1,num_internal  

      if (diag_path%occupations(r) == 0) then
          
          cycle levels
          
      elseif(diag_path%occupations(r) == 1) then
          
          n1 = n1 + 1
          head_occupations(n1+1) = r
          diag_path%element1 = diag_path%element1 + pp(index2(r,r))
          
          
          !// NOW LOOP OVER SINGLES P < R
          do singles = 1, n1-1
              p = head_occupations(1+singles) 
              diag_path%element1 = diag_path%element1 + pprr(index2m1(r,p))
          enddo
          
          !// NOW LOOP OVER DOUBLES
          do doubles = 1, n2
              p = head_occupations(number_singles + 2 + doubles)
              diag_path%element1 = diag_path%element1 + two*pprr(index2m1(r,p)) - prpr(index2m1(r,p))
          enddo
          
          cycle levels
          
      elseif(diag_path%occupations(r) == 2) then
          
          n2 = n2 + 1
          head_occupations(number_singles + 2 + n2) = r
          diag_path%element1 = diag_path%element1 + two*pp(index2(r,r)) + pppp(r)

          !// LOOP OVER SINGLES
          do singles = 1,n1
              p = head_occupations(1 + singles)
              diag_path%element1 = diag_path%element1 + two*pprr(index2m1(r,p)) - prpr(index2m1(r,p))
          enddo
          
          !// LOOP OVER DOUBLES
          do doubles = 1,n2-1
              p = head_occupations(number_singles + 2 + doubles)
              diag_path%element1 = diag_path%element1 + two*(two*(pprr(index2m1(r,p)))-prpr(index2m1(r,p)))
          enddo
          
          cycle levels
          
      endif
      
  enddo levels

!  intdebug = intdebug + diag_path%element1

 
  !// LETS RECAP WHAT WE JUST DID HERE.  FOR EACH PARTIAL PATH WE MADE AN ARRAY.  THE FIRST
  !// ELEMENT IN THE ARRAY IS THE NUMBER OF SINGLES IN THE PATH.  FOR EACH SINGLE, THERE 
  !// FOLLOWS AN ELEMENT GIVING THE ORBITAL LABEL OF THE FIRST SINGLE, SECOND SINGLE, ETC.  
  !// THE NEXT ELEMENT IS GOING TO BE THE NUMBER OF DOUBLES FOLLOWED BY AN ELEMENT FOR EACH
  !// DOUBLE GIVING THE LABEL OF EACH ORBITAL.  
  
  !// LETS PUT IN THE NUMBER OF DOUBLES
  head_occupations(number_singles + 2) = n2
  
  internal_weight = sum(diag_path%arc_weights(0:num_internal))+1
  
  !// LETS CALL THIS ROUTINE ELSE WE ARE DONE 
  
  if (internal_elec == num_elec-2) then

!      if (virtual_truncation_flag /=1) then 
!          call ext_diag_nm2_aa(diag_path)
!          call ext_diag_nm2_ab(diag_path)
!      else if (virtual_truncation_flag ==1) then
!          if (direct_ci_flag == 933 .and. local_ortho_mos == 0) then    
!          call ext_diag_nm2_aa_trunc_direct(diag_path)
!          call ext_diag_nm2_ab_trunc_direct(diag_path)
!          else if (direct_ci_flag == 933 .and. local_ortho_mos == 1) then
          call ext_diag_nm2_aa_trunc_di_lmo(diag_path)
          call ext_diag_nm2_ab_trunc_di_lmo(diag_path,loc_scr)
 !         else    
 !         call ext_diag_nm2_aa_truncation(diag_path)
 !         call ext_diag_nm2_ab_truncation(diag_path)
 !         endif
 !     endif
      
      !// BY THIS TIME I HAVE Hab and Haa for EACH INTERNAL CSF
      
      internal_weight = sum(diag_path%arc_weights(0:num_internal))+1
      if (skip_this_internal(internal_weight-1,num_elec-2)) then 
      else 
          
          interaction_index = 2 
          
          !// GET THE DIMENSION OF SPIN SPACE.
          internal_singles = diag_path%num_singles  
          dim2 = fsn(internal_singles+2)
          
          !// GET THE AB_ADDRESS AND AA_ADDRESS
          ab_address = internal_index_vector2(internal_weight)
          aa_address = internal_index_vector3(internal_weight)
          
          !// WRITE OUT ALL THE INFO NECESSARY TO DATA FILE
          if (ab_address > 0) then 
              count = 0 
              do i=1,dim2
                  do j=1,i 
                      count = count + 1 
                      write(iodiag1) interaction_index
                      write(iodiag1) ab_address,aa_address,internal_weight,i,j, &
                           Hab(count),Haa(count),internal_singles

!                      intdebug = intdebug + Hab(count) + Haa(count)
                  enddo
              enddo
              
          endif
          
      endif
      
      
  elseif (internal_elec == num_elec-1) then
 !     if (virtual_truncation_flag /=1) then 
 !         call ext_diag_nm1(diag_path)
 !     else if (virtual_truncation_flag == 1) then
 !         if (direct_ci_flag == 933 .and. local_ortho_mos == 0) then    
!          call ext_diag_nm1_truncation_direct(diag_path)
 !         else if (direct_ci_flag == 933 .and. local_ortho_mos == 1) then 
          call ext_diag_nm1_trunc_di_lmo(diag_path,loc_scr)
 !         else
 !         call ext_diag_nm1_truncation(diag_path)
 !         endif
 !     endif
      
      !// BY THIS TIME I HAVE Ha for EACH INTERNAL CSF
      internal_weight = sum(diag_path%arc_weights(0:num_internal))+1
      if (skip_this_internal(internal_weight-1,num_elec-1)) then 
      else 
          
          interaction_index = 1 
          
          !// GET THE DIMENSION OF SPIN SPACE.
          internal_singles = diag_path%num_singles  
          dim1 = fsn(internal_singles+1)
          
          !// GET THE A_ADDRESS
          a_address = internal_index_vector1(internal_weight)
          
          !// WRITE OUT ALL THE INFO NECESSARY TO DATA FILE
          if (a_address > 0) then 
              count = 0 
              do i=1,dim1
                  do j=1,i 
                      count = count + 1 
                      write(iodiag1) interaction_index
                      write(iodiag1) a_address,internal_weight,i,j,Ha1(count),internal_singles

                  enddo
              enddo
          endif
      endif
      
  elseif (internal_elec == num_elec) then
      
      call external_diag_valence(diag_path)
      
      !// BY THIS TIME I HAVE valence for EACH INTERNAL CSF
      internal_weight = sum(diag_path%arc_weights(0:num_internal))+1
      if (skip_this_internal(internal_weight-1,num_elec)) then 
      else 
          
          interaction_index = 0 
          
          !// GET THE DIMENSION OF SPIN SPACE.
          internal_singles = diag_path%num_singles  
          dim0 = fsn(internal_singles)
          
          !// GET THE V_ADDRESS 
          v_address = internal_index_vector0(internal_weight)
          
          !// WRITE OUT ALL THE INFO NECESSARY TO DATA FILE
          if (v_address > 0) then 
              count = 0 
              do i=1,dim0
                  do j=1,i 
                      count = count + 1 
                      write(iodiag1) interaction_index
                      write(iodiag1) v_address,internal_weight,i,j,Valence(count),internal_singles
                     
                  enddo
              enddo
          endif
      endif
  endif
  
end subroutine build_head_element
  
         
!**************************************************************
subroutine external_diag_valence(diag_path)
  
  !// WE NEED TO CALCULATE THE SUMS:
  !// NP(PP) + D(NP,2)(PP|PP) + SUM NP*NR*((PP|RR) - (1/2)(1-D(NP*NR,1))(PR|PR))
  !//                            R
  !// HERE WE WILL STORE THE RESULT IN ELEMENT2 OF THE ORBITAL_PATH DERIVED TYPE.  
  !// THIS ROUTINE IS FOR VALENCE STATES.  
  
  use graph_var_mod
  use spin_var_mod
  use ci_utilities_mod
  use integral_storage_mod
  use utilities_mod
!  use integral_transform_mod,ONLY: ijij_int
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path),intent(in)::diag_path
  
  integer::r                   !// LOOP VARIABLE
  integer::p                   !// USED TO HOLD AN ORBITAL LABEL
  integer::spin_dim            !// DIMENSION OF THE SPIN SPACE
  integer::internal_weight     !// WEIGHT OF INTERNAL PART OF PATH
  integer::internal_singles    !// NUMBER OF SINGLES IN THE INTERNAL SPACE
  
  integer::v_address           !// ADDRESS OF VALENCE STATES
  integer::pbar,rbar           !// THE INDICES USED IN REPRESENTATIONS MATRICES US(pbar,rbar)
  integer::sf1,sf2             !// INDEXING SPIN FUNCTIONS
  integer::head_singles        !// HEAD SINGLES
  integer::v_start_1,v_end_1   !// STARTING AND ENDING POSITIONS FOR FIRST SPIN FUNCTION
  integer::v_start_2,v_end_2   !// STARTING AND ENDING POSITIONS FOR SECOND SPIN FUNCTION  
  integer::transposition_label !// TRANSPOSITION LABEL
  integer::integral_label      !// INTEGRAL LABEL INDEX
  integer::row_offset          !// ROW OFFSET FOR TRANSPOSITIONS  

  integer::count
  
  real(real8)::full_val        !// FULL ELEMENT FOR VALENCE STATE 
  real(real8)::transposition_element  !// TRANSPOSITION ELEMENTS
  real(real8),parameter::zero = real(0.0,real8),two = real(2.0,real8)              !// NUMBERS 0.0 AND 2.0
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  head_singles = head_occupations(1)
  
  internal_singles = diag_path%num_singles
  internal_weight = sum(diag_path%arc_weights(0:num_internal))+1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// CHECK TO SEE IF WE SHOULD SKIP THIS INTERNAL......
  if (skip_this_internal(internal_weight-1,num_elec)) return  
  
  v_address = internal_index_vector0(internal_weight)
  
  if (v_address < 0) return 
  
  !!!!!!!!!!!!!!!!!!!
  !//
  !// VALENCE STATE 
  !//
  !!!!!!!!!!!!!!!!!!!
  
  spin_dim = fsn(internal_singles)

  count = 0 
  
  !// NO PURELY EXTERNAL PART TO THE MATRIX ELEMENT IN THIS CASE,
  !// ALSO, NO INTERNAL/EXTERNAL CROSS TERMS
  do sf1 = 1,spin_dim
  
      v_start_1 = v_address + (sf1-1)
      v_end_1   = v_start_1 
     
      row_offset = sf1*(sf1-1)/2
      
      do sf2 = 1,sf1 
      
          v_start_2 = v_address + (sf2-1)
          v_end_2   = v_start_2
          
          !// ADD THE RELEVANT CONTRIBUTIONS FOR CASES WHERE WE ARE
          !// DIAGONAL IN THE SPIN SPACE         
          
          !// ADD IN THE PURELY INTERNAL PART OF THE TRANSPOSITION
          !// SUM.  WE DO THIS HERE, RATHER THAN IN THE CALLING
          !// ROUTINE BECAUSE WE HAVE ARRAYS ALLOCATED
     
          full_val = zero      
          do rbar = 2, head_singles
             r = head_occupations(1 + rbar)
             do pbar = 1,rbar-1
                p = head_occupations(1 + pbar)
                transposition_label = index2m1(rbar, pbar) 
                transposition_element = transpositions(transposition_label,row_offset+sf2)

                     integral_label = index2m1(r,p)
                     full_val = full_val +&
                             transposition_element*prpr(integral_label)

             enddo
          enddo
          
          !// NO PURELY EXTERNAL TERM FOR THIS CASE 
          
          !// THERE IS NO ADD INTERNAL/EXTERNAL CROSS TERMS FOR 
          !// VALENCE STATES
         
          if (sf1 == sf2) full_val =  full_val + diag_path%element1
  
          !// STORE FULL DIAGONAL ELEMENTS FOR DAVIDSON
          !// PRECONDITIONING
          if (sf1 == sf2) then
              diagonal_elements(v_start_1) = full_val
          endif
      
          !// NOW, WRITE MATRIX ELEMENTS TO DISK
          !// CODE TO WRITE GOES HERE....
          !if (virtual_truncation_flag /=1) then 
          !   write(iodiag) v_address,sf1,sf2,0
          !   write(iodiag) full_val,v_start_1,v_start_2  
          !endif 
 
          !// WE NEED THIS 
          count = count + 1 
          Valence(count) = full_val

      enddo
  
  enddo    
  
end subroutine external_diag_valence
  
!*************************************************************
Subroutine ext_diag_nm2_ab_trunc_di_lmo(diag_path,loc_scr)

  !// WE NEED TO CALCULATE THE SUMS:
  !// NP(PP) + D(NP,2)(PP|PP) + SUM NP*NR*((PP|RR) - (1/2)(1-D(NP*NR,1))(PR|PR))
  !//                            R
  !// HERE WE WILL STORE THE RESULT IN ELEMENT2 OF THE ORBITAL_PATH DERIVED TYPE.  
  !// THIS ROUTINE IS FOR STATES HAVING TWO VIRTUALS SINGLY OCCUPIED
  
  use graph_var_mod
  use spin_var_mod
  use ci_utilities_mod
  use integral_storage_mod
  use utilities_mod
  use locist_var_mod,only:locist_scratch
!  use integral_transform_mod ,ONLY: iijj_int,ijij_int,pprr_new,prpr_new,pppp_new
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path),intent(in)::diag_path
  type(locist_scratch),intent(inout)::loc_scr
  
  integer::p                   !// USED TO HOLD AN ORBITAL LABEL
  integer::r                   !// LOOP VARIABLE
  integer::spin_dim            !// DIMENSION OF THE SPIN SPACE
  integer::internal_weight     !// WEIGHT OF INTERNAL PART OF PATH
  integer::internal_singles    !// NUMBER OF SINGLES IN THE INTERNAL SPACE
  integer::allocatestatus
  integer::deallocatestatus
  
  integer::ab_address          !// ADDRESS OF AB 
  integer::a,b                 !// INDEXING EXTERNAL ORBITALS  
  integer::ab_count            !// COUNT FOR EXTERNAL ORBITALS 
  integer::pbar,rbar           !// THE INDICES USED IN REPRESENTATIONS MATRICES US(pbar,rbar)
  integer::sf1,sf2             !// INDEXING SPIN FUNCTIONS
  integer::head_singles        !// HEAD SINGLES           
  integer::ab_start_1,ab_end_1 !// STARTING AND ENDING POSITIONS FOR FIRST SPIN FUNCTION
  integer::ab_start_2,ab_end_2 !// STARTING AND ENDING POSITIONS FOR SECOND SPIN FUNCTION   
  integer::transposition_label !// TRANSPOSITION LABEL
  integer::integral_label      !// INTEGRAL LABEL INDEX
  integer::row_offset          !// ROW OFFSET FOR TRANSPOSITIONS    
  integer::count         
  integer::num_virtC2      
  integer::ione       
  integer::idum,idum2
  integer,dimension(num_orbitals)::virt_store
  
  real(real8)::transposition_element  !// TRANSPOSITION ELEMENTS
  real(real8),parameter::zero = real(0.0,real8),two = real(2.0,real8)               !// NUMBERS 0.0 AND 2.0
  real(real8)::temp_nm2_ab  

  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  head_singles = head_occupations(1)
  
  internal_singles = diag_path%num_singles
  internal_weight = sum(diag_path%arc_weights(0:num_internal))+1
 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// CHECK TO SEE IF WE SHOULD SKIP THIS INTERNAL......

  if (skip_this_internal(internal_weight-1,num_elec-2)) return  
  
  ab_address = internal_index_vector2(internal_weight)
  
  if (ab_address < 0) return 
  
  !// NOTE THE CHANGE BELOW 
  loc_scr%num_virt = num_allowed_virtuals(internal_weight,"d")
  virt_store = 0
  call get_virtuals(internal_weight,"d",virt_store)

  num_virtC2 = loc_scr%num_virt*(loc_scr%num_virt-1)/2
  ione = num_orbitals*(num_orbitals-1)/2

!  write(6,*)"I'm into ext_diag_nm2_ab_trunc_direct"

  allocate(diag_ab(num_virtC2),stat=allocatestatus)
  call allocatecheck(allocatestatus,"diag_ab ")

  allocate(full_matrix_ab(num_virtC2),stat=allocatestatus)
  call allocatecheck(allocatestatus,"full_ab ")
! NOTICE FULL_MATRIX_AB AND DIAG_AB ARE OF THE SAME DIMENSION

  !!!!!!!!!!!!!!!!!!!
  !//
  !// TWO VIRTUALS SINGLY OCCUPIED
  !//
  !!!!!!!!!!!!!!!!!!!
  
  spin_dim = fsn(internal_singles + 2)
  
  !// BUILD PURELY EXTERNAL PART OF THE MATRIX ELEMENTS
  !// WHICH DOESN'T DEPEND ON SPIN INDEX
  ab_count = 0        
  do a = 1,loc_scr%num_virt
     idum = virt_store(a)
     do b = 1,a-1
        idum2 = virt_store(b)
        ab_count = ab_count + 1
        integral_label = index2m1(idum,idum2)
  
!call check_integral(pprr(integral_label),idum,idum,idum2,idum2)

        diag_ab(ab_count) = pprr(integral_label) + pp(index2(idum,idum)) + pp(index2(idum2,idum2))
      enddo
  enddo
  
  !// SET UP DIAGONAL PARTS OF Uss(S-1,S)
  !E_vector(1:spin_dim) = E_matrix_as_vector(diag_path%num_singles + 2)

  count = 0
 
  do sf1 = 1,spin_dim
  
      ab_start_1 = ab_address + (sf1-1)*num_virtC2
      ab_end_1   = ab_start_1 + num_virtC2 - 1
     
      row_offset = sf1*(sf1-1)/2
      
      do sf2 = 1,sf1 
      
          ab_start_2 = ab_address + (sf2-1)*num_virtC2
          ab_end_2   = ab_start_2 + num_virtC2 - 1
          
          !// ADD THE RELEVANT CONTRIBUTIONS FOR CASES WHERE WE ARE
          !// DIAGONAL IN THE SPIN SPACE         
          
  
          !// ADD IN THE PURELY INTERNAL PART OF THE TRANSPOSITION
          !// SUM.  WE DO THIS HERE, RATHER THAN IN THE CALLING
          !// ROUTINE BECAUSE WE HAVE ARRAYS ALLOCATED
          
          full_matrix_ab = zero
          temp_nm2_ab = zero 
             do rbar = 2, head_singles
               r = head_occupations(1 + rbar)
                 do pbar = 1,rbar-1
                    p = head_occupations(1 + pbar)
                transposition_label = index2m1(rbar, pbar) 
                transposition_element = transpositions(transposition_label,row_offset+sf2)


                integral_label = index2m1(r,p)


!                call check_integral(prpr(integral_label),p,r,p,r)

                full_matrix_ab = full_matrix_ab +&
                             transposition_element*prpr(integral_label)
                        temp_nm2_ab  = temp_nm2_ab  + &
                             transposition_element*prpr(integral_label)
             enddo
          enddo

          !// ADD PURELY EXTERNAL TERM COMING FROM TRANSPOSITION
          transposition_label = index2m1(internal_singles+1,internal_singles+2)
          transposition_element = transpositions(transposition_label,row_offset+sf2)
          ab_count = 0
          do a = 1,loc_scr%num_virt
             idum = virt_store(a)
             do b = 1,a-1
                idum2 = virt_store(b)
                ab_count = ab_count + 1
                integral_label = index2m1(idum,idum2)


!                call check_integral( prpr(integral_label), idum , idum2, idum, idum2 )

                full_matrix_ab(ab_count) = full_matrix_ab(ab_count) + & 
                                           transposition_element*prpr(integral_label)
              enddo
          enddo            
          
          !// NOW ADD INTERNAL/EXTERNAL CROSS TERMS
          call add_ab_int_ext_cross_trunc_lmo(diag_path,row_offset,sf1,sf2,loc_scr%num_virt,virt_store(1:loc_scr%num_virt))
  
         
          !// ADD THE CONTRIBUTIONS FROM PURELY EXTERNAL AND PURELY INTERNAL PART.
          if (sf1 == sf2) full_matrix_ab = full_matrix_ab + diag_ab + diag_path%element1
          if (sf1 == sf2) temp_nm2_ab = temp_nm2_ab + diag_path%element1
          
! DIAG_PATH%ELEMENT1 FROM BUILD_HEAD_ELEMENT
          
          !// STORE FULL DIAGONAL ELEMENTS FOR DAVIDSON
          !// PRECONDITIONING
          if (sf1 == sf2) then
              diagonal_elements(ab_start_1:ab_end_1) = full_matrix_ab
          endif

          !// LETS STORE Hab 
          count = count + 1  
          Hab(count) = temp_nm2_ab 
! Hab IS OF SPIN DIMENSION          
      enddo
  
  enddo    

  deallocate(diag_ab,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"diag_ab ")

  deallocate(full_matrix_ab,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"full_ab ")

end subroutine ext_diag_nm2_ab_trunc_di_lmo
!*************************************************************
Subroutine add_ab_int_ext_cross_trunc_lmo(diag_path,row_offset,sf1,sf2,num_virt,virt_store)
  
  !// NOW WE JOIN UP THE INTERNAL AND EXTERNAL SPACE CONTRIBUTIONS.
  !// HERE THE EXTERNAL SPACE IS HAVING TWO VIRTUALS SINGLY OCCUPIED.
  !// TO JOIN THE CONTRIBUTIONS WE NEED TO ADD THE INTERNAL AND
  !// EXTERNAL SPACE DIRECTLY AND WE NEED TO DO THE INTERNAL/
  !// EXTERNAL SPACE CROSS TERMS
  
  use ci_utilities_mod    
  use graph_var_mod
  use spin_var_mod
  use integral_storage_mod
  use utilities_mod
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// ALOT OF VARIABLES HERE.  WE ARE PUTTING TOGETHER ALL THE INFORMATION FRO ALL
  !// THE PREVIOUS ROUTINES
  
  type(orbital_path),intent(in)::diag_path
  
  integer::head_singles          !// ALSO A LOOP CONTROL VARIABLE, NUMBER OF SINGLES IN HEAD PATH 
  integer::head_doubles          !// NUMBER OF DOUBLES IN A HEAD PATH
  integer::pbar                  !// ORBITAL POSITION INDICES; USED AS LOOP CONTROL VARIABLES
  integer::p                     !// ORBITAL INDICES
  integer::transposition_label   !// LABEL FOR THE TRANSOSITION WE NEED
  integer::integral_label        !// INDEX FOR THE INTEGRAL WE NEED 
  integer::internal_singles      !// NUMBER OF SINGLES IN INTERNAL SPACE
  integer::a,b                   !// INDEXING EXTERNAL ORBITALS  
  integer::ab_count              !// COUNT FOR EXTERNAL ORBITALS 
  integer,intent(in)::sf1,sf2               !// INDEXING SPIN FUNCTIONS
  integer,intent(in)::row_offset            !// ROW OFFSET FOR TRANSPOSITIONS  
  integer,intent(in)::num_virt
  integer::idum,idum2
  integer,dimension(num_orbitals),intent(in)::virt_store
  
  real(real8),parameter::zero = real(0.0, real8),two = real(2.0, real8),four  = real(4.0, real8)    !// NUMBERS 0.0,2.0 and 4.0
  real(real8)::transposition_element !// TRANSPOSITION ELEMENT 
  real(real8)::ppaa,ppbb,papa,pbpb   !// INTEGRALS  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// WE NEED TO REAL QUICKLY GET THE TAIL_SINGLES, TAIL_DOUBLES, ETC.
  head_singles = head_occupations(1)
  head_doubles = head_occupations(2 + head_singles)
  
  internal_singles = diag_path%num_singles
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !// NOW WE DO TAIL SINGLES/HEAD SINGLES FOR NP*NR(PP|RR) AND 
  !// US(PBAR,RBAR)(PR|PR)
  do pbar = 1, head_singles
      p = head_occupations(1 + pbar)
      ab_count = 0
      do a = 1,num_virt
         idum = virt_store(a)
         do b = 1,a-1
            idum2 = virt_store(b)
            ab_count = ab_count + 1 
          
              !!!!!!!
              !//
              !// FIRST DO A (A > B)
              !//
              !!!!!!!

              transposition_label = index2m1(internal_singles+2,pbar)
              transposition_element = transpositions(transposition_label,row_offset+sf2)

              integral_label = index2m1(p,idum)

!              call check_integral( prpr(integral_label), p, idum , p , idum)
!              call check_integral( pprr(integral_label), p, p, idum , idum)

              full_matrix_ab(ab_count) = full_matrix_ab(ab_count) +&
                                         transposition_element*prpr(integral_label)
              if (sf1 == sf2) full_matrix_ab(ab_count) = full_matrix_ab(ab_count) + pprr(integral_label)

              
              !!!!!!!
              !//   
              !// NOW DO B  (B < A)
              !//
              !!!!!!!

              transposition_label = index2m1(internal_singles+1,pbar)
              transposition_element = transpositions(transposition_label,row_offset+sf2)


! THIS IS THE PART YOU CAN DO SCREENING.... IF BELOW THRESHOLD ... SKIP...

              integral_label = index2m1(p,idum2)
            
!              call check_integral( prpr(integral_label), p , idum2, p, idum2)
!              call check_integral( pprr(integral_label), p, p, idum2, idum2)


              full_matrix_ab(ab_count) = full_matrix_ab(ab_count) +&
                                         transposition_element*prpr(integral_label)
                 
              if (sf1 == sf2) full_matrix_ab(ab_count) = full_matrix_ab(ab_count) + & 
                                                         pprr(integral_label)
          enddo
      enddo           
  enddo    
  
  !// NOW DO THE TAIL SINGLES/HEAD DOUBLES 2(PP|RR) - (PR|PR)
  if (sf1 == sf2) then 
      do pbar = 1, head_doubles
          p = head_occupations(2+head_singles+pbar)
          ab_count = 0
          do a = 1,num_virt
             idum = virt_store(a)
             do b = 1,a-1
                idum2 = virt_store(b)
               ab_count = ab_count + 1 
          
                  !!!!!!!
                  !//
                  !// FIRST DO A (A > B)
                  !//
                  !!!!!!!

                  integral_label = index2m1(p,idum)
                  ppaa = pprr(integral_label)
                  papa = prpr(integral_label)

!                  call check_integral(papa, p, idum, p, idum)
!                  call check_integral(ppaa, p, p, idum, idum) 


                   full_matrix_ab(ab_count) = full_matrix_ab(ab_count) + two*ppaa - papa

              
                  !!!!!!!
                  !//   
                  !// NOW DO B  (B < A)
                  !//
                  !!!!!!!

                  integral_label = index2m1(p,idum2)
                  ppbb = pprr(integral_label)
                  pbpb = prpr(integral_label)

!                  call check_integral( pbpb, p, idum2, p, idum2)
!                  call check_integral( ppbb, p, p, idum2, idum2)

                  full_matrix_ab(ab_count) = full_matrix_ab(ab_count) + two*ppbb - pbpb

              enddo
          enddo           
      enddo
  endif            
  
end subroutine add_ab_int_ext_cross_trunc_lmo
!**************************************************
 subroutine ext_diag_nm2_aa_trunc_di_lmo(diag_path)

  !// WE NEED TO CALCULATE THE SUMS:
  !// NP(PP) + D(NP,2)(PP|PP) + SUM NP*NR*((PP|RR) - (1/2)(1-D(NP*NR,1))(PR|PR))
  !//                            R
  !// HERE WE WILL STORE THE RESULT IN ELEMENT2 OF THE ORBITAL_PATH DERIVED TYPE.  
  !// THIS ROUTINE IS FOR STATES HAVING ONE VIRTUAL DOUBLY OCCUPIED
  
  use graph_var_mod
  use spin_var_mod
  use ci_utilities_mod
  use integral_storage_mod
  use utilities_mod
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path),intent(in)::diag_path
  
  integer::p                   !// USED TO HOLD AN ORBITAL LABEL
  integer::r                   !// LOOP VARIABLE
  integer::spin_dim            !// DIMENSION OF THE SPIN SPACE
  integer::internal_weight     !// WEIGHT OF INTERNAL PART OF PATH
  integer::internal_singles    !// NUMBER OF SINGLES IN THE INTERNAL SPACE
  integer::allocatestatus
  integer::deallocatestatus
  
  integer::aa_address          !// ADDRESS OF AA
  integer::a                   !// INDEXING EXTERNAL ORBITAL  A
  integer::count             !// COUNT FOR EXTERNAL ORBITAL A
  integer::pbar,rbar           !// THE INDICES USED IN REPRESENTATIONS MATRICES US(pbar,rbar)
  integer::sf1,sf2             !// INDEXING SPIN FUNCTIONS
  integer::head_singles        !// HEAD SINGLES
  integer::aa_start_1,aa_end_1 !// STARTING AND ENDING POSITIONS FOR FIRST SPIN FUNCTION
  integer::aa_start_2,aa_end_2 !// STARTING AND ENDING POSITIONS FOR SECOND SPIN FUNCTION  
  integer::transposition_label !// TRANSPOSITION LABEL
  integer::integral_label      !// INTEGRAL LABEL INDEX
  integer::row_offset          !// ROW OFFSET FOR TRANSPOSITIONS  
  integer::num_virt,num_virtC2        
  integer::idum
  integer,dimension(num_orbitals)::virt_store
  
  real(real8)::transposition_element  !// TRANSPOSITION ELEMENTS
  real(real8),parameter::zero = real(0.0,real8),two = real(2.0,real8)               !// NUMBERS 0.0 AND 2.0
  real(real8)::temp_nm2_aa
 
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  head_singles = head_occupations(1)
  
  internal_singles = diag_path%num_singles
  internal_weight = sum(diag_path%arc_weights(0:num_internal))+1

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// CHECK TO SEE IF WE SHOULD SKIP THIS INTERNAL......

  if (skip_this_internal(internal_weight-1,num_elec-2)) return  
  
  aa_address = internal_index_vector3(internal_weight)
  
  if (aa_address < 0) return
 
  !// LETS DO THIS FOR TRUNCATION TO GET THE RIGHT INTEGRALS 
  !// COMPUTE TOTAL NUMBER OF ORBITALS .. THIS NEED NOT BE THE SAME AS IT DEPENDS ON TRUNCATION

  num_virt = num_allowed_virtuals(internal_weight,"d")
  virt_store = 0
  call  get_virtuals(internal_weight,"d",virt_store)
  num_virtC2 = num_virt*(num_virt-1)/2
  !ione = num_orbitals*(num_orbitals-1)/2

  if (num_virt <= 0) return 

  allocate(diag_aa(num_virt),stat=allocatestatus)
  call allocatecheck(allocatestatus,"diag_aa ")

  allocate(full_matrix_aa(num_virt),stat=allocatestatus)
  call allocatecheck(allocatestatus,"full_aa ")

  !!!!!!!!!!!!!!!!!!!
  !//
  !// ONE VIRTUAL DOUBLY OCCUPIED
  !//
  !!!!!!!!!!!!!!!!!!!
  
  spin_dim = fsn(internal_singles)
 
  count = 0 
  
  !// BUILD PURELY EXTERNAL PART OF THE MATRIX ELEMENTS
  !// WHICH DOESN'T DEPEND ON SPIN INDEX
  do a = 1,num_virt
     idum = virt_store(a)

!     call check_integral( pppp(idum), idum, idum, idum, idum)

     diag_aa(a) = two*pp(index2(idum,idum)) + pppp(idum)
  enddo
  

!          write(*,*) "Diagonal elements before this particular step = ", diagonal_elements
      

  do sf1 = 1,spin_dim
  
      aa_start_1 = aa_address + (sf1-1)*num_virt
      aa_end_1   = aa_start_1 + num_virt - 1
     
      row_offset = sf1*(sf1-1)/2
      
      do sf2 = 1,sf1 
      
          aa_start_2 = aa_address + (sf2-1)*num_virt
          aa_end_2   = aa_start_2 + num_virt - 1
          
          !// ADD THE RELEVANT CONTRIBUTIONS FOR CASES WHERE WE ARE
          !// DIAGONAL IN THE SPIN SPACE         
          
          !// ADD IN THE PURELY INTERNAL PART OF THE TRANSPOSITION
          !// SUM.  WE DO THIS HERE, RATHER THAN IN THE CALLING
          !// ROUTINE BECAUSE WE HAVE ARRAYS ALLOCATED
          
          full_matrix_aa = zero
          temp_nm2_aa = zero 
          
          do rbar = 2, head_singles
              r = head_occupations(1 + rbar)
              do pbar = 1,rbar-1
                  p = head_occupations(1 + pbar)
                  transposition_label = index2m1(rbar, pbar) 
                  transposition_element = transpositions(transposition_label,row_offset+sf2)
                
                  integral_label = index2m1(r,p)

!                  call check_integral(prpr(integral_label), r,p,r,p)

                  full_matrix_aa = full_matrix_aa +&
                       transposition_element*prpr(integral_label)
                  temp_nm2_aa = temp_nm2_aa +   &
                       transposition_element*prpr(integral_label)
!                  write(*,*) "Full matrix aa = ", full_matrix_aa
                  flush(6)
              enddo
          enddo

 !         write(*,*) "Full matrix aa (prepost) = ", full_matrix_aa
                    
          !// NO PURELY EXTERNAL TERM FOR USs  
          
          !// NOW ADD INTERNAL/EXTERNAL CROSS TERMS

          call add_aa_int_ext_cross_trunc_lmo(diag_path,sf1,sf2,num_virt,virt_store(1:num_virt))

          !// ADD THE CONTRIBUTIONS FROM PURELY EXTERNAL AND PURELY INTERNAL 
          if (sf1 == sf2) full_matrix_aa = full_matrix_aa + diag_aa + diag_path%element1
          if (sf1 == sf2) temp_nm2_aa = temp_nm2_aa +  diag_path%element1

!          write(*,*) "Diag aa = ", diag_aa
 !         write(*,*) "Diag path = ", diag_path%element1

  !        write(*,*) "Full matrix aa (post) = ", full_matrix_aa
          
          !// STORE FULL DIAGONAL ELEMENTS FOR DAVIDSON
          !// PRECONDITIONING
          if (sf1 == sf2) then
              diagonal_elements(aa_start_1:aa_end_1) = full_matrix_aa
          endif


!          write(*,*) "Diagonal elements at this particular step = ", diagonal_elements
      
          
          !// LETS STORE Haa 
          count = count + 1 
          Haa(count) = temp_nm2_aa 
          
      enddo
      
  enddo
  
  deallocate(diag_aa,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"diag_aa ")

  deallocate(full_matrix_aa,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"full_aa ")
  
end subroutine ext_diag_nm2_aa_trunc_di_lmo
!*************************************************************
Subroutine add_aa_int_ext_cross_trunc_lmo(diag_path,sf1,sf2,num_virt,virt_store)
  
  !// NOW WE JOIN UP THE INTERNAL AND EXTERNAL SPACE CONTRIBUTIONS.
  !// HERE THE EXTERNAL SPACE IS HAVING ONE VIRTUAL DOUBLY OCCUPIED.
  !// TO JOIN THE CONTRIBUTIONS WE NEED TO ADD THE INTERNAL AND
  !// EXTERNAL SPACE DIRECTLY AND WE NEED TO DO THE INTERNAL/
  !// EXTERNAL SPACE CROSS TERMS
  
  use ci_utilities_mod    
  use graph_var_mod
  use spin_var_mod
  use integral_storage_mod
  use utilities_mod
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// ALOT OF VARIABLES HERE.  WE ARE PUTTING TOGETHER ALL THE INFORMATION FRO ALL
  !// THE PREVIOUS ROUTINES
  
  type(orbital_path),intent(in)::diag_path
  
  integer::head_singles          !// ALSO A LOOP CONTROL VARIABLE, NUMBER OF SINGLES IN HEAD PATH 
  integer::head_doubles          !// NUMBER OF DOUBLES IN A HEAD PATH
  integer::pbar                  !// ORBITAL POSITION INDICES; USED AS LOOP CONTROL VARIABLES
  integer::p                     !// ORBITAL INDICES
  integer::integral_label        !// INDEX FOR THE INTEGRAL WE NEED 
  integer::internal_singles      !// NUMBER OF SINGLES IN INTERNAL SPACE
  integer::a                     !// INDEXING EXTERNAL ORBITAL A
  integer,intent(in)::sf1,sf2               !// INDEXING SPIN FUNCTIONS  
  integer,intent(in)::num_virt
  integer::idum
  integer,dimension(num_orbitals),intent(in)::virt_store
  real(real8)::ppaa,papa
  
  real(real8),parameter::zero = real(0.0, real8),two = real(2.0, real8),four = real(4.0, real8)         !// NUMBERS 0.0,2.0 and 4.0
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// WE NEED TO REAL QUICKLY GET THE TAIL_SINGLES, TAIL_DOUBLES, ETC.
  head_singles = head_occupations(1)
  head_doubles = head_occupations(2 + head_singles)
  
  internal_singles = diag_path%num_singles
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// NOW WE DO TAIL DOUBLES/HEAD SINGLES FOR NP*NR(PP|RR) AND 
  !// US(PBAR,RBAR)(PR|PR)
  if (sf1 == sf2) then 
  do pbar = 1, head_singles
      p = head_occupations(1 + pbar)
      do a = 1,num_virt
         idum = virt_store(a)

         integral_label = index2m1(p,idum)
         ppaa = pprr(integral_label)
         papa = prpr(integral_label)

!         call check_integral( ppaa, p, p, idum, idum)
!         call check_integral( papa, p, idum, p, idum)

         full_matrix_aa(a) = full_matrix_aa(a)+ two*ppaa - papa 
      enddo
  enddo    
  endif 

  !// NOW DO THE TAIL DOUBLES/HEAD DOUBLES 2(PP|RR) - (PR|PR)
  if (sf1 == sf2) then 
      do pbar = 1, head_doubles
          p = head_occupations(2+head_singles+pbar)
          do a = 1,num_virt
             idum = virt_store(a)

             integral_label = index2m1(p,idum)
             ppaa = pprr(integral_label)
             papa = prpr(integral_label)

!             call check_integral( ppaa, p, p, idum, idum)
!             call check_integral( papa, p, idum ,p, idum)

             full_matrix_aa(a) = full_matrix_aa(a) + four*ppaa - two*papa

          enddo           
      enddo
  endif            
  
end subroutine add_aa_int_ext_cross_trunc_lmo
!*************************************************************
subroutine ext_diag_nm1_trunc_di_lmo(diag_path,loc_scr)

  !// WE NEED TO CALCULATE THE SUMS:
  !// NP(PP) + D(NP,2)(PP|PP) + SUM NP*NR*((PP|RR) - (1/2)(1-D(NP*NR,1))(PR|PR))
  !//                            R
  !// HERE WE WILL STORE THE RESULT IN ELEMENT2 OF THE ORBITAL_PATH DERIVED TYPE.  
  !// THIS ROUTINE IS FOR STATES HAVING ONE VIRTUAL SINGLY OCCUPIED
  
  use graph_var_mod
  use spin_var_mod
  use ci_utilities_mod
  use integral_storage_mod
  use utilities_mod
  use locist_var_mod,only:locist_scratch
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(orbital_path),intent(in)::diag_path
  type(locist_scratch),intent(inout)::loc_scr
  
  integer::p                   !// USED TO HOLD AN ORBITAL LABEL
  integer::r                   !// LOOP VARIABLE
  integer::spin_dim            !// DIMENSION OF THE SPIN SPACE
  integer::internal_weight     !// WEIGHT OF INTERNAL PART OF PATH
  integer::internal_singles    !// NUMBER OF SINGLES IN THE INTERNAL SPACE
  integer::allocatestatus
  integer::deallocatestatus
  
  integer::a_address           !// ADDRESS OF A
  integer::a                   !// INDEXING EXTERNAL ORBITAL  A
  integer::count               !// COUNT 
  integer::pbar,rbar           !// THE INDICES USED IN REPRESENTATIONS MATRICES US(pbar,rbar)
  integer::sf1,sf2             !// INDEXING SPIN FUNCTIONS
  integer::head_singles        !// HEAD SINGLES
  integer::a_start_1,a_end_1   !// STARTING AND ENDING POSITIONS FOR FIRST SPIN FUNCTION
  integer::a_start_2,a_end_2   !// STARTING AND ENDING POSITIONS FOR SECOND SPIN FUNCTION  
  integer::transposition_label !// TRANSPOSITION LABEL
  integer::integral_label      !// INTEGRAL LABEL INDEX
  integer::row_offset          !// ROW OFFSET FOR TRANSPOSITIONS      
  integer::idum
  integer,dimension(num_orbitals) :: virt_store
  
  real(real8)::transposition_element  !// TRANSPOSITION ELEMENTS
  real(real8),parameter::zero = real(0.0,real8),two = real(2.0,real8)               !// NUMBERS 0.0 AND 2.0
  
  real(real8)::temp_nm1 
 
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  head_singles = head_occupations(1)
  
  internal_singles = diag_path%num_singles
  internal_weight = sum(diag_path%arc_weights(0:num_internal))+1

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// CHECK TO SEE IF WE SHOULD SKIP THIS INTERNAL......
  if (skip_this_internal(internal_weight-1,num_elec-1)) return  
  
  a_address = internal_index_vector1(internal_weight)
  
  if (a_address < 0) return 
  
  !// LETS DO THIS FOR TRUNCATION TO GET THE RIGHT INTEGRALS 
  !// COMPUTE TOTAL NUMBER OF ORBITALS .. THIS NEED NOT BE THE SAME AS IT DEPENDS ON TRUNCATION

  virt_store = 0
  loc_scr%num_virt = num_allowed_virtuals(internal_weight,"s")
  call get_virtuals(internal_weight,"s",virt_store)

  allocate(diag_a(loc_scr%num_virt),stat=allocatestatus)
  call allocatecheck(allocatestatus,"diag_a  ")

  allocate(full_matrix_a(loc_scr%num_virt),stat=allocatestatus)
  call allocatecheck(allocatestatus,"full_a  ")
  !!!!!!!!!!!!!!!!!!!
  !//
  !// ONE VIRTUAL SINGLY OCCUPIED
  !//
  !!!!!!!!!!!!!!!!!!!
  
  spin_dim = fsn(internal_singles+1)

  count = 0 
  
  !// BUILD PURELY EXTERNAL PART OF THE MATRIX ELEMENTS
  !// WHICH DOESN'T DEPEND ON SPIN INDEX
  diag_a = zero 
  do a = 1,loc_scr%num_virt
     idum = virt_store(a)
     diag_a(a) = pp(index2(idum,idum)) 
  enddo
  
  do sf1 = 1,spin_dim
  
      a_start_1 = a_address + (sf1-1)*loc_scr%num_virt
      a_end_1   = a_start_1 + loc_scr%num_virt - 1
     
      row_offset = sf1*(sf1-1)/2
      
      do sf2 = 1,sf1 
      
          a_start_2 = a_address + (sf2-1)*loc_scr%num_virt
          a_end_2   = a_start_2 + loc_scr%num_virt - 1
          
          !// ADD THE RELEVANT CONTRIBUTIONS FOR CASES WHERE WE ARE
          !// DIAGONAL IN THE SPIN SPACE         
          
          !// ADD IN THE PURELY INTERNAL PART OF THE TRANSPOSITION
          !// SUM. 
  
          full_matrix_a = zero 
          temp_nm1 = zero 
          do rbar = 2, head_singles
             r = head_occupations(1 + rbar)
             do pbar = 1,rbar-1
                p = head_occupations(1 + pbar)
                transposition_label = index2m1(rbar, pbar) 
                transposition_element = transpositions(transposition_label,row_offset+sf2)
                   
                integral_label = index2m1(r,p)
                
!                call check_integral( prpr(integral_label), r,p,r,p)

                        full_matrix_a = full_matrix_a + transposition_element*prpr(integral_label)
                        temp_nm1 = temp_nm1 + transposition_element*prpr(integral_label)

             enddo
          enddo
          
          !// ADD THE FOLLOWING CONTRIBUTIONS FOR DIAGONAL IN SPIN SPACE
          if (sf1 == sf2) temp_nm1 = temp_nm1 + diag_path%element1
          !// TEMP_NM1 CONTAINS ONLY CONRIBUTIONS NOW INVOLVING INTERNAL TERMS
          
          !// NO PURELY EXTERNAL TERM FOR THIS CASE INVOLVING TRANSPOSITIONS (USs)
          
          !// NOW ADD INTERNAL/EXTERNAL CROSS TERMS
          call add_nm1_int_ext_cross_trunc_lmo(diag_path,row_offset,sf1,sf2,loc_scr%num_virt,virt_store(1:loc_scr%num_virt) )
         
          !// ADD THE CONTRIBUTIONS FROM PURELY EXTERNAL AND PURELY PART  
          if (sf1 == sf2) full_matrix_a = full_matrix_a + diag_a + diag_path%element1 
          
          !// STORE FULL DIAGONAL ELEMENTS FOR DAVIDSON
          !// PRECONDITIONING
  
          if (sf1 == sf2) then
              diagonal_elements(a_start_1:a_end_1) = full_matrix_a
          endif

          !// COPY TEMP_NM1 TO Ha1
          count = count + 1 
          Ha1(count) = temp_nm1

!          intdebug = intdebug + temp_nm1 + sum(diag_a(:))
          
      enddo
  enddo    
  
  deallocate(diag_a,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"diag_a  ")
  
  deallocate(full_matrix_a,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"full_a  ")
    
end subroutine ext_diag_nm1_trunc_di_lmo
!**************************************************************
Subroutine add_nm1_int_ext_cross_trunc_lmo(diag_path,row_offset,sf1,sf2,num_virt,virt_store)
  
  !// NOW WE JOIN UP THE INTERNAL AND EXTERNAL SPACE CONTRIBUTIONS.
  !// HERE THE EXTERNAL SPACE IS HAVING ONE VIRTUAL SINGLY OCCUPIED.
  !// TO JOIN THE CONTRIBUTIONS WE NEED TO ADD THE INTERNAL AND
  !// EXTERNAL SPACE DIRECTLY AND WE NEED TO DO THE INTERNAL/
  !// EXTERNAL SPACE CROSS TERMS
  
  use ci_utilities_mod    
  use graph_var_mod
  use spin_var_mod
  use integral_storage_mod
  use utilities_mod
!  use integral_transform_mod ,ONLY: iijj_int,ijij_int,pprr_new,prpr_new,pppp_new
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// ALOT OF VARIABLES HERE.  WE ARE PUTTING TOGETHER ALL THE INFORMATION FRO ALL
  !// THE PREVIOUS ROUTINES
  
  type(orbital_path),intent(in)::diag_path
  
  integer::head_singles          !// ALSO A LOOP CONTROL VARIABLE, NUMBER OF SINGLES IN HEAD PATH 
  integer::head_doubles          !// NUMBER OF DOUBLES IN A HEAD PATH
  integer::pbar                  !// ORBITAL POSITION INDICES; USED AS LOOP CONTROL VARIABLES
  integer::p                     !// ORBITAL INDICES
  integer::transposition_label   !// LABEL FOR THE TRANSOSITION WE NEED
  integer::integral_label        !// INDEX FOR THE INTEGRAL WE NEED 
  integer::internal_singles      !// NUMBER OF SINGLES IN INTERNAL SPACE
  integer::a                     !// INDEXING EXTERNAL ORBITALS   
  integer,intent(in)::sf1,sf2               !// INDEXING SPIN FUNCTIONS
  integer,intent(in)::row_offset            !// ROW OFFSET FOR TRANSPOSITIONS  
  integer,intent(in)::num_virt
  integer::idum
  integer,dimension(num_orbitals),intent(in)::virt_store

  
  real(real8),parameter::zero = real(0.0, real8),two = real(2.0, real8),four = real(4.0, real8)         !// NUMBERS 0.0,2.0 and 4.0
  real(real8)::transposition_element !// TRANSPOSITION ELEMENT 
  real(real8)::ppaa,papa             !// INTEGRALS  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// WE NEED TO REAL QUICKLY GET THE TAIL_SINGLES, TAIL_DOUBLES, ETC.
  head_singles = head_occupations(1)
  head_doubles = head_occupations(2 + head_singles)
  
  internal_singles = diag_path%num_singles
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !// NOW WE DO TAIL SINGLES/HEAD SINGLES FOR NP*NR(PP|RR) AND 
  !// US(PBAR,RBAR)(PR|PR)

  
  do pbar = 1, head_singles
      p = head_occupations(1 + pbar)
      do a = 1,num_virt
         idum = virt_store(a)

         transposition_label = index2m1(internal_singles+1,pbar)
         transposition_element = transpositions(transposition_label,row_offset+sf2)

         integral_label = index2m1(idum,p)
         ppaa = pprr(integral_label)

!         call check_integral(ppaa, idum,idum,p,p)
!         call check_integral( prpr(integral_label), idum, p, idum, p)

              
         full_matrix_a(a) = full_matrix_a(a) +&
         transposition_element*prpr(integral_label)
                  
         if (sf1 == sf2) full_matrix_a(a) = full_matrix_a(a) + ppaa  

      enddo           
  enddo    

  
  !// NOW DO THE TAIL SINGLES/HEAD DOUBLES 2(PP|RR) - (PR|PR)
  if (sf1 == sf2) then 
      do pbar = 1, head_doubles
         p = head_occupations(2+head_singles+pbar)
         do a = 1,num_virt  
            idum = virt_store(a)        

            integral_label = index2m1(p,idum)
            ppaa = pprr(integral_label)
            papa = prpr(integral_label)

!            call check_integral( ppaa, p, p, idum, idum)
!            call check_integral( papa, p, idum, p, idum)

                
            full_matrix_a(a) = full_matrix_a(a) + two*ppaa - papa
          enddo           
      enddo
  endif            
  
  
end subroutine add_nm1_int_ext_cross_trunc_lmo

!*************************************************************
  !> \brief Clears the current diagonal elements
  subroutine delete_diagonal_elements()
      implicit none
      ! we don't include diagonal_elements here as we usually deallocate immediately after use 
      ! since its the size of the wavefunction
      
      call try_deallocate_int_1D(head_occupations, "head_occ")
      call try_deallocate_int_1D(tail_occupations, "tail_occ")
      call try_deallocate_int_1D(E_vector, "E_vector")
      
      call try_deallocate_real_1D(diag_head, "diag_head")
      call try_deallocate_real_1D(diag_tail, "diag_tail")
      call try_deallocate_real_1D(diag_ab, "diag_ab")
      call try_deallocate_real_1D(full_matrix_ab, "full_matrix_ab")
      call try_deallocate_real_1D(diag_aa, "diag_aa")
      call try_deallocate_real_1D(full_matrix_aa, "full_matrix_aa")
      call try_deallocate_real_1D(diag_a, "diag_a")
      call try_deallocate_real_1D(full_matrix_a, "full_matrix_a")
      call try_deallocate_real_1D(Hab, "Hab")
      call try_deallocate_real_1D(Haa, "Haa")
      call try_deallocate_real_1D(Ha1, "Ha1")
      call try_deallocate_real_1D(Valence, "Valence")
      
      rewind(iodiag1)
  end subroutine

end module diag_element_mod
  
