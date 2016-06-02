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
!!*****************************************************************
!> \brief THIS MODULE IS A REINCARNATION OF THE TWO
!> SEGMENT MODULE WHICH INCLUDES SEPARATION OF THE EXTERNAL
!> SPACE
!>
!> \date 2000-2002
!> \author Derek Walter
!> \author Arun Venkatnathan
!*****************************************************************

! In order to reduce symbol size the following abbreviations are used
! David Krisiloff 12/22/2010
! internal   --> intern
! complement --> compl
! external   --> extern
!
! I'm also dropping purely from all the function names I can


module two_seg_mod_3

  use global_var_mod
  use molecule_var_mod
  use integral_storage_mod
  use get_integrals_mod
  use ci_utilities_mod
  use decider_symbols
  use tree_search_mod
  use new_tree_search_mod
  use locist_mod
  use scepmaker_mod
  use io_unit_numbers
  use cholesky_structs
  use two_seg_var_mod
  use spin_mod
  use spin_var_mod
  use twostwo_seg_info_replacement
  use pseudo_seg_info_mod
  use sgga_helpers
  use blocked_locks_mod

  use fortran_timing
#ifdef TIGER_USE_OMP
  use omp_lib
#endif

  !debugging matmul issues
  use math_utils
  
  implicit none
  
  private
  
  public :: two_seg_driver,one_intern_seg_compl_vec_lmo,intern_two_seg_compl_vec_lmo

#ifdef TIGER_USE_OMP
  !// THESE ARRAYS ARE *ONLY* USED TO COMMUNICATE INFORMATION (OFFSETS ETC)
  !// IN THE OPENMP CASE WITHIN THIS MODULE. DO NOT USE FOR ANYTHING ELSE!
  !// ALSO, DO NOT PUT ANY OTHER GLOBAL MODULE VARIABLES HERE WITHOUT 
  !// PROPER JUSTIFICATION!
  integer,dimension(:),allocatable,private::omp_offsets_twostwo_rec
  integer,dimension(:),allocatable,private::omp_offsets_twosingles_rec
  integer,dimension(:),allocatable,private::omp_offsets_ipjp_rec
#endif
  
  integer, private, parameter, dimension(5,2,2)::full_loops=&
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
       
  type purelyInternalScratch
     integer,dimension(:),allocatable::virt_lam_allow,virt_mu_allow
  end type purelyInternalScratch

  contains


  !*****************************************************************
  subroutine two_seg_driver(civec, sigmavec, iteration,cho_data,loc_scr,twoVars)
  
    use locist_var_mod,only:locist_scratch
    use three_four_seg_var_mod
    use time_var_mod
  
    implicit none

    !// THIS SUBROUTINE IS THE MAIN DRIVER FOR TREATMENT OF THE
    !// TWO SEGMENT LOOPS WITH SEPARATION OF THE EXTERNAL SPACE.  

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(orbital_path)::lambda_path
    type(cholesky_data)::cho_data
    type(locist_scratch)::loc_scr
    type(twoModVars)::twoVars
    type(blockedLockVectorType) :: civec, sigmavec 

    integer::allocatestatus             !// FOR DYNAMIC MEMORY ALLOCATION
    integer::deallocatestatus
    integer::iteration
    type(clock) :: timer

#ifdef TIGER_USE_OMP
    integer::numthreads
    numthreads = numberOfThreads
#else
    integer,parameter::numthreads = 1
#endif
    

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(lambda_path%arc_weights(0:num_orbitals), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%2sg_wei")
    lambda_path%arc_weights = 0

    allocate(lambda_path%occupations(0:num_orbitals), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%lam_occ")
    lambda_path%occupations = 0

    allocate(lambda_path%singles(0:num_orbitals), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%lam_sin")
    lambda_path%singles = 0

    allocate(lambda_path%constraints(3,6), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%lam_con")
    lambda_path%constraints = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// GET THE PESKY ONE ELECTRON INTEGRALS
    call getpp

    !// ALLOCATION OF SCEPPER MATRICES  
    allocate(twoVars%scep_ci_1(num_external,num_external),&
         stat = allocatestatus)
    call allocatecheck(allocatestatus,"scep_ci1")             
    twoVars%scep_ci_1 = 0

    !// ALLOCATION OF SCEPPER MATRICES  
    allocate(twoVars%scep_ci_2(num_external,num_external),&
         stat = allocatestatus)
    call allocatecheck(allocatestatus,"scep_ci2")             
    twoVars%scep_ci_2 = 0

!!!!!!!!!!!!!!!!!
    !// 
    !// PURELY EXTERNAL TWO SEGMENT LOOPS
    !//
!!!!!!!!!!!!!!!!!

    twoVars%iter_count = iteration

    if (valence_ci_flag .ne. 1 .and. reference_ci_flag .ne. 1) then   
       call start_clock(timer)
       call purely_external_two_seg_cho_2(civec%v, sigmavec%v, cho_data,lambda_path,loc_scr,twoVars) ! faster
       call print_clock(timer, "  + purely external")
    endif


!!!!!!!!!!!!!!!!!
    !// 
    !// PURELY INTERNAL TWO SEGMENT LOOPS
    !//
!!!!!!!!!!!!!!!!!
    call start_clock(timer)
    call  purely_internal_two_seg_2_cho(civec%v, sigmavec%v, sigmavec%l, cho_data,lambda_path,twoVars)
    call print_clock(timer, "  + purely internal")

!!!!!!!!!!!!!!!!!
    !// 
    !// TWO SEGMENT LOOPS WITH ONE SEGMENT
    !// IN THE INTERNAL SPACE
    !//
!!!!!!!!!!!!!!!!! 

    if (valence_ci_flag .ne. 1 .and. reference_ci_flag .ne. 1) then
       call start_clock(timer)
       call one_internal_seg_2_cho(civec, sigmavec, cho_data,lambda_path,loc_scr,twoVars)
       call print_clock(timer, "  + one external")
    endif

    !// CLEAN UP
    deallocate(twoVars%scep_ci_1, stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"scep_ci1")               

    deallocate(twoVars%scep_ci_2, stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"scep_ci2")               

    deallocate(lambda_path%arc_weights, stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"%2sg_wei")

    deallocate(lambda_path%occupations, stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"%lam_occ")

    deallocate(lambda_path%singles, stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"%lam_sin")

    deallocate(lambda_path%constraints, stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"%lam_con")

  end subroutine two_seg_driver

  !*****************************************************************
  subroutine add_in_Gij(Kab, lambda_dim, Gij, loop_type)

    !// THIS SUBROUTINE JUST ADDS IN Gij TO THE DIAGONAL ELEMENTS OF THE
    !// Kab MATRIX.  

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,intent(in)::lambda_dim                             !// NUMBER OF DIAGONAL ELEMENTS
    integer::i                                      !// LOOP CONTROL VARIABLE
    integer,intent(in)::loop_type                              !// TYPE OF LOOP

    real(real8),dimension(:,:),intent(inout)::Kab         !// THE MATRIX WE ARE BULKING UP
    real(real8),intent(in)::Gij                                !// WHAT WE ADD TO THE MATRIX
    
!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (loop_type == 4) then
       do i = 1, lambda_dim
          Kab(i,i) = Kab(i,i) - Gij
       enddo
    else
       do i = 1, lambda_dim
          Kab(i,i) = Kab(i,i) + Gij
       enddo
    endif

  end subroutine add_in_Gij


  !*****************************************************************
  subroutine purely_external_two_seg_cho_2(civec, sigmavec, cho_data,lambda_path,loc_scr,twoVars)

    !// THE FIRST STEP IN BUILDING THE EXTERNAL SPACE CONTRIBUTIONS
    !// IS TO SEARCH DOWN TO THE INTERNAL EXTERNAL SPACE BOUNDARY.
    !// THUS, THIS IS A VERY SIMPLE ROUTINE

    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use three_four_seg_var_mod

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(orbital_path)::lambda_path
    integer::allocatestatus,deallocatestatus
    integer::a,b,ab_count,i
    integer::offset    
    integer::idum     
    integer::iiab_count,iaib_count
    real(kind=real8), dimension(:) :: civec, sigmavec
    
    real(kind=real8), external::ddot ! BLAS1 (dot product)
    type(cholesky_data)::cho_data
    type(locist_scratch)::loc_scr
    type(twoModVars)::twoVars

    type(graph_search_state) :: graph
#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num() + 1
#else
    integer,parameter::threadID = 1
#endif
    
    allocate(twoVars%aibi_debug(2*num_internal*num_external*(num_external+1)/2), &
       stat=allocatestatus)
    call allocatecheck(allocatestatus,"purely_external_two_seg_cho_2")
    twoVars%aibi_debug = 0


!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// INITIALIZE DERIVED TYPE ARRAYS
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%num_singles = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    iiab_count = 0
    iaib_count = 0

    do i = 1, num_internal

      !********************************************************************************
      ! For the iiab

      offset = (i-1)*num_external*(num_external+1)/2
      ab_count = 0

      do a = num_internal+1,num_orbitals
         do b = num_internal+1, a

            ab_count = ab_count + 1

            idum = a*(a-1)/2+b
            idum = cho_data%mo_ind_inv(idum)
            if (idum .eq. 0) cycle

            if (ignorable_pair(b,a) ) cycle

            iiab_count = iiab_count + 1
            
            call for_double_buf_readElement(cd_iiab_no,iiab_count,twoVars%aibi_debug(offset+ab_count),threadID)

         enddo
      enddo
      
      !*********************************************************************************
      ! For the iaia,iaib

      offset = num_internal*num_external*(num_external+1)/2 + (i-1)*num_external*(num_external+1)/2
      ab_count = 0

      do a = 1, num_external
         do b = 1,a

            ab_count = ab_count + 1

            if (ignorable_pair(a+num_internal,i) ) cycle
            if (ignorable_pair(b+num_internal,i) ) cycle
            
            iaib_count = iaib_count + 1
                
            call for_double_buf_readElement(cd_iaib_no,iaib_count,twoVars%aibi_debug(offset+ab_count),threadID)
            
         enddo
      enddo

    enddo

    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%num_singles = 0

    call init_tree_search(graph, 0, 0, num_internal, num_elec-2)
    do while ( get_internal_path(lambda_path, graph))
      call extern_two_seg_compl_vec_lmo(civec, sigmavec, lambda_path,loc_scr,twoVars)
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
      call extern_two_seg_compl_vec_lmo(civec, sigmavec, lambda_path,loc_scr,twoVars)
    end do

    deallocate(twoVars%aibi_debug,stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"purely_external_two_seg_cho_2")
    
  end subroutine purely_external_two_seg_cho_2
  !*****************************************************************
  subroutine make_gij_internal(lambda_path, level1, level2, Gij_internal, offset)

    !// IN THIS SUBROUTINE WE MAKE THE CONTRIBUTION TO Gij FROM THE
    !// PORTION OF A PATH WHICH RESIDES ENTIRELY IN THE INTERNAL SPACE.
    !// THE FORMULA IS:
    !//     Gij_internal = sum' { Np(ij|pp) - delta(Np,2)(ip|jp) }
    !// WHERE THE SUM IS OFVER ORBITALS p AND THE PRIME
    !// MEANS DON'T INCLUDE ORBITALS i AND j IN THE SUM.  Np IS THE 
    !// OCCUPATION OF ORBITAL p.  
    !//
    !// WARNING:  NOTE THAT WE DON'T INCLUDE THE ONE
    !// ELECTRON INTEGRAL HERE.  THIS IS SO THAT IF WE HAVE TO CALL THIS
    !// ROUTINE TWICE (AS IN THE CASE OF PURELY INTERNAL LOOPS) WE DON'T
    !// ADD IT IN TWICE.  ALSO, THIS ROUTINE ASSUMES THAT YOU INITIALIZE
    !// Gij_internal TO ZERO BEFORE CALLING. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(orbital_path)::lambda_path

    integer::path_level                 !// FOR LOOPING OVER LEVELS
    integer::level1, level2             !// LOOP LEVELS
    integer::offset                     !// FOR ACCESSING INTEGRALS

    real(real8)::Gij_internal           !// THIS IS WHAT WE COMPUTE
    real(real8),parameter::two=real(2.0, real8)                    !// NUMBER 2

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do path_level = 1, num_internal

       if (path_level == level1.or.path_level == level2) cycle

       if (lambda_path%occupations(path_level) == 1) then
          Gij_internal = Gij_internal + integral_buffer(1,path_level+offset)
       elseif (lambda_path%occupations(path_level) == 2) then
          Gij_internal = Gij_internal + two*integral_buffer(1,path_level+offset) &
               -     integral_buffer(2,path_level+offset)
       endif

    enddo

  end  subroutine make_gij_internal

  !****************************************************************
  subroutine compute_Kab_2(num_singles,orb_sing,level1, level2, loop_type,&
       ibar, jbar, lambda_singles, mu_singles,&
       Kab, size, offset,ipjp_rec)

    !// FOR THE {13}, {35}, {53}, AND {57} LOOPS THE MATRIX ELEMENTS ALL
    !// INVOLVE A SUM OVER SINGLY OCCUPIED ORBITALS MULTIPLIED WHERE YOU
    !// MULTIPLY AN INTEGRAL BY A TRANSPOSITION.  HERE WE TREAT THIS SUM.
    !// THE FORMULAS ARE AS FOLLOWS:
    !// 
    !//    {13} = SUM { DELTA(Np,1) (ibar, pbar) (IP|JP) }
    !//    {35} = SUM { DELTA(Np,1) (jbar, pbar) (IP|JP) }
    !//    {53} = SUM { DELTA(Np,1) (ibar, pbar) (IP|JP) }
    !//    {57} = SUM { DELTA(Np,1) ((ibar, pbar) + 1) (IP|JP) }

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    integer,intent(in)::ipjp_rec

    integer::threadID
    integer::level1                         !// STARTING LEVEL FOR SUM
    integer::level2                         !// ENDING LEVEL FOR SUM
    integer::loop_type                      !// TYPE OF THE LOOP
    integer::ibar,jbar                      !// POSITION INDICES
    integer::lambda_singles                 !// NUMBERS OF SINGLES
    integer::mu_singles
    integer::j,k                          !// USED FOR LOOPING OVER LEVELS
    integer::pbar                           !// KEEPS TRACK OF POSITION INDICES
    integer::lambda_dim, mu_dim             !// DIMENSIONS OF MATRICES
    integer::size                           !// SIZE OF TRANSPOSITION SUM
    integer::transposition_label            !// LABELS THE TRANSPOSITION SUM
    integer::row_offset
    integer::orbital
    integer::offset                         !// FOR INTEGRAL ACCESS
    integer::num_singles
    integer,dimension(num_singles):: orb_sing

    real(real8)::the_integral               !// POINTER TO A LPVOID BUNGHOLE CLASS
    real(real8)::transposition_element      !// ELEMENT OF A TRASPOSITION MATRIX
    real(real8), dimension(size,size)::Kab          !// SEE FORMULAS IN INTRO
    
#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pbar = 0
    lambda_dim = fsn(lambda_singles)
    mu_dim = fsn(mu_singles)

    Kab = real(0.0,real8)

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// LOOP OVER SINGLY OCCUPIED ORBITALS IN THE INTERNAL SPACE
    do pbar = 1, num_singles

       !// GET THE ORBITAL LABEL
       orbital = orb_sing(pbar)

       !// CYCLE IF THIS IS ONE OF THE LOOP LEVELS
       if (orbital == level1) cycle
       if (orbital == level2) cycle

       call for_double_buf_readElement(cd_ipjp_no,ipjp_rec*num_internal+offset+orbital,the_integral,threadID)
       
!!!!!!!!
       !//
       !// {35} LOOPS NEED (JBAR, PBAR)
       !//
!!!!!!!!
       if (loop_type == 2) then

          if (pbar /= jbar) then

             transposition_label = index2m1(pbar,jbar)

             do j = 1, lambda_dim
                row_offset = j*(j-1)/2
                do k = 1, j-1
                   transposition_element = &
                        transpositions(transposition_label, row_offset + k)
                   Kab(j,k) = Kab(j,k) + &
                        the_integral*transposition_element
                   Kab(k,j) = Kab(k,j) + &
                        the_integral*transposition_element                        
                enddo
                transposition_element = &
                     transpositions(transposition_label, row_offset + j)

                Kab(j,j) = Kab(j,j) + &
                     the_integral*transposition_element
             enddo

          else

             !// PBAR == JBAR => (JBAR,PBAR) = IDENTITY
             do j = 1, lambda_dim
                Kab(j,j) = Kab(j,j)+the_integral
             enddo


          endif

!!!!!!!!!
          !//    
          !// HERE WE NEED (ibar, pbar) TO TREAT EVERYTHING BUT
          !// THE {35} LOOPS.   IN THE CASE OF THE {57} LOOPS
          !// WE WILL ADD IN THE IDENTITY AFTER
          !//
!!!!!!!!!
       else


          if (pbar /= ibar) then

             transposition_label = index2m1(pbar, ibar)

             do j = 1, lambda_dim
                row_offset = j*(j-1)/2
                do k = 1, j-1
                   transposition_element = &
                        transpositions(transposition_label, row_offset + k)
                   Kab(j,k) = Kab(j,k) + &
                        the_integral*transposition_element
                   Kab(k,j) = Kab(k,j) + &
                        the_integral*transposition_element    
                enddo
                transposition_element = &
                     transpositions(transposition_label, row_offset + j)
                Kab(j,j) = Kab(j,j) + &
                     the_integral*transposition_element
             enddo

          else

             !// PBAR == IBAR => (IBAR,PBAR) = IDENTITY
             do j = 1, lambda_dim
                Kab(j,j) = Kab(j,j)+the_integral
             enddo

          endif

       endif

!!!!!!!!
       !//
       !// ADD IN THE INDENTITY FOR LOOP {57}
       !//
!!!!!!!!
       !// ADD IN THE INDENTITY FOR LOOP {57}
       if (loop_type == 4) then
          do j = 1, lambda_dim
             Kab(j,j) = Kab(j,j)+the_integral
          enddo
       endif

    enddo

  end subroutine compute_Kab_2
  !*****************************************************************

  !*****************************************************************
  subroutine one_internal_seg_2_cho(civec, sigmavec, cho_data,lambda_path,loc_scr,twoVars)

    !// MAIN DRIVER ROUTINE FOR TWO SEGMENT LOOPS WITH
    !// ONE SEGMENT IN THE INTERNAL SPACE
    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use three_four_seg_var_mod

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(blockedLockVectorType) :: civec, sigmavec
    type(orbital_path)::lambda_path
    type(cholesky_data)::cho_data
    type(locist_scratch)::loc_scr
    type(twoModVars)::twoVars

    integer::i                         !// LABELS LOOP LEVELS

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// INITIALIZE DERIVED TYPE ARRAYS
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%num_singles = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rewind (unit=335) 
    rewind (unit=336)

    !// LOOP OVER THE INTERNAL WHERE THE LOOP IS GOING TO BE FORMED
    do i = 1, num_internal

       twoVars%i_ind = i


       !// GET INTEGRALS
       call  get_ipap_cho(cho_data,i,twoVars%iter_count)

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
!       if (virtual_truncation_flag /=1) then
!          lambda_path%constraints(1,6) = 12 !// NEXT TASK
!       else if (virtual_truncation_flag ==1 .and. local_ortho_mos == 0) then 
!          lambda_path%constraints(1,6) = 42 !// NEXT TASK
!       else if (virtual_truncation_flag ==1 .and. local_ortho_mos == 1) then
          lambda_path%constraints(1,6) = LMO_1_INT_COMP!// NEXT TASK
!       endif

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
       call broken_constrained_search(lambda_path,0,0,loc_scr=loc_scr,a=civec,b=sigmavec)


!!!!!!!!!!!!!!!!!
       !//
       !// NOW THE {5} SEGMENT
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
!       if (virtual_truncation_flag /=1) then
!          lambda_path%constraints(1,6) = 12 !// NEXT TASK
!       else if (virtual_truncation_flag ==1 .and. local_ortho_mos == 0) then
!          lambda_path%constraints(1,6) = 42 !// NEXT TASK
!       else if (virtual_truncation_flag == 1 .and. local_ortho_mos == 1) then
          lambda_path%constraints(1,6) = LMO_1_INT_COMP !// NEXT TASK
!       endif

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
       call broken_constrained_search(lambda_path,0,0,loc_scr=loc_scr,a=civec,b=sigmavec)

    enddo

  end subroutine one_internal_seg_2_cho
  !*****************************************************************
  subroutine one_info_simplifier(loop_type,start_elec,internal_singles,ibar,jbar, & 
       lambda_singles,mu_singles, lambda_dim, mu_dim, & 
       lambda_weight,mu_weight,& 
       v_address,a_address,ab_address,aa_address,& 
       Px)

    !// THIS SUBROUTINE GETS THE COUPLING COEFFICIENT INFORMATION
    !// FOR THE LOOPS WITH ONE INTERNAL SEGMENT AND PACKAGES THEM
    !// FOR THE VECTORIZED MULTIPLICATION                                  

    use spin_mod,only:contracted_cycle_stateless
    use spin_var_mod,only:cycles,sft,spin_matrix

    implicit none

    !// THIS ROUTINE COMPUTES COUPLING COEFFICIENTS, DIMENSIONS OF LAMBDA, MU AND 
    !// ADDRESSES FOR VALENCE,N-1 and N-2 STATES.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,intent(in)::loop_type                                                
    integer,intent(in)::start_elec 
    integer,intent(in)::internal_singles                                          
    integer,intent(in)::ibar,jbar
    integer,intent(out)::lambda_singles,mu_singles
    integer,intent(out)::lambda_dim,mu_dim
    integer,intent(in)::lambda_weight,mu_weight
    integer,intent(out)::v_address,a_address,ab_address,aa_address

    real(real8),dimension(:,:),intent(inout)::Px 

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (start_elec == num_elec - 1) then 

       if (loop_type == 1) then 

          lambda_singles = internal_singles + 1 
          mu_singles     = lambda_singles 

          lambda_dim     = fsn(lambda_singles)
          mu_dim         = lambda_dim 

          v_address      = internal_index_vector0(mu_weight)
          a_address      = internal_index_vector1(lambda_weight)

          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)

       elseif (loop_type == 2) then 

          lambda_singles = internal_singles + 1 
          mu_singles     = internal_singles - 1 

          lambda_dim     = fsn(lambda_singles)
          mu_dim         = fsn(mu_singles)

          v_address      = internal_index_vector0(mu_weight)
          a_address      = internal_index_vector1(lambda_weight)

          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)

       endif

    elseif (start_elec == num_elec - 2) then 

       if (loop_type == 1) then 

          lambda_singles = internal_singles + 2  
          mu_singles     = lambda_singles 

          lambda_dim     = fsn(lambda_singles)
          mu_dim         = lambda_dim 

          a_address      = internal_index_vector1(mu_weight)
          ab_address     = internal_index_vector2(lambda_weight)
          aa_address     = internal_index_vector3(lambda_weight)

          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)

       elseif (loop_type == 2) then 

          lambda_singles = internal_singles + 2 
          mu_singles     = internal_singles  

          lambda_dim     = fsn(lambda_singles)
          mu_dim         = fsn(mu_singles)

          a_address      = internal_index_vector1(mu_weight)
          ab_address     = internal_index_vector2(lambda_weight)
          aa_address     = internal_index_vector3(lambda_weight)

          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)

       endif

    endif

  end subroutine one_info_simplifier

  !*****************************************************************
  !*****************************************************************

  !*****************************************************************
  subroutine add_in_Gij_vec(Ka, lambda_dim, Gij, loop_type,KaGijind)

    !// THIS SUBROUTINE JUST ADDS IN Gij TO THE DIAGONAL ELEMENTS OF THE
    !// Ka MATRIX.  

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,intent(in)::lambda_dim                             !// NUMBER OF DIAGONAL ELEMENTS
    !integer::length                                 !// SIZE OF MATRIX
    integer,intent(in)::loop_type                              !// TYPE OF LOOP
    integer,intent(in)::KaGijind
    real(real8), dimension(:,:,:),intent(inout)::Ka !// THE MATRIX WE ARE BULKING UP
    real(real8),dimension(:),intent(in)::Gij         !// WHAT WE ADD TO THE MATRIX
    
    integer::i,a                                    !// LOOP CONTROL VARIABLE

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (loop_type == 4) then
       do a = 1,KaGijind
          do i = 1, lambda_dim
             Ka(i,i,a) = Ka(i,i,a) - Gij(a)
          enddo
       enddo
    else
       do a = 1,KaGijind
          do i = 1, lambda_dim
             Ka(i,i,a) = Ka(i,i,a) + Gij(a)
          enddo
       enddo
    endif

  end subroutine add_in_Gij_vec



  !*****************************************************************
  subroutine purely_internal_two_seg_2_cho(civec, sigmavec, sigmavec_lock, cho_data,lambda_path,twoVars)

    !// IN THIS ROUTINE WE TREAT THE PURELY INTERNAL TWO SEGMENT
    !// LOOPS.  THIS WILL ALSO INCLUDE THE SPECIAL {26} LOOPS. THIS SUBROUTINE
    !// IS IN THE VECTORIZED MODE

    use get_integrals_mod,only:prpr,pppp
    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use three_four_seg_var_mod

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(orbital_path)::lambda_path
    type(cholesky_data)::cho_data
    type(twoModVars)::twoVars
    real(kind=real8), dimension(:) :: civec, sigmavec
    type(blockedLockType) :: sigmavec_lock

    integer::i,j                        !// LABELS LOOP LEVELS
    integer::loop                       !// LABELS THE LOOP
    integer::idum,idum2,idum3,ij_1
    integer::counter,twostwo_rec,twosingles_rec
    integer::allocatestatus

    type(purelyInternalScratch),dimension(:),allocatable::scr

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef TIGER_USE_OMP
    integer::numthreads
    integer::threadID
    numthreads = numberOfThreads
#else
    integer,parameter::threadID = 1
    integer,parameter::numthreads = 1
#endif
    !// INITIALIZE DERIVED TYPE ARRAYS
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%num_singles = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !// GET THE INTEGRALS FOR TREATING THE {26} LOOPS.  THIS IS A SMALL
    !// STORAGE REQUIREMENT WE HOLD THE ENTIRE THING IN CORE

    !  YOU ALSO NEED (ij|ij)

    twoVars%twostwo_rec = 0
    twoVars%twosingles_rec = 1 ! please note that this must be one!

    idum3 = 0
    
    if(twoVars%iter_count .eq. 1) then

#ifdef TIGER_USE_OMP
    ! allocate omp stuff
    if(.not.allocated(omp_offsets_ipjp_rec)) then
      allocate(omp_offsets_ipjp_rec(num_internal),omp_offsets_twostwo_rec(num_internal), &
      omp_offsets_twosingles_rec(num_internal),stat=allocatestatus)
      call allocatecheck(allocatestatus,"omp_stuff_twoseg")
    endif
#endif
    
    counter = 0
    do i = 1, num_internal
    
#ifdef TIGER_USE_OMP
       omp_offsets_ipjp_rec(i) = counter
       omp_offsets_twostwo_rec(i) = twoVars%twostwo_rec
       omp_offsets_twosingles_rec(i) = twoVars%twosingles_rec
#endif
        
       do j = 1, i-1
       
          idum = i*(i-1)/2+j
          idum2 = (i-1)*(i-2)/2+j

          if (prpr(idum2)*pppp(i) .lt. integral_threshold) cycle
          if (prpr(idum2)*pppp(j) .lt. integral_threshold) cycle 

          !// GET INTEGRALS
          counter = counter + 1
          call make_ipjp(cho_data,i,j,counter)
          
             do loop = 1,5

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
                lambda_path%constraints(1,6) = LMO_INT_2_COMP !// NEXT TASK
                lambda_path%constraints(2,3) = num_internal+1

                !// FOR BUILDING THE MU PATH
                lambda_path%level1 = j

                !// LOOP LEVELS
                lambda_path%constraints(2,1:2) = (/j,i/)

                !// LAMBDA, MU PATHS
                lambda_path%constraints(1,1:2) = &
                     full_loops(loop,1,1:2)

                lambda_path%constraints(3,1:2) = &
                     full_loops(loop,2,1:2)

                !// NOW DO THE SEARCH
                lambda_path%num_singles = 0
                call broken_constrained_search(lambda_path,0,0,twoVars=twoVars)

             enddo
       enddo
    enddo

    endif
    
    counter = 0
    twostwo_rec = 0
    twosingles_rec = 1
    
    ! allocate our scratch space
    allocate(scr(numthreads),stat=allocatestatus)
    call allocatecheck(allocatestatus,"purelyInternalScratch")
    do i = 1, numthreads
       allocate(scr(i)%virt_lam_allow(num_external),scr(i)%virt_mu_allow(num_external),stat=allocatestatus)
       call allocatecheck(allocatestatus,"purelyInternalScratchSpaces")
    enddo
    
    !// LOOP OVER PAIRS OF INTERNALS I > J
    !$omp parallel do &
    !$omp schedule(static) &
    !$omp default(none) &
    !$omp private(idum,idum2,lambda_path,ij_1,counter,twostwo_rec,twosingles_rec,threadID) &
    !$omp shared(civec,sigmavec,sigmavec_lock,num_internal,prpr,pppp,integral_threshold, &
    !$omp omp_offsets_ipjp_rec,omp_offsets_twostwo_rec, scr, &
    !$omp omp_offsets_twosingles_rec)
    do i = 1, num_internal
    
#ifdef TIGER_USE_OMP
       counter = omp_offsets_ipjp_rec(i)
       twostwo_rec = omp_offsets_twostwo_rec(i)
       twosingles_rec = omp_offsets_twosingles_rec(i)
       threadID = OMP_get_thread_num()+1
#endif
        
       do j = 1, i-1
          idum = i*(i-1)/2+j
          idum2 = (i-1)*(i-2)/2+j

          if (prpr(idum2)*pppp(i) .lt. integral_threshold) cycle
          if (prpr(idum2)*pppp(j) .lt. integral_threshold) cycle 

          counter = counter + 1

          ij_1 = idum
          do while(ij_1 == idum )
             call pure_int_2_seg_comp_vec_lmo_res(civec, sigmavec, sigmavec_lock, ij_1,twostwo_rec,twosingles_rec,scr(threadID)%virt_lam_allow,scr(threadID)%virt_mu_allow,counter)
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    ! no deallocate necessary, fortran must do it

  end subroutine purely_internal_two_seg_2_cho

  !*****************************************************************
  !*****************************************************************
  !****************************************************************
  subroutine intern_two_seg_compl_vec_lmo(lambda_path,twoVars)

    !// IN THIS ROUTINE WE BUILD THE COMPLEMENT TO THE TWO
    !// SEGMENT LOOPS WHICH RESIDE ENTIRELY IN THE INTERNAL SPACE
    !// THIS IS A LONG COMPLICATED AND ARDUOUS ROUTINE.  NOT FOR
    !// THE FAINT OF HEART.  

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!

    use spin_mod 
    !use get_integrals_mod,ONLY:prpr
    use locist_var_mod,only:locist_scratch
    
    implicit none

    type(orbital_path),intent(inout)::lambda_path
    type(twoModVars),intent(inout)::twoVars

    integer::ibar,jbar               !// ORBITAL POSITION INDICES
    integer::internal_singles        !// NUMBER OF SINGLES IN INTERNAL PATH OF PATH
    integer::current_vertex          !// USED IN BUILDING MU PATH
    integer::step_type               !// KEEPS TRACK OF STEP IN MU PATH CONSTRUCTION
    integer::path_elecs              !// CUMULATIVE OCCUPATION IN MU PATH
    integer::top_level               !// LEVEL OF FIRST LOOP SEGMENT
    integer::mu_levels               !// LABELS LEVELS IN MU PATH
    integer::constraint_count        !// KEEP TRACK OF CONSTRAINTS IN BUILDING MU PATH
    integer::start_elec              !// OCCUPATION WHERE INTERNAL CSF MEETS UP WITH EXTERNAL SPACE
    integer::s1,s2,s3,s4             !// FOR THE CALL TO GET NUMBERS OF SINGLES
    integer::loop_type               !// LABELS THE LOOP
    integer::i,j                     !// LABELS THE INTERNAL LOOP LEVELS

    integer::index_ab                !// Nv CHOSE 2
    integer::i_dum     

    real(real8)::Gij_internal        !// INTERNAL PART OF Gij
    real(real8),parameter::zero=real(0.0, real8)        
    
    !//////////// EXTRA VARIABLES
    
#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num()+1
#else
    integer,parameter::threadID = 1
#endif

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// STORE THE WEIGHT OF THE LAMBDA_PATH
    
    lambda_path%weight = sum(lambda_path%arc_weights(0:num_internal))
    start_elec = sum(lambda_path%occupations(0:num_internal))

    if (skip_this_internal(lambda_path%weight,start_elec)) return 

    index_ab = num_external*(num_external-1)/2

    jbar = 0 
    ibar = 0
    Gij_internal = zero 

    loop_type = lambda_path%loop_type
    internal_singles = lambda_path%num_singles
    j = lambda_path%constraints(2,1)
    i = lambda_path%constraints(2,2)


!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// INITIALIZE VARIABLES FOR BUILDING MU PATH
    top_level = lambda_path%level1
    constraint_count = 1
    lambda_path%rt_loop_weight = 0
    path_elecs = sum(lambda_path%occupations(0:top_level-1))


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


    !// FOR THE RT_LOOP_WEIGHT ADD IN THE WEIGHTS OF THE RELEVANT PARTS
    !// OF THE HEAD AND TAIL PATHS.  THIS WILL END UP BEING THE MU PATH WEIGHT
    lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
         sum(lambda_path%arc_weights(0:lambda_path%level1-1))

    if (skip_this_internal(lambda_path%rt_loop_weight,path_elecs)) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !//
    !// GET SINGLES AND CONVERT TO IBAR, JBAR
    !// BUT DON'T WASTE OUR TIME IF WE HAVE A 
    !// {26} LOOP
    !//
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (loop_type /= 5) then
       call get_numbers_of_singles(s1,s2,s3,s4,lambda_path)

       if (loop_type == 1.or.loop_type == 4) then
          ibar = s1 + s2
          jbar = s1 + 1
       else
          ibar = s1 + s2
          jbar = s1 
       endif

    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !//
    !// COMPUTE Gij INTERNAL, AGAIN, DON'T
    !// WASTE OUR TIME IF THIS IS A {26} LOOP
    !//
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (loop_type /= 5) then

       Gij_internal = zero
       call make_gij_internal(lambda_path, i, j, Gij_internal,0)
       Gij_internal = Gij_internal + pp(index2(i,j))

       !// ADD IN (IJ|JJ) FOR {53} AND {75}
       if (loop_type >= 3) then
          Gij_internal = Gij_internal + integral_buffer(1,j)
       endif

       !// FINALLY, ADD (IJ|II) FOR {75} AND {35}
       if (loop_type == 4.or.loop_type == 2) then
          Gij_internal = Gij_internal + integral_buffer(1,i)
       endif

    endif


    twoVars%twostwo_rec = twoVars%twostwo_rec + 1
    i_dum = i*(i-1)/2+j  
    
    call twostwo_seg_writeData(twoVars%twostwo_rec,pseudo_twostwo_info,threadID,i_dum,lambda_path%weight,lambda_path%rt_loop_weight,ibar,jbar,internal_singles,&
         Gij_internal,loop_type,start_elec,i,j)

    if (internal_singles .gt. 0) then
       call for_int_buf_writearray(twosingles_no, twoVars%twosingles_rec, internal_singles, lambda_path%singles(1:internal_singles), threadID)
       twoVars%twosingles_rec = twoVars%twosingles_rec + internal_singles
    endif

  end subroutine intern_two_seg_compl_vec_lmo
  !*****************************************************************
  subroutine pure_int_2_seg_comp_vec_lmo_res(civec, sigmavec, sigmavec_lock, ij_1,twostwo_rec,twosingles_rec,virt_lam_allow,&
                                             virt_mu_allow, ipjp_rec)

    !// IN THIS ROUTINE WE BUILD THE COMPLEMENT TO THE TWO
    !// SEGMENT LOOPS WHICH RESIDE ENTIRELY IN THE INTERNAL SPACE
    !// THIS IS A LONG COMPLICATED AND ARDUOUS ROUTINE.  NOT FOR
    !// THE FAINT OF HEART.  
    !// THE CURRENT (HOPEFULLY SOON OLD) VERSION OF THIS WAS A
    !// NIGHTMARE WAITING TO HAPPEN. THE COMPLETE NEGLECT OF LOOP
    !// FUSION IS ASTONISHING! (JMD)

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!

    use spin_mod 
    use get_integrals_mod,ONLY:prpr
    
    implicit none
   
    real(kind=real8), dimension(:) :: civec, sigmavec
    type(blockedLockType) :: sigmavec_lock
    integer,intent(inout)::twostwo_rec,twosingles_rec
    integer,intent(in)::ipjp_rec
    integer,dimension(:),intent(inout)::virt_lam_allow,virt_mu_allow
    integer,intent(inout)::ij_1
   
    integer::ibar,jbar               !// ORBITAL POSITION INDICES
    integer::internal_singles        !// NUMBER OF SINGLES IN INTERNAL PATH OF PATH
    integer::allocatestatus          !// FOR DYNAMIC MEMORY
    integer::deallocatestatus
    integer::lambda_singles          !// SINGLES IN LAMBDA AND MU PATH
    integer::mu_singles
    integer::lambda_dim,mu_dim       !// DIMENSIONS FOR LAMBDA AND MU CONFIGURATIONS
    integer::lambda_weight           !// WEIGHTS OF LAMBDA AND MU PATH
    integer::mu_weight              

    integer::start_elec              !// OCCUPATION WHERE INTERNAL CSF MEETS UP WITH EXTERNAL SPACE
    integer::loop_type               !// LABELS THE LOOP
    integer::i,j                     !// LABELS THE INTERNAL LOOP LEVELS

    integer::v1_address,v2_address   !// ADDRESSES FOR VALENCE-VALENCE INTERACTIONS
    integer::a1_address,a2_address   !// ADDRESSES FOR N-1/N-1 INTERACTIONS
    integer::v1_spin,v2_spin         !// INDICES FOR VALENCE STATE SPINS
    integer::a1_spin,a2_spin         !// INDICES FOR N-1 STATE SPINS
    integer::v1_start,v2_start       !// STARTING INDEX FOR CI AND SIGMA VECTOR FOR VALENCE STATES
    integer::a1_start,a2_start       !// STARTING INDEX FOR CI AND SIGMA VECTOR FOR N-1  STATES 
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
    integer::mode_1,mode_2           !// MODES FOR SINGLET AND TRIPLET COUPLINGS

    integer::lam_length1,mu_length1
    integer::lam_length2,mu_length2,lam_length2C2,mu_length2C2
    integer::i_dum,j_dum,ij_2

    real(real8)::Gij_internal        !// INTERNAL PART OF Gij
    real(real8),parameter::zero=real(0.0, real8),two=real(2.0, real8)      !// THE NUMBERS ZERO AND TWO
    real(real8),parameter::sqrt2=real(sqrt(real(2.0, real8)), real8)          !// SQRT2 
    
    real(real8)::the_integral        !// THE INTEGRAL FOR LOOP TYPE 5 
    integer::matsize
    real(real8),dimension(:,:),allocatable::Px                    !// COUPLING COEFFFICIENT MATRIX 
    real(real8),dimension(:,:),allocatable::Kab                   !// PARTIAL MATRIX ELEMENTS
    real(real8),dimension(:,:),allocatable::tmpmat                   !// tmpmat
    !real(real8),dimension(:),allocatable::Gij_vector              !// Gij_internal + integral  

    !//////////// EXTRA VARIABLES

    integer::common_virt_length,lam_offset,mu_offset,i_x,j_x,i_y,j_y,lam_offset2,mu_offset2
    !>todo we MUST also check if we want to size these allocations correctly...
    integer,dimension(:),allocatable::virt_lam_long,virt_lam_pos,virt_mu_long,virt_mu_pos,&
         common_virt_lam_pos,common_virt_mu_pos
    integer,dimension(:),allocatable::orb_sing
    real(real8),dimension(:),allocatable::sig_nm2_lam,sig_nm2_mu
    real(real8),dimension(:),allocatable::sig_nm2_lam_d,sig_nm2_mu_d,sig_nm1_lam,sig_nm1_mu
    
    
    real(real8),dimension(:),allocatable::tmp_arr1,tmp_arr2
    
#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num()+1
#else
    integer,parameter::threadID = 1
#endif

    allocate(sig_nm2_lam_d(num_external),sig_nm2_mu_d(num_external),sig_nm1_lam(num_external),sig_nm1_mu(num_external), &
    sig_nm2_lam(num_external*(num_external-1)/2),sig_nm2_mu(num_external*(num_external-1)/2),virt_lam_long(num_orbitals), &
         virt_lam_pos(num_orbitals),virt_mu_long(num_orbitals),virt_mu_pos(num_orbitals),orb_sing(num_orbitals),&
         common_virt_lam_pos(num_orbitals),common_virt_mu_pos(num_orbitals),stat=allocatestatus)
    call allocatecheck(allocatestatus,"the first bunch")

    ! START AUTOGENERATED INITIALIZATION 
    i_x = 0
    i_y = 0
    j_y = 0
    ab_mu_end = 0
    mu_dim = 0
    aa_mu_start = 0
    aa1_spin = 0
    v2_start = 0
    v2_spin = 0
    mode_2 = 0
    aa2_spin = 0
    v1_start = 0
    ab_lam_end = 0
    ab_mu_start = 0
    aa_mu_end = 0
    j_x = 0
    a1_address = 0
    a1_start = 0
    ab_lam_add = 0
    mu_length1 = 0
    aa_lam_start = 0
    a2_address = 0
    v2_address = 0
    a1_end = 0
    lam_offset = 0
    a2_end = 0
    aa_lam_end = 0
    mu_offset = 0
    aa_mu_add = 0
    a2_start = 0
    ab_mu_add = 0
    aa_lam_add = 0
    a2_spin = 0
    a1_spin = 0
    ! END AUTOGENERATED INITIALIZATION

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// ALLOCATE MATRIX ELEMENT STORAGE ARRAY
    matsize=fsn(open_shells)
    allocate(Kab(matsize,matsize),Px(matsize,matsize),tmpmat(matsize,matsize),stat=allocatestatus)
    call allocatecheck(allocatestatus,"PxKab")

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    i_dum = -1
    twostwo_rec = twostwo_rec + 1
    ij_2 = -1000               ! in theory the next i/o statement could fail resulting 
    ! in ij_2 being uninitialized
    ! here we'll initialize it and to a weird value to help debug if this is
    ! ever an issue
    call twostwo_seg_readData(twostwo_rec,pseudo_twostwo_info,threadID,ij_2,lambda_weight,mu_weight,ibar,jbar,internal_singles,&
         Gij_internal,loop_type,start_elec,i,j,i_dum)

    if (ij_1 /= ij_2) then
       twostwo_rec = twostwo_rec -1
       ij_1 = 0
       return
    endif

    if (i_dum /= 0) then ! jmd: this check was previously in the code doing absolutely NIL, I assume it was supposed to do the same
                         ! as the similar check in abcd_seg_read, checking a return/error condition flag. so now we are using this
                         ! for that.
       ij_1 = 0
       twostwo_rec = twostwo_rec -1
       return
    endif

    orb_sing = 0
    if (internal_singles .gt. 0) then
       call for_int_buf_readarray(twosingles_no, twosingles_rec, internal_singles, orb_sing(1:internal_singles), threadID)
       twosingles_rec = twosingles_rec + internal_singles
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// 
    !// NOW LOOP OVER SETS OF VIRTUALS,
    !// COMPUTE MATRIX ELEMENTS, AND THEN
    !// DO THE MULTIPLICATION
    !//
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// THREE CASES HERE: ONE FOR THE D,S,AND V VERTICES
    if (start_elec == num_elec-2) then

!!!!!!!!!!!!!!!!!
       !// 
       !// TWO VIRTUALS SINGLY OCCUPIED  AND ONE VIRTUAL DOUBLY OCCUPIED            
       !// pure_int_2_seg_comp_vec_lmo_res
!!!!!!!!!!!!!!!!!

       !// GET DIMENSIONS
       if (loop_type == 2 .or. loop_type == 3) then
          lambda_singles = internal_singles +2  
          mu_singles = lambda_singles - 2
       else
          lambda_singles = internal_singles + 2
          mu_singles = lambda_singles 
       endif

       lambda_dim = fsn(lambda_singles)
       mu_dim = fsn(mu_singles)        

       !// NOW WE CAN COMPUTE PX PART OF COUPLING COEFFICIENTS    
       if (loop_type == 2 .or. loop_type == 3) then 
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)
       else
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)
       endif


       ab_lam_add = internal_index_vector2(lambda_weight +1)
       ab_mu_add  = internal_index_vector2(mu_weight +1)

       aa_lam_add = internal_index_vector3(lambda_weight +1)
       aa_mu_add  = internal_index_vector3(mu_weight +1)



       if (lambda_dim /= 0.and.mu_dim /= 0) then

          if (ab_lam_add > 0 .and. ab_mu_add > 0) then  

             lam_length2  = num_allowed_virtuals(lambda_weight+1,"D") 
             lam_length2C2 = lam_length2*(lam_length2-1)/2
             mu_length2   = num_allowed_virtuals(mu_weight+1,"D")
             mu_length2C2   = mu_length2*(mu_length2-1)/2 


             ! Don't use i,j ! 
             call get_virtuals(lambda_weight+1,"D",virt_lam_allow)
             call get_virtuals(mu_weight+1,"D",virt_mu_allow)

             ! Making Changes for the lmo version

             !common_virt = 0
             virt_lam_long = 0
             virt_lam_pos = 0
             virt_mu_long = 0
             virt_mu_pos = 0

             ! Keep the virt_allowed in long format and keep track of position
             do i_dum = 1, lam_length2
                j_dum = virt_lam_allow(i_dum)
                virt_lam_long(j_dum) = j_dum
                virt_lam_pos(j_dum) = i_dum
             enddo

             do i_dum = 1, mu_length2
                j_dum = virt_mu_allow(i_dum)
                virt_mu_long(j_dum) = j_dum
                virt_mu_pos(j_dum) = i_dum
             enddo

             ! Get the common virtuals between lam,mu and also the reverse tracking position

             common_virt_length = 0
             common_virt_lam_pos = 0
             common_virt_mu_pos = 0
             do i_dum = 1, lam_length2
                j_dum = virt_lam_allow(i_dum)
                if ( virt_lam_long(j_dum) == virt_mu_long(j_dum) ) then
                   common_virt_length = common_virt_length + 1
                   !common_virt(common_virt_length) = j_dum
                   common_virt_lam_pos(common_virt_length) = virt_lam_pos(j_dum)
                   common_virt_mu_pos(common_virt_length) = virt_mu_pos(j_dum)
                endif
             enddo


             !// NOW HANDLE THE {26} MULTIPLICATION IF THIS IS
             !// A {26} LOOP. HANDLE OTHER LOOPS LATER 
             if (loop_type == 5) then


                the_integral = prpr(index2m1(i,j))

                do aa1_spin  = 1,lambda_dim 

                   mode_1 = -1
                   if (aa1_spin <= fsn(lambda_singles-2)) mode_1 = 1

                   ab_lam_start = ab_lam_add   + lam_length2C2*(aa1_spin - 1)
                   ab_lam_end   = ab_lam_start + lam_length2C2- 1

                   ab_mu_start = ab_mu_add   + mu_length2C2*(aa1_spin - 1)
                   ab_mu_end   = ab_mu_start + mu_length2C2- 1
                   !
                   aa_lam_start = aa_lam_add   + lam_length2*(aa1_spin - 1)
                   aa_lam_end   = aa_lam_start + lam_length2- 1

                   aa_mu_start  = aa_mu_add   + mu_length2*(aa1_spin - 1)
                   aa_mu_end    = aa_mu_start + mu_length2- 1

                   !///////////////////////////////// For the LMO version //////////////////////////////////////////////////////////////////


                   sig_nm2_lam = 0.0D0
                   sig_nm2_lam_d = 0.0D0
                   sig_nm2_mu = 0.0D0
                   sig_nm2_mu_d = 0.0D0

                   if (mode_1 == 1) then
                      do i_dum = 1, common_virt_length
                         i_x = common_virt_lam_pos(i_dum)
                         i_y = common_virt_mu_pos(i_dum)
                         lam_offset = (i_x-1)*(i_x-2)/2 
                         mu_offset = (i_y-1)*(i_y-2)/2 + ab_mu_start -1 ! Include the (ab_mu_start-1) here
                         lam_offset2 = (i_x-1)*(i_x-2)/2 + ab_lam_start-1 ! Include the (ab_lam_start-1) here
                         mu_offset2 = (i_y-1)*(i_y-2)/2 

                         do j_dum = 1,i_dum-1
                            j_x = common_virt_lam_pos(j_dum)
                            j_y = common_virt_mu_pos(j_dum)
                            sig_nm2_lam(lam_offset+j_x) = the_integral*civec(mu_offset+j_y)
                            sig_nm2_mu(mu_offset2+j_y) = the_integral*civec(lam_offset2+j_x)
                         enddo
                         sig_nm2_lam_d(i_x) = the_integral*civec(aa_mu_start+i_y-1)
                         sig_nm2_mu_d(i_y) = the_integral*civec(aa_lam_start+i_x-1)
                      enddo

                   else

                      do i_dum = 1, common_virt_length
                         i_x = common_virt_lam_pos(i_dum)
                         i_y = common_virt_mu_pos(i_dum)
                         lam_offset = (i_x-1)*(i_x-2)/2
                         mu_offset = (i_y-1)*(i_y-2)/2 + ab_mu_start -1 ! Include the (ab_mu_start-1) here
                         lam_offset2 = (i_x-1)*(i_x-2)/2 + ab_lam_start-1 ! Include the (ab_lam_start-1) here
                         mu_offset2 = (i_y-1)*(i_y-2)/2

                         do j_dum = 1,i_dum-1
                            j_x = common_virt_lam_pos(j_dum)
                            j_y = common_virt_mu_pos(j_dum)
                            sig_nm2_lam(lam_offset+j_x) = the_integral*civec(mu_offset+j_y)
                            sig_nm2_mu(mu_offset2+j_y) = the_integral*civec(lam_offset2+j_x)
                         enddo
                      enddo

                   endif

                   !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef TIGER_USE_OMP
                   call setLocks(sigmavec_lock,ab_lam_start,ab_lam_end)
#endif
                   sigmavec(ab_lam_start:ab_lam_end) = sigmavec(ab_lam_start:ab_lam_end) + sig_nm2_lam(1:lam_length2C2)
#ifdef TIGER_USE_OMP
                   call unsetLocks(sigmavec_lock,ab_lam_start,ab_lam_end)
#endif

                   if (aa_lam_add > 0) then
#ifdef TIGER_USE_OMP
                       call setLocks(sigmavec_lock,aa_lam_start,aa_lam_end)
#endif
                       sigmavec(aa_lam_start:aa_lam_end) = sigmavec(aa_lam_start:aa_lam_end) + sig_nm2_lam_d(1:lam_length2)
#ifdef TIGER_USE_OMP
                       call unsetLocks(sigmavec_lock,aa_lam_start,aa_lam_end)
#endif
                   endif


#ifdef TIGER_USE_OMP
                   call setLocks(sigmavec_lock,ab_mu_start,ab_mu_end)
#endif
                   sigmavec(ab_mu_start:ab_mu_end) = sigmavec(ab_mu_start:ab_mu_end) + sig_nm2_mu(1:mu_length2C2)
#ifdef TIGER_USE_OMP
                   call unsetLocks(sigmavec_lock,ab_mu_start,ab_mu_end)
#endif

                   if (aa_mu_add > 0) then
#ifdef TIGER_USE_OMP
                       call setLocks(sigmavec_lock,aa_mu_start,aa_mu_end)
#endif
                       sigmavec(aa_mu_start:aa_mu_end) = sigmavec(aa_mu_start:aa_mu_end) + sig_nm2_mu_d(1:mu_length2)
#ifdef TIGER_USE_OMP
                       call unsetLocks(sigmavec_lock,aa_mu_start,aa_mu_end)
#endif
                   endif
                enddo


             elseif (loop_type /= 5) then

                !// FORM INTERNAL PART OF Kab.  NOTICE HOW WE ONLY
                !// DO THIS IF WE KNOW WE ARE GOING TO NEED IT
                call compute_Kab_2(internal_singles,orb_sing(1:internal_singles),i, j,loop_type, &
                     ibar, jbar, &
                     lambda_singles, mu_singles, &
                     Kab,size(Kab,1),0,ipjp_rec)

                !//THE FOUR INTEGRALS WHICH ARE INVOLVED IN THIS INTERACTION ARE
                !//INCLUDED IN IJAB ROUTINE. 
                !//((aa|ij) + (bb|ij))A(Px) + (ai|aj)A(Py(s)) + (bi|bj)A(Py(s-1))
                !//(2(aa|ij) - (ai|aj))A(Px) 


                !// NOW, WE ADD IN Gij TO THE DIAGONAL ELEMENTS OF Kab ACCORDING TO LOOP_TYPE
                if (loop_type == 4) then
                   do aa1_spin  = 1,lambda_dim 
                      Kab(aa1_spin,aa1_spin)= Kab(aa1_spin,aa1_spin) - Gij_internal 
                   enddo
                else
                   do aa1_spin  = 1,lambda_dim 
                      Kab(aa1_spin,aa1_spin)= Kab(aa1_spin,aa1_spin) + Gij_internal 
                   enddo
                endif

                !// DO THE MATMUL OF COUPLING COEFFICIENTS WITH Ka TO OBTAIN MATRIX ELEMENTS
                call dgemm('N','N',lambda_dim,mu_dim,lambda_dim,1.d0,Kab,matsize,Px,matsize,0.d0,tmpmat,matsize)
                do aa1_spin  = 1,lambda_dim 
                
                   ab_lam_start = ab_lam_add   + lam_length2C2*(aa1_spin - 1)
                   ab_lam_end   = ab_lam_start + lam_length2C2- 1

                   aa_lam_start = aa_lam_add   + lam_length2*(aa1_spin - 1)
                   aa_lam_end   = aa_lam_start + lam_length2- 1

                   mode_1 = -1
                   if (aa1_spin <= fsn(lambda_singles-2)) mode_1 = 1

                   sig_nm2_lam = 0.0D0
                   sig_nm2_lam_d = 0.0D0
                   do aa2_spin  = 1,mu_dim 

                      !
                      ab_mu_start = ab_mu_add   + mu_length2C2*(aa2_spin - 1)
                      ab_mu_end   = ab_mu_start + mu_length2C2- 1

                      aa_mu_start  = aa_mu_add   + mu_length2*(aa2_spin - 1)
                      aa_mu_end    = aa_mu_start + mu_length2- 1

                      mode_2 = -1
                      if (aa2_spin <= fsn(mu_singles-2)) mode_2 = 1

                      !///////////////////////////////// For the LMO version //////////////////////////////////////////////////////////////////

                      if ( mode_2 == 1) then    ! DEBUG: Changed mode_1 into mode_2 because aa_mu_start is the critical vairable here. CMK

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2 
                            mu_offset = (i_y-1)*(i_y-2)/2 + ab_mu_start-1 ! Include (ab_mu_start-1) here 

                            do j_dum = 1,i_dum-1
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_lam(lam_offset+j_x) = sig_nm2_lam(lam_offset+j_x) + tmpmat(aa1_spin,aa2_spin)*civec(mu_offset+j_y)
                            enddo
                            sig_nm2_lam_d(i_x) = sig_nm2_lam_d(i_x) + tmpmat(aa1_spin,aa2_spin)*civec(aa_mu_start+i_y-1)
                         enddo

                      else

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2
                            mu_offset = (i_y-1)*(i_y-2)/2 + ab_mu_start-1 ! Include (ab_mu_start-1) here

                            do j_dum = 1,i_dum-1
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_lam(lam_offset+j_x) = sig_nm2_lam(lam_offset+j_x) + tmpmat(aa1_spin,aa2_spin)*civec(mu_offset+j_y)
                            enddo
                         enddo

                      endif



                      !///////////////////////////////// For the LMO version //////////////////////////////////////////////////////////////////
                      !TODO AND GUESS WHAT? LOOP FUSION!


                      sig_nm2_mu = 0.0D0
                      sig_nm2_mu_d = 0.0D0

                      if (mode_1 == 1) then ! DEBUG: Changed mode_2 into mode_1 because aa_lam_start is the critical vairable here. CMK

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2 + ab_lam_start-1 ! Include (ab_lam_start-1) here
                            mu_offset = (i_y-1)*(i_y-2)/2 

                            do j_dum = 1, i_dum-1 
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_mu(mu_offset+j_y) = tmpmat(aa1_spin,aa2_spin)*civec(lam_offset+j_x)
                            enddo
                            sig_nm2_mu_d(i_y) = tmpmat(aa1_spin,aa2_spin)*civec(aa_lam_start+i_x-1)
                         enddo

                      else

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2 + ab_lam_start-1 ! Include (ab_lam_start-1) here
                            mu_offset = (i_y-1)*(i_y-2)/2

                            do j_dum = 1, i_dum-1
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_mu(mu_offset+j_y) = tmpmat(aa1_spin,aa2_spin)*civec(lam_offset+j_x)
                            enddo
                         enddo

                      endif

                      !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef TIGER_USE_OMP
                      call setLocks(sigmavec_lock,ab_mu_start,ab_mu_end)
#endif
                      sigmavec(ab_mu_start:ab_mu_end) = sigmavec(ab_mu_start:ab_mu_end) + sig_nm2_mu(1:mu_length2C2)
#ifdef TIGER_USE_OMP
                      call unsetLocks(sigmavec_lock,ab_mu_start,ab_mu_end)
#endif

                      if (aa_mu_add > 0) then
#ifdef TIGER_USE_OMP
                              call setLocks(sigmavec_lock,aa_mu_start,aa_mu_end)
#endif
                              sigmavec(aa_mu_start:aa_mu_end) = sigmavec(aa_mu_start:aa_mu_end) + sig_nm2_mu_d(1:mu_length2)
#ifdef TIGER_USE_OMP
                              call unsetLocks(sigmavec_lock,aa_mu_start,aa_mu_end)
#endif
                      endif


                   enddo

#ifdef TIGER_USE_OMP
                   call setLocks(sigmavec_lock,ab_lam_start,ab_lam_end)
#endif
                   sigmavec(ab_lam_start:ab_lam_end) = sigmavec(ab_lam_start:ab_lam_end) + sig_nm2_lam(1:lam_length2C2)
#ifdef TIGER_USE_OMP
                   call unsetLocks(sigmavec_lock,ab_lam_start,ab_lam_end)
#endif

                   if (aa_lam_add > 0) then
#ifdef TIGER_USE_OMP
                           call setLocks(sigmavec_lock,aa_lam_start,aa_lam_end)
#endif
                           sigmavec(aa_lam_start:aa_lam_end) = sigmavec(aa_lam_start:aa_lam_end) + sig_nm2_lam_d(1:lam_length2)
#ifdef TIGER_USE_OMP
                           call unsetLocks(sigmavec_lock,aa_lam_start,aa_lam_end)
#endif
                   endif

                enddo

             endif

          endif

       endif

    elseif (start_elec == num_elec -1) then

!!!!!!!!!!!!!!!!!!!!
       !//                 
       !// ONE VIRTUAL SINGLY OCCUPIED             
       !//                                         
!!!!!!!!!!!!!!!!!!!!

       !// GET DIMENSIONS
       if (loop_type == 2 .or. loop_type == 3) then
          lambda_singles = internal_singles + 1
          mu_singles = lambda_singles - 2
       else
          lambda_singles = internal_singles + 1
          mu_singles = lambda_singles 
       endif

       lambda_dim = fsn(lambda_singles)
       mu_dim = fsn(mu_singles)

       !// GET PX PART OF COUPLING COEFFICIENTS
       if (loop_type == 2 .or. loop_type == 3) then 
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)
       else
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)
       endif

       a1_address = internal_index_vector1(lambda_weight +1)
       a2_address = internal_index_vector1(mu_weight +1)

       !// MAKE SURE DIMENSIONS ARE NON ZERO    
       if (lambda_dim /= 0.and.mu_dim /= 0) then


          if (a1_address > 0 .and. a2_address > 0 ) then 

             lam_length1   = num_allowed_virtuals(lambda_weight+1,"S") 
             mu_length1    = num_allowed_virtuals(mu_weight+1,"S") 


             call get_virtuals(lambda_weight+1,"S",virt_lam_allow)
             call get_virtuals(mu_weight+1,"S",virt_mu_allow)

             ! Making Changes for the lmo version

             !common_virt = 0
             virt_lam_long = 0
             virt_lam_pos = 0
             virt_mu_long = 0
             virt_mu_pos = 0

             ! Keep the virt_allowed in long format and keep track of position
             do i_dum = 1, lam_length1
                j_dum = virt_lam_allow(i_dum)
                virt_lam_long(j_dum) = j_dum
                virt_lam_pos(j_dum) = i_dum
             enddo

             do i_dum = 1, mu_length1
                j_dum = virt_mu_allow(i_dum)
                virt_mu_long(j_dum) = j_dum
                virt_mu_pos(j_dum) = i_dum
             enddo

             ! Get the common virtuals between lam,mu and also the reverse tracking position

             common_virt_length = 0
             common_virt_lam_pos = 0
             common_virt_mu_pos = 0
             do i_dum = 1, lam_length1
                j_dum = virt_lam_allow(i_dum)
                if ( virt_lam_long(j_dum) == virt_mu_long(j_dum) ) then
                   common_virt_length = common_virt_length + 1
                   !common_virt(common_virt_length) = j_dum
                   common_virt_lam_pos(common_virt_length) = virt_lam_pos(j_dum)
                   common_virt_mu_pos(common_virt_length) = virt_mu_pos(j_dum)
                endif
             enddo

             !          
             !// NOW HANDLE THE {26} MULTIPLICATION IF THIS IS
             !// A {26} LOOP. HANDLE OTHER LOOPS LATER 

             if (loop_type == 5) then

                the_integral = prpr(index2m1(i,j))

                do a1_spin  = 1,lambda_dim 

                   a1_start = a1_address + lam_length1*(a1_spin - 1)
                   a1_end   = a1_start   + lam_length1- 1

                   a2_start = a2_address + mu_length1*(a1_spin - 1)
                   a2_end   = a2_start   + mu_length1- 1

                   !///////////////////////////////// For the LMO version //////////////////////////////////////////////////////////////////

                   sig_nm1_lam= 0.0D0
                   sig_nm1_mu = 0.0D0
                   do i_dum = 1, common_virt_length
                      i_x = common_virt_lam_pos(i_dum)
                      i_y = common_virt_mu_pos(i_dum)
                      sig_nm1_lam(i_x) = the_integral*civec(a2_start+i_y-1)
                      sig_nm1_mu(i_y)  = the_integral*civec(a1_start+i_x-1)
                   enddo

                   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef TIGER_USE_OMP
                   call setLocks(sigmavec_lock,a1_start,a1_end)
#endif
                   sigmavec(a1_start:a1_end) = sigmavec(a1_start:a1_end) + sig_nm1_lam(1:lam_length1)
#ifdef TIGER_USE_OMP
                   call unsetLocks(sigmavec_lock,a1_start,a1_end)
                   call setLocks(sigmavec_lock,a2_start,a2_end)
#endif
                   sigmavec(a2_start:a2_end) = sigmavec(a2_start:a2_end) + sig_nm1_mu(1:mu_length1)
#ifdef TIGER_USE_OMP
                   call unsetLocks(sigmavec_lock,a2_start,a2_end)
#endif
                enddo

             elseif (loop_type /= 5) then

                !// FORM INTERNAL PART OF Kab.  NOTICE HOW WE ONLY
                !// DO THIS IF WE KNOW WE ARE GOING TO NEED IT
                call compute_Kab_2(internal_singles,orb_sing(1:internal_singles),i, j,loop_type, &
                     ibar, jbar, &
                     lambda_singles, mu_singles, &
                     Kab, size(Kab,1),0,ipjp_rec)


                !//THE TWO INTEGRALS WHICH ARE INVOLVED IN THIS INTERACTION ARE
                !//INCLUDED IN IJAB ROUTINE. 
                !//(aa|ij)A(Px) + (ai|aj)A(Py(s))

                !// NOW, WE ADD IN Gij TO THE DIAGONAL ELEMENTS OF Kab ACCORDING TO LOOP_TYPE
                if (loop_type == 4) then
                   do aa1_spin  = 1,lambda_dim 
                      Kab(aa1_spin,aa1_spin)= Kab(aa1_spin,aa1_spin) - Gij_internal 
                   enddo
                else
                   do aa1_spin  = 1,lambda_dim 
                      Kab(aa1_spin,aa1_spin)= Kab(aa1_spin,aa1_spin) + Gij_internal 
                   enddo
                endif

                !// DO THE MATMUL OF COUPLING COEFFICIENTS WITH Ka TO OBTAIN MATRIX ELEMENTS
                call dgemm('N','N',lambda_dim,mu_dim,lambda_dim,1.d0,Kab,matsize,Px,matsize,0.d0,tmpmat,matsize)

                do a1_spin  = 1,lambda_dim 

                   !
                   a1_start = a1_address + lam_length1*(a1_spin - 1)
                   a1_end   = a1_start   + lam_length1- 1


                   !                     if (index_vector_pao1(lambda_path%weight +1) < 0) cycle 
                   if (internal_index_vector1(lambda_weight +1) < 0) cycle

                   sig_nm1_lam= 0.0D0
                   do a2_spin  = 1,mu_dim

                      a2_start = a2_address + mu_length1*(a2_spin - 1)
                      a2_end   = a2_start   + mu_length1- 1

                      !                        if (index_vector_pao1(lambda_path%rt_loop_weight +1) < 0) cycle 
                      if (internal_index_vector1(mu_weight +1) < 0) cycle


                      !///////////////////////////////// For the LMO version //////////////////////////////////////////////////////////////////

                      sig_nm1_mu = 0.0d0
                      do i_dum = 1, common_virt_length
                         i_x = common_virt_lam_pos(i_dum)
                         i_y = common_virt_mu_pos(i_dum)
                         sig_nm1_lam(i_x) = sig_nm1_lam(i_x) + tmpmat(a1_spin,a2_spin)*civec(a2_start+i_y-1)
                         sig_nm1_mu(i_y) = tmpmat(a1_spin,a2_spin)*civec(a1_start+i_x-1)
                      enddo

                      !////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef TIGER_USE_OMP
                      call setLocks(sigmavec_lock,a2_start,a2_end)
#endif
                      sigmavec(a2_start:a2_end) = sigmavec(a2_start:a2_end) + sig_nm1_mu(1:mu_length1)
#ifdef TIGER_USE_OMP
                      call unsetLocks(sigmavec_lock,a2_start,a2_end)
#endif
                   enddo
                   
#ifdef TIGER_USE_OMP
                   call setLocks(sigmavec_lock,a1_start,a1_end)
#endif
                   sigmavec(a1_start:a1_end) = sigmavec(a1_start:a1_end) + sig_nm1_lam(1:lam_length1)
#ifdef TIGER_USE_OMP
                   call unsetlocks(sigmavec_lock,a1_start,a1_end)
#endif
                enddo

             endif

          endif

       endif


    elseif (start_elec == num_elec) then

!!!!!!!!!!!!!!!!!!
       !//               
       !// VALENCE STATE; NO EXTERNALS OCCUPIED        
       !//                                             
!!!!!!!!!!!!!!!!!!


       !// GET DIMENSIONS
       if (loop_type == 2 .or. loop_type == 3) then
          lambda_singles = internal_singles 
          mu_singles = lambda_singles - 2
       else
          lambda_singles = internal_singles 
          mu_singles = lambda_singles 
       endif

       lambda_dim = fsn(lambda_singles)
       mu_dim = fsn(mu_singles)

       !// GET PX PART OF COUPLING COEFFICIENT
       if (loop_type == 2 .or. loop_type ==3) then 
          call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)
       else
          call contracted_cycle_stateless(jbar,ibar,0,0,0,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)
       endif

       v1_address = internal_index_vector0(lambda_weight+1)
       v2_address = internal_index_vector0(mu_weight+1)


       !// MAKE SURE DIMENSIONS ARE NON ZERO    
       if (lambda_dim /= 0 .and. mu_dim /= 0) then

          if (v1_address > 0 .and. v2_address > 0 ) then 


             !// FORM INTERNAL PART OF Kab.  NOTICE HOW WE ONLY
             !// DO THIS IF WE KNOW WE ARE GOING TO NEED IT
             if (loop_type /= 5) then
                call compute_Kab_2(internal_singles,orb_sing(1:internal_singles),i, j,loop_type, &
                     ibar, jbar, &
                     lambda_singles, mu_singles, &
                     Kab, size(Kab,1),0,ipjp_rec)
             endif

             !// NOW HANDLE THE {26} MULTIPLICATION IF THIS IS
             !// A {26} LOOP.  THE OTHER LOOPS WILL BE HANDLED
             !// BELOW.
             if (loop_type == 5) then

                the_integral =  prpr(index2m1(i,j))
                
                allocate(tmp_arr1(lambda_dim),tmp_arr2(lambda_dim),stat=allocatestatus)
                call allocatecheck(allocatestatus,"tmp_arr12_5")
                tmp_arr1 = 0.0
                tmp_arr2 = 0.0

                v1_start = v1_address
                v2_start = v2_address
                do v1_spin  = 1,lambda_dim
                   tmp_arr1(v1_spin) = tmp_arr1(v1_spin) + the_integral*civec(v2_start)
                   tmp_arr2(v1_spin) = tmp_arr2(v1_spin) + the_integral*civec(v1_start)
                   v1_start = v1_start + 1
                   v2_start = v2_start + 1
                enddo
                
#ifdef TIGER_USE_OMP
                !call setSigmaLocks(v1_address,v1_address+lambda_dim-1)
                call setLocks(sigmavec_lock,v1_address,v1_address+lambda_dim-1)
#endif
                sigmavec(v1_address:v1_address+lambda_dim-1) = sigmavec(v1_address:v1_address+lambda_dim-1) + tmp_arr1
#ifdef TIGER_USE_OMP
                !call unsetSigmaLocks(v1_address,v1_address+lambda_dim-1)
                call unsetLocks(sigmavec_lock,v1_address,v1_address+lambda_dim-1)
                !call setSigmaLocks(v2_address,v2_address+lambda_dim-1)
                call setLocks(sigmavec_lock,v2_address,v2_address+lambda_dim-1)
#endif          
                sigmavec(v2_address:v2_address+lambda_dim-1) = sigmavec(v2_address:v2_address+lambda_dim-1) + tmp_arr2
#ifdef TIGER_USE_OMP
                !call unsetSigmaLocks(v2_address,v2_address+lambda_dim-1)
                call unsetLocks(sigmavec_lock,v2_address,v2_address+lambda_dim-1)
#endif
                
                deallocate(tmp_arr1,tmp_arr2,stat=deallocatestatus)
                call deallocatecheck(deallocatestatus,"tmp_arr12_5")
             else

                call add_in_Gij(Kab, lambda_dim, &
                     Gij_internal,&
                     loop_type)

                !// DO THE MATMUL OF COUPLING COEFFICIENTS WITH Ka TO OBTAIN MATRIX ELEMENTS
                call dgemm('N','N',lambda_dim,mu_dim,lambda_dim,1.d0,Kab,matsize,Px,matsize,0.d0,tmpmat,matsize)
                !Kab(1:lambda_dim,1:lambda_dim) = tmpmat(1:lambda_dim,1:lambda_dim)
                
                allocate(tmp_arr1(lambda_dim),tmp_arr2(mu_dim),stat=allocatestatus)
                call allocatecheck(allocatestatus,"tmp_arr12_5")
                tmp_arr1 = 0.0
                tmp_arr2 = 0.0

                v1_start = v1_address
                do v1_spin  = 1,lambda_dim
                   v2_start = v2_address
                   do v2_spin  = 1,mu_dim
                      tmp_arr1(v1_spin) = tmp_arr1(v1_spin) + tmpmat(v1_spin,v2_spin)*civec(v2_start)
                      tmp_arr2(v2_spin) = tmp_arr2(v2_spin) + tmpmat(v1_spin,v2_spin)*civec(v1_start)
                      v2_start = v2_start + 1
                   enddo
                   v1_start = v1_start + 1
                enddo
#ifdef TIGER_USE_OMP
                !call setSigmaLocks(v1_address,v1_address+lambda_dim-1)
                call setLocks(sigmavec_lock,v1_address,v1_address+lambda_dim-1)
#endif
                sigmavec(v1_address:v1_address+lambda_dim-1) = sigmavec(v1_address:v1_address+lambda_dim-1) + tmp_arr1
#ifdef TIGER_USE_OMP
                !call unsetSigmaLocks(v1_address,v1_address+lambda_dim-1)
                call unsetLocks(sigmavec_lock,v1_address,v1_address+lambda_dim-1)
                !call setSigmaLocks(v2_address,v2_address+mu_dim-1)
                call setLocks(sigmavec_lock,v2_address,v2_address+mu_dim-1)
#endif          
                sigmavec(v2_address:v2_address+mu_dim-1) = sigmavec(v2_address:v2_address+mu_dim-1) + tmp_arr2
#ifdef TIGER_USE_OMP
                !call unsetSigmaLocks(v2_address,v2_address+mu_dim-1)
                call unsetLocks(sigmavec_lock,v2_address,v2_address+mu_dim-1)
#endif
                
                deallocate(tmp_arr1,tmp_arr2,stat=deallocatestatus)
                call deallocatecheck(deallocatestatus,"tmp_arr12_e")

             endif

          endif

       endif

    endif


    !// CLEAN UP BEFORE EXITING

    deallocate(Px,Kab,tmpmat,stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"PxKab")

  end subroutine pure_int_2_seg_comp_vec_lmo_res
  !*****************************************************************
  !****************************************************************
  !*****************************************************************

  !*****************************************************************
  !*****************************************************************
  !****************************************************************
  Subroutine  extern_two_seg_compl_vec_lmo(civec, sigmavec, lambda_path,loc_scr,twoVars)
  
    use locist_var_mod,only:locist_scratch

    !// IN THIS SUBROUTINE WE BUILD THE TWO SEGMENT LOOPS WHICH RESIDE
    !// ENTIRELY IN THE EXTERNAL SPACE AND JOIN THEM UP TO THE PATHS
    !// WHICH GO FROM THE HEAD OF THE GRAPH TO THE INTERNAL-EXTERNAL
    !// SPACE BORDER.  THE TREATMENT HERE USES A FAST VECTORIZED MODE, BUT
    !// ALL OF THE OLD CODE IS HERE IN CASE YOU NEED TO DO SOME DEBUGGING
    !// OF THIS LATER ON.

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(orbital_path),intent(in)::lambda_path
    type(locist_scratch)::loc_scr
    type(twoModVars)::twoVars
    real(kind=real8), dimension(:) :: civec, sigmavec 

    
    real(real8),parameter::zero = real(0.0, real8)                       !// THE NUMBER ZERO
    real(real8),parameter::one = real(1.0, real8)                        !// THE NUMBER ONE
    real(real8),parameter::sqrt2  = real(sqrt(real(2.0,real8)),real8)                     !// SQUARE ROOT OF TWO
    real(real8),dimension(:,:),allocatable::Hab35                !// MATRIX ELEMENTS FOR {35} AND {53}

    integer::start_elec             !// NUMBER OF ELECTRONS IN THE INTERNAL SPACE
    integer::a,b                   !// LABELS THE LOOP LEVELS
    integer::internal_weight        !// WEIGHT OF INTERNAL PART OF PATH
    integer::internal_singles       !// NUMBER OF SINGLES IN INTERNAL PART OF PATH
    integer::internal_address       !// ADDRESS OF INTERNAL SPACE CSF IN INTERNAL INDEX
    integer::lambda_singles         !// NUMBER OF SINGLES IN PATHS
    integer::mu_singles
    integer::lambda_dim             !// DIMENSIONS OF MU PATH
    integer::mu_dim
    integer::lambda_singlet_dim     !// DIMENSIONS OF SPIN SPACE FOR SINGLETS
    integer::mu_singlet_dim
    integer::lambda_address         !// PATH ADDRESSES
    integer::mu_address
    integer::allocatestatus         !// FOR DYNAMIC MEMORY
    integer::deallocatestatus
    integer::num_virtual            !// NUMBER OF VIRTUAL ORBITALS
    integer::ab_count               !// FOR LOOPING OVER A,B
    integer::spin_lambda,spin_mu    !// FOR LOOPING OVER SPIN FUNCTIONS
    integer::mode_1,mode_2          !// FOR STORING THE SYMMETRIC/ANTISYMMETRIC SCEP STYLE

    integer::lam_length,lam_lengthC2       
    integer::mu_length,mu_lengthC2         
    integer::aa_address
    integer::ab_lam_start
    integer::ab_mu_start
    integer::aa_lam_start
    integer::aa_mu_start
    integer::tmp
    integer::status
    
    allocate(Hab35(num_external,num_external),stat=status)
    call allocatecheck(status,"extern_two_seg_compl_vec_lmo")
    
    ! START AUTOGENERATED INITIALIZATION 
    aa_address = 0
    mu_singlet_dim = 0
    ab_mu_start = 0
    aa_lam_start = 0
    mode_1 = 0
    mode_2 = 0
    mu_length = 0
    lambda_address = 0
    spin_mu = 0
    lambda_singles = 0
    mu_singles = 0
    mu_dim = 0
    a = 0
    b = 0
    lambda_singlet_dim = 0
    spin_lambda = 0
    lambda_dim = 0
    lam_lengthc2 = 0
    mu_address = 0
    mu_lengthc2 = 0
    ! END AUTOGENERATED INITIALIZATION 
    
!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    start_elec = sum(lambda_path%occupations(0:num_internal))
    internal_weight = sum(lambda_path%arc_weights(0:num_internal))

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (skip_this_internal(internal_weight, start_elec)) return

    internal_singles = lambda_path%num_singles
    num_virtual = num_orbitals-num_internal

    !////////////////////SEE WE HAVE NOW CHANGED THE DIMENSION OF GIJ_MATRIX
    !////////////////////ARUN WATCH OUT FOR DIMENSIONS OVERFLOW:
    allocate(twoVars%Gij_matrix(num_virtual*(num_virtual+1)/2), stat = allocatestatus)
    call allocatecheck(allocatestatus, "Gij_matr")
    twoVars%Gij_matrix = 0

    !// THERE ARE TWO CASES HERE.  FIRST, IF THE PATH ENDS AT NUM_ELEC-1 WE MAY
    !// BUILD ONLY THE {13} LOOP.  THE OTHER, MORE COMPLICATED, CASE IS WHEN THE
    !// INTERNAL PATH ENDS AT NUM_ELEC-2.  THIS MAKES THINGS A BITCH.  IN THIS CASE,
    !// WE CAN BUILD THE {75}, {35}, AND {53} LOOPS.  ADDITIONALLY, WE CAN BUILD THE {13}
    !// LOOP.

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

       lam_length = num_allowed_virtuals(internal_weight+1,"S")

       mu_length = lam_length

       !// MAKE THE GIJ MATRIX
       call make_gij_ab_vec_pao(lambda_path,loc_scr%virt_lam_allow,twoVars,internal_weight+1,lam_length,"S")

       internal_address = internal_index_vector1(internal_weight+1)
           if (internal_address > 0) then

       !// SET UP SINGLES AND DIMENSIONS
       lambda_singles     = internal_singles + 1
       mu_singles         = lambda_singles
       lambda_dim         = fsn(lambda_singles)
       mu_dim             = lambda_dim
       !lambda_singlet_dim = fsn(lambda_singles-2)
       !mu_singlet_dim     = lambda_singlet_dim
       !lines commented out due to possible out of bounds array ... besides these values aren't used in 
       ! in this if block anyway    David Krisiloff 6/2011


       !// COMPUTE MATRIX ELEMENTS AND DO MULTIPLICATIONS
       if (lambda_dim > 0) then
          call compute_Hab_vec_pao(lambda_path,loc_scr%virt_lam_allow,twoVars%Hab,twoVars%aibi_debug,twoVars%Gij_matrix, &
               lambda_singles,mu_singles,lambda_singles,internal_weight+1,lam_length,"S")

          !OMP TODO perhaps?
          do spin_lambda = 1, lambda_dim
             !                  if (internal_index_vector1(internal_weight+1) < 0) cycle

             lambda_address = internal_address + (spin_lambda-1)*lam_length

             do spin_mu = 1, spin_lambda

                !                      if (internal_index_vector1(internal_weight+1) < 0) cycle

                mu_address = internal_address + (spin_mu-1)*mu_length

                tmp = (max(spin_mu,spin_lambda)*(max(spin_mu,spin_lambda)-1)/2+min(spin_mu,spin_lambda))
                
                call dspmv('U',lam_length,one,twoVars%Hab(1,tmp),&
                     civec(mu_address),1,one,&
                     sigmavec(lambda_address),1)


                !//THE IF STATEMENT IS REQUIRED ..SEE ..ANALYTIC ..
                !// WE MAY HAVE TO USE THE TRANSPOSE OF THE INTEGRAL MATRIX HERE
                !// CHECK HERE .. PLEASE
                if (spin_mu /= spin_lambda) then

                   ! technically, this should be the transpose of Hab35, but since it is symmetric, this is not needed
                   call dspmv('U',lam_length,one,twoVars%Hab(1,tmp),&
                        civec(lambda_address),1,one,&
                        sigmavec(mu_address),1)


                endif

             enddo
          enddo
       endif

       deallocate(twoVars%Hab,stat=deallocatestatus)
       call deallocatecheck(deallocatestatus,"Hab     ")

       endif

!!!!!!!!!!!!!!!!!!!!!!!!
       !//
       !// NOW THE NUM_ELEC-2 CASE
       !//
!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (start_elec == num_elec - 2) then startelec
!!!!!!!!!!!!!!!!!!!!!!


       lam_length = num_allowed_virtuals(internal_weight+1,"D")

       mu_length = lam_length

       !// MAKE THEI GIJ MATRIX
       call make_gij_ab_vec_pao(lambda_path,loc_scr%virt_lam_allow,twoVars,internal_weight+1,lam_length,"D")

       internal_address = internal_index_vector2(internal_weight+1)
       aa_address       = internal_index_vector3(internal_weight+1)

       if (internal_address > 0) then 

       !// SET UP SINGLES AND DIMENSIONS
       lambda_singles     = internal_singles + 2
       mu_singles         = lambda_singles
       lambda_dim         = fsn(lambda_singles)
       mu_dim             = lambda_dim
       lambda_singlet_dim = fsn(lambda_singles-2)
       mu_singlet_dim     = lambda_singlet_dim


       !// COMPUTE MATRIX ELEMENTS AND DO MULTIPLICATIONS
       if (lambda_dim > 0) then

          call compute_Hab_vec_pao(lambda_path,loc_scr%virt_lam_allow,twoVars%Hab,twoVars%aibi_debug,twoVars%Gij_matrix, &
              lambda_singles,mu_singles,lambda_singles-1,internal_weight+1,lam_length,"D")


          lam_lengthC2 = lam_length*(lam_length-1)/2
          mu_lengthC2 = mu_length*(mu_length-1)/2

!!!!!!!!!!!!!!!!!!!
          !//
          !// {13}, {35}, {53} LOOPS
          !//
!!!!!!!!!!!!!!!!!!!

          !bs$bomp TODO perhaps?
          do spin_lambda = 1, lambda_dim

             ab_lam_start = internal_address + (spin_lambda-1)*lam_lengthC2
             aa_lam_start = aa_address + (spin_lambda-1)*lam_length

             mode_1 = -1
             if (spin_lambda <= lambda_singlet_dim) mode_1 = 1

             if (mode_1 .eq. 1) then
                 call scepper_diag_p1(twoVars%scep_ci_1,civec,ab_lam_start,&
                  lam_length,civec,aa_lam_start)
             elseif (mode_1 .eq. -1) then
                 call scepper_diag_m1(twoVars%scep_ci_1,civec,ab_lam_start,&
                  lam_length)
             endif

             loc_scr%sigma_nm2_lambda = 0.0d0
             do spin_mu = 1, spin_lambda

                ab_mu_start = internal_address + (spin_mu-1)*mu_lengthC2
                aa_mu_start = aa_address + (spin_mu-1)*mu_length

                mode_2 = -1
                if (spin_mu <= mu_singlet_dim) mode_2 = 1

                if (mode_2 .eq. 1) then
                    call scepper_diag_p1(twoVars%scep_ci_2,civec,ab_mu_start,&
                     mu_length,civec,aa_mu_start)
                elseif (mode_2 .eq. -1) then
                    call scepper_diag_m1(twoVars%scep_ci_2,civec,ab_mu_start,&
                     mu_length)
                endif
                     
                tmp = (max(spin_mu,spin_lambda)*(max(spin_mu,spin_lambda)-1)/2+min(spin_mu,spin_lambda))
                ab_count = 0
                do a=1,lam_length
                   do b=1,a
                      ab_count = ab_count + 1
                      Hab35(a,b) = twoVars%Hab(ab_count,tmp)
                      !Hab35(b,a) = Hab35(a,b) jmd: not needed
                   enddo
                enddo
                ! TODO remove Hab altogether. it is idiotic to begin with!


                !// DO THE MULTIPLICATIONS
                !// SEE EQUATIONS .. FOR DIMENSIONS..
                
                ! so there is no dspmm specified in BLAS (and they probably have good reasons for that)
                ! still it feels like cheating needing to scep -> dsymm -> unscepp...

                call dsymm('L','L',lam_length,mu_length,one,Hab35,num_external,&
                     twoVars%scep_ci_2,num_external,&
                     one,loc_scr%sigma_nm2_lambda,num_external)

                if (spin_mu /= spin_lambda) then
                   !// SEE EQUATIONS .. FOR DIMENSIONS..
                   !// technically, Hab35 should be transposed. Since Hab35 is symmetric, this is not needed
                   call dsymm('L','L',mu_length,lam_length,one,Hab35,num_external,&
                        twoVars%scep_ci_1,num_external,&
                        zero,loc_scr%sigma_nm2_mu,num_external)
                   if (mode_2 .eq. 1) then
                       call unscepper_diag_p1(loc_scr%sigma_nm2_mu,&
                        sigmavec,ab_mu_start,mu_length,sigmavec,aa_mu_start)
                   elseif (mode_2 .eq. -1) then
                       call unscepper_diag_m12(loc_scr%sigma_nm2_mu,&
                        sigmavec,ab_mu_start,mu_length)
                   endif

                endif
             enddo
             
             if (mode_1 .eq. 1) then
                 call unscepper_diag_p1(loc_scr%sigma_nm2_lambda,&
                   sigmavec,ab_lam_start,lam_length,sigmavec,aa_lam_start)
             elseif (mode_1 .eq. -1) then
                 call unscepper_diag_m12(loc_scr%sigma_nm2_lambda,&
                   sigmavec,ab_lam_start,lam_length)
             endif

          enddo

       endif

       deallocate(twoVars%Hab,stat=deallocatestatus)
       call deallocatecheck(deallocatestatus,"Hab     ")
       

       endif ! check for aa_address > 0

    endif startelec

    !// CLEAN UP
    deallocate(twoVars%Gij_matrix,Hab35,stat=status)
    call deallocatecheck(status,"extern_two_seg_compl_vec_lmo")

  end subroutine  extern_two_seg_compl_vec_lmo
  !*****************************************************************

  !******************************************************************
  subroutine compute_Hab_vec_pao(lambda_path,virt_lam_allow,Hab,aibi_debug,Gij_matrix,lambda_singles, mu_singles,s,weight,virt_length,state)

    !// FOR THE {13}, {35}, {53}, AND {57} LOOPS THE MATRIX ELEMENTS ALL
    !// INVOLVE A SUM OVER SINGLY OCCUPIED ORBITALS WHERE YOU
    !// MULTIPLY AN INTEGRAL BY A TRANSPOSITION.  HERE WE TREAT THIS SUM.
    !// THE FORMULAS ARE AS FOLLOWS:
    !// 
    !//    {13} = SUM { DELTA(Np,1) (ibar, pbar) (IP|JP) }
    !//    {35} = SUM { DELTA(Np,1) (jbar, pbar) (IP|JP) }
    !//    {53} = SUM { DELTA(Np,1) (ibar, pbar) (IP|JP) }
    !//    {57} = SUM { DELTA(Np,1) ((ibar, pbar) + 1) (IP|JP) }
    !// 
    !// THIS ROUTINE COMPUTES A KAB MATRIX FOR ALL A AND B TO USE
    !// IN THE VECTORIZED TREATMENT.


!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    use locist_var_mod,only:locist_scratch

    implicit none

    type(orbital_path),intent(in)::lambda_path
    integer,dimension(:),intent(inout)::virt_lam_allow
    real(real8),dimension(:,:),intent(inout),allocatable::Hab
    real(real8),dimension(:),intent(in)::aibi_debug
    real(real8),dimension(:),intent(in)::Gij_matrix

    integer,intent(in)::lambda_singles                 !// NUMBERS OF SINGLES
    integer,intent(in)::mu_singles
    integer,intent(in)::s
    character,intent(in)::state  
    integer,intent(in)::weight
    integer,intent(in)::virt_length                    !// LENGTH OF THE VIRTUALS
    
    integer:: j,k                          !// USED FOR LOOPING OVER LEVELS
    integer::a,b                            !// VIRTUALS
    integer::pbar                           !// KEEPS TRACK OF POSITION INDICES
    integer::lambda_dim                    !// DIMENSIONS OF MATRICES
    integer::transposition_label            !// LABELS THE TRANSPOSITION SUM
    integer::row_offset
    integer::orbital
    integer::offset                         !// FOR INTEGRAL ACCESS
    integer::allocatestatus
    integer::a1,b1 
    integer::idum,idum2
    real(real8)::transposition_element      !// ELEMENT OF A TRASPOSITION MATRIX
    integer::tmp

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lambda_dim = fsn(lambda_singles)
    !mu_dim = lambda_dim

    ! we assume, that the current code calling pattern holds true. Namely, lambda_singles
    ! and mu_singles are the same value
    if(lambda_singles /= mu_singles) then
       write(ioOutput,*) "ERROR: Assumption that mu_singles and lambda_singles are the same does NOT hold true."
       flush(ioOutput)
       stop
    endif

    allocate(Hab(virt_length*(virt_length+1)/2,((lambda_dim*lambda_dim-lambda_dim)/2+lambda_dim)),&
         stat=allocatestatus)
    call allocatecheck(allocatestatus,"Hab     ")
    Hab = real(0.0,real8)

    !// GET THE ALLOWED VIRTUALS 
    call get_virtuals(weight,state,virt_lam_allow)

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// LOOP OVER SINGLY OCCUPIED ORBITALS IN THE INTERNAL SPACE

    idum = num_internal*num_external*(num_external+1)/2
    !$omp parallel &
    !$omp default(none) &
    !$omp private(row_offset,tmp,transposition_element,a1,b1,idum2) &
    !$omp shared(lambda_dim,transposition_label,virt_length,Hab,aibi_debug,num_internal,offset, &
    !$omp transpositions,virt_lam_allow,ignorable_pair,orbital,Gij_matrix,lambda_path,idum,s,num_external)
    do pbar = 1, lambda_path%num_singles

       !// GET THE ORBITAL LABEL
       orbital = lambda_path%singles(pbar)

       offset = idum + (orbital-1)*num_external*(num_external+1)/2

       !// HERE WE NEED (ibar, pbar) TO TREAT EVERYTHING BUT
       !// THE {35} LOOPS.   IN THE CASE OF THE {57} LOOPS
       !// WE WILL ADD IN THE IDENTITY AFTER
       if (pbar /= s) then

          transposition_label = index2m1(pbar, s)
          
          !$omp do schedule(static)
          do j = 1, lambda_dim
             row_offset = j*(j-1)/2
             do k = 1, j-1
                tmp = j*(j-1)/2+k
                
                transposition_element = transpositions(transposition_label, row_offset + k)
                
                do a = 1,virt_length 
                   a1 = virt_lam_allow(a) 
                   if (ignorable_pair(a1,orbital) ) cycle
                   a1 = a1-num_internal
                   idum2 = a*(a-1)/2
                   do b = 1, a
                      idum2= idum2 + 1
                      b1 = virt_lam_allow(b) 
                      if (ignorable_pair(b1,orbital)) cycle
                      b1 = b1-num_internal
                      
                      Hab(idum2,tmp) = Hab(idum2,tmp) + aibi_debug(index2(a1,b1)+offset)*transposition_element
                   enddo
                enddo
             enddo

             transposition_element = transpositions(transposition_label, row_offset + j)
             !TODO THINK ABOUT FUSING THOSE TWO LOOPS. CAUSE... IT SUCKS!!!


             tmp = j*(j-1)/2+j
             do a = 1,virt_length
                a1 = virt_lam_allow(a)
                if (ignorable_pair(a1,orbital) ) cycle
                a1 = a1 - num_internal
                idum2 = a*(a-1)/2
                do b = 1, a
                   idum2 = idum2 + 1
                   b1 = virt_lam_allow(b) 
                   if (ignorable_pair(b1,orbital)) cycle
                   b1 = b1 - num_internal

                   Hab(idum2,tmp) = Hab(idum2,tmp) + aibi_debug(index2(a1,b1)+offset)*transposition_element
                enddo
             enddo

          enddo
          !$omp end do

       else

          !// HERE WE NEED IDENTITY MATRIX
          !$omp do schedule(static)
          do j = 1, lambda_dim
             tmp = j*(j-1)/2+j
             do a = 1,virt_length 
                a1 = virt_lam_allow(a) 
                if (ignorable_pair(a1,orbital) ) cycle
                a1 = a1 - num_internal
                idum2 = a*(a-1)/2
                do b = 1,a
                   idum2 = idum2 + 1
                   b1 = virt_lam_allow(b) 
                   if (ignorable_pair(b1,orbital) ) cycle
                   b1 = b1 - num_internal

                   Hab(idum2,tmp) = Hab(idum2,tmp) + aibi_debug(index2(a1,b1)+offset)

                enddo
             enddo
          enddo
          !$omp end do

       endif

    enddo

    !// ADD IN GIJ
    !$omp do schedule(static)
    do j = 1, lambda_dim
       tmp = j*(j-1)/2+j
       Hab(:,tmp) = Hab(:,tmp) + Gij_matrix(1:virt_length*(virt_length+1)/2)      
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine compute_Hab_vec_pao
  !*****************************************************************
  subroutine make_gij_ab_vec_pao(lambda_path,virt_lam_allow,twoVars,weight,virt_length,state)

    !// IN THIS SUBROUTINE WE MAKE THE CONTRIBUTION TO Gij FROM THE
    !// PORTION OF A PATH WHICH RESIDES ENTIRELY IN THE INTERNAL SPACE
    !// OR IS BEFORE THE FIRST LOOP SEGMENT.
    !// THE FORMULA IS:
    !//     Gij_internal = sum' { Np(ij|pp) - delta(Np,2)(ip|jp) }
    !// WHERE SUM IS OVER ORBITALS p AND THE PRIME
    !// MEANS DON'T INCLUDE ORBITALS i AND j IN THE SUM.  Np IS THE 
    !// OCCUPATION OF ORBITAL p.  UNLIKE MAKE_GIJ_INTERNAL, WE FORM
    !// THE MATRIX FOR ALL POSSIBLE A AND B.  THE RESULT GETS STORED
    !// IN GIJ_MATRIX, AND WE DO INCLUDE THE ONE ELECTRON INTEGRAL  

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use locist_var_mod,only:locist_scratch
    
    implicit none

    type(orbital_path),intent(in)::lambda_path
    integer,dimension(:),intent(inout)::virt_lam_allow
    type(twoModVars),intent(inout)::twoVars
    integer,intent(in)::virt_length                !/  LENGTH OF VIRTUALS FOR EACH INTERNAL CSF 
    integer,intent(in)::weight
    character,intent(in)::state 

    integer::path_level                 !// FOR LOOPING OVER LEVELS
    integer::a,b,ab_count               !// FOR LOOPING OVER A AND B
    integer::a1,b1                      !// FOR LOOPING A1 AND B1 
    integer::offset,offset2,idum,idum2

    real(real8),parameter::two = real(2.0, real8)                    !// NUMBER 2

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_virtuals(weight,state,virt_lam_allow)


    !// FIRST DO THE TWO ELECTRON INTEGRALS
    twoVars%Gij_matrix = real(0.0,real8)
    idum = num_internal*num_external*(num_external+1)/2

    do path_level = 1, num_internal

       if (lambda_path%occupations(path_level) == 1) then

          !          ab_count = 0
          offset = (path_level-1)*num_external*(num_external+1)/2

          do a = 1,virt_length 
             a1 = virt_lam_allow(a) - num_internal
             !             if (sep_orb(a1,path_level) == 0) cycle
             !             a1 = a1-num_internal

             idum2 = a*(a-1)/2
             do b = 1,a
                b1 = virt_lam_allow(b) - num_internal
                !                 if (sep_orb(b1,path_level) == 0) cycle
                !                 b1 = b1-num_internal

                ab_count = idum2 + b
                twoVars%Gij_matrix(ab_count) = twoVars%Gij_matrix(ab_count) +&
                     twoVars%aibi_debug(offset+index2(a1,b1))
             enddo
          enddo

       elseif (lambda_path%occupations(path_level) == 2) then

          !          ab_count = 0
          offset = (path_level-1)*num_external*(num_external+1)/2
          offset2 = idum + offset

          do a = 1,virt_length 
             a1 = virt_lam_allow(a) - num_internal 
             !             if (sep_orb(a1,path_level) == 0) cycle
             !             a1 = a1-num_internal

             idum2 = a*(a-1)/2
             do b = 1,a
                b1 = virt_lam_allow(b) - num_internal
                !                 if (sep_orb(b1,path_level) == 0) cycle
                !                 b1 = b1-num_internal

                ab_count = idum2+b
                twoVars%Gij_matrix(ab_count) = twoVars%Gij_matrix(ab_count) + &
                     two*twoVars%aibi_debug(index2(a1,b1)+offset) &
                     -twoVars%aibi_debug(index2(a1,b1)+offset2)

             enddo
          enddo

       endif

    enddo

    !// ADD IN THE ONE ELECTRON INTEGRAL
    ab_count = 0
    do a = 1,virt_length 
       a1 = virt_lam_allow(a) 
       do b = 1,a
          b1 = virt_lam_allow(b)
          ab_count = ab_count + 1
          twoVars%Gij_matrix(ab_count) = twoVars%Gij_matrix(ab_count) + pp(index2(a1,b1))
       enddo
    enddo


  end subroutine make_gij_ab_vec_pao
  !*****************************************************************

  !****************************************************************
  Subroutine one_intern_seg_compl_vec_lmo(civec, sigmavec, lambda_path,loc_scr)
  
    use locist_var_mod,only:locist_scratch

    !// IN THIS ROUTINE WE BUILD THE COMPLEMENT TO THE TWO SEGMENT LOOPS
    !// WHERE ONE SEGMENT RESIDES IN THE INTERNAL SPACE

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    type(orbital_path),intent(inout)::lambda_path
    type(locist_scratch),intent(inout)::loc_scr
    real(kind=real8), dimension(:) :: civec, sigmavec

    integer::ibar,jbar              !// ORBITAL POSITION INDICES
    integer::internal_singles       !// NUMBER OF SINGLES IN INTERNAL PATH OF PATH
    integer::allocatestatus         !// FOR DYNAMIC MEMORY
    integer::deallocatestatus
    integer::lambda_singles         !// SINGLES IN LAMBDA AND MU PATH
    integer::mu_singles
    integer::lambda_dim,mu_dim      !// DIMENSIONS FOR LAMBDA AND MU CONFIGURATIONS
    integer::lambda_weight          !// WEIGHTS OF LAMBDA AND MU PATH
    integer::mu_weight             
    integer::current_vertex         !// USED IN BUILDING MU PATH
    integer::step_type              !// KEEPS TRACK OF STEP IN MU PATH CONSTRUCTION
    integer::path_elecs             !// CUMULATIVE OCCUPATION IN MU PATH
    integer::top_level              !// LEVEL OF FIRST LOOP SEGMENT
    integer::mu_levels              !// LABELS LEVELS IN MU PATH
    integer::constraint_count       !// KEEP TRACK OF CONSTRAINTS IN BUILDING MU PATH
    integer::start_elec             !// OCCUPATION WHERE INTERNAL CSF MEETS UP WITH EXTERNAL SPACE
    integer::a                  !// USED TO INDEX EXTERNAL ORBITALS
    integer::s1,s2,s3,s4            !// FOR THE CALL TO GET NUMBERS OF SINGLES
    integer::loop_type              !// LABELS THE LOOP
    integer::i,j                    !// LABELS THE INTERNAL LOOP LEVELS
    integer::i_dum,j_dum,k_dum,l_dum


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
    integer::mode_1
    integer::case
    integer::common_virt_length

    integer::lam_length,mu_length,lam_lengthC2,mu_lengthC2

    real(real8)::zero, two                              !// THE NUMBERS ZERO AND TWO
    real(real8)::rdum
    real(real8),dimension(:,:,:),allocatable::Ka        !// PARTIAL MATRIX ELEMENTS
    real(real8),dimension(:,:),allocatable::Px          !// COUPLING COEFFICENTS
    real(real8),dimension(:),allocatable::Gij_vector    !// STORES INTERNAL GIJ        

    !////// EXTRA DECLARATIONS

    integer,dimension(num_orbitals)::common_virt,virt_lam_long,virt_lam_pos,virt_mu_long,virt_mu_pos,&
         common_virt_lam_pos,common_virt_mu_pos                 


    real(real8),dimension(:,:),allocatable:: tmp_storage        ! storing results from matmul
    real(real8),dimension(:,:,:),allocatable::tmp_ten1,tmp_ten2
    integer :: tmp_size              ! size of the storage
    integer :: tmp_size2

    real(real8),external::ddot
    
#ifdef TIGER_USE_OMP
    integer::threadID,numthreads
    threadID = OMP_get_thread_num()+1
    numthreads = numberOfThreads
#else
    integer,parameter::threadID = 1
    integer,parameter::numthreads = 1
#endif
    ! START AUTOGENERATED INITIALIZATION 
    ab_start = 0
    lam_length = 0
    lambda_singles = 0
    k_dum = 0
    lambda_weight = 0
    s4 = 0
    s3 = 0
    s1 = 0
    lam_lengthc2 = 0
    rdum = 0.0
    mu_length = 0
    a_end = 0
    ab_end = 0
    current_vertex = 0
    v_start = 0
    case = 0
    lambda_dim = 0
    j_dum = 0
    a_address = 0
    i_dum = 0
    common_virt_length = 0
    aa_start = 0
    a_start = 0
    mode_1 = 0
    ab_address = 0
    step_type = 0
    aa_end = 0
    a_spin = 0
    mu_lengthc2 = 0
    v_end = 0
    v_spin = 0
    mu_singles = 0
    s2 = 0
    mu_dim = 0
    mu_weight = 0
    ibar = 0
    a = 0
    l_dum = 0
    v_address = 0
    ! END AUTOGENERATED INITIALIZATION 


!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zero = real(0.0, real8)
    two = real(2.0, real8)
    
    loop_type = lambda_path%loop_type
    internal_singles = lambda_path%num_singles
    i = lambda_path%constraints(2,1)
    j = i
    start_elec = sum(lambda_path%occupations(0:num_internal))

    !// STORE THE WEIGHT OF THE LAMBDA_PATH
    lambda_path%weight = sum(lambda_path%arc_weights(0:num_internal))

    if (skip_this_internal(lambda_path%weight,start_elec)) return

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

    !// FOR THE RT_LOOP_WEIGHT ADD IN THE WEIGHTS OF THE RELEVANT PARTS
    !// OF THE HEAD AND TAIL PATHS.  THIS WILL END UP BEING THE MU PATH WEIGHT
    lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
         sum(lambda_path%arc_weights(0:lambda_path%level1-1))

    if (skip_this_internal(lambda_path%rt_loop_weight,path_elecs)) return

    !// ALLOCATE SCRATCH SPACE
    tmp_size = fsn(open_shells)
    allocate(Px(tmp_size,tmp_size),stat=allocatestatus)
    call allocatecheck(allocatestatus,"one_intern_seg_compl_vec_lmo")
    !Px = 0

!!!!!!!!!!!!!!!
    !//
    !// GET JBAR
    !//
!!!!!!!!!!!!!!!                             
    call get_numbers_of_singles(s1,s2,s3,s4,lambda_path)
    jbar = s1+1
    if(loop_type == 2) jbar = s1

!!!!!!!!!!!!!!!
    !//
    !//  FORM EXTERNAL CONTRIBUTIONS
    !//  AND DO MULTIPLICATIONS
    !//
!!!!!!!!!!!!!!!                             
    if (start_elec == num_elec-1) then    

       !// HERE WE HANDLE THE FOLLOWING TERMS (a|G|i)A(Px) + Ka  DESCRIBING INTERACTIONS 
       !/ /BETWEEN VALENCE AND N-1 STATES AND VICE VERSA. 
       !// THE COUPLING COEFFICIENTS FOR LOOP_TYPE 1 IS A1 AND FOR LOOP_TYPE 2 IS A2. 
       !// THE VARIABLE CASE IS USED FOR COMPUTING Ka,ADD_GIJ_VEC         

       if (loop_type ==1) then 
          case = 1
       elseif (loop_type == 2) then 
          case = 3
       endif

       ibar = internal_singles + 1
       call one_info_simplifier(loop_type,start_elec,internal_singles,ibar,jbar, & 
            lambda_singles,mu_singles, lambda_dim, mu_dim, & 
            lambda_path%weight+1, lambda_path%rt_loop_weight+1, & 
            v_address,a_address,ab_address,aa_address,Px)

!!!!!!!!!!!!!!!!!!!
       !//               
       !// BUILD {13} and {53} LOOPS 
       !//               
!!!!!!!!!!!!!!!!!!!

       if (lambda_dim /=0 .and. mu_dim /=0 ) then 

          if (a_address > 0 .and. v_address > 0) then
          
             if(internal_index_vector1(lambda_path%weight+1) >= 0) then

             !// MAKE GIJ FOR ALL A: 

             lambda_weight = lambda_path%weight+1 
             lam_length  = num_allowed_virtuals(lambda_weight,"S")

             tmp_size2 = max(lambda_dim,mu_dim)
             allocate(Ka(tmp_size2,tmp_size2,lam_length),Gij_vector(lam_length),tmp_storage(lambda_dim,mu_dim),stat=allocatestatus)
             call allocatecheck(allocatestatus,"one_intern_seg_compl_vec_lmo II")
             call make_gij_a_pao(lambda_path,loc_scr%virt_lam_allow,i,Gij_vector,loop_type,lambda_weight,lam_length,"S")

             !// MAKE Ka . THIS DEPENDS ON THE VALUE OF CASE WHICH DEPENDS ON THE LOOP_TYPE 
             call compute_Ka_vec_pao(lambda_path,loc_scr%virt_lam_allow,i,case,&
                  ibar, jbar, &
                  lambda_singles, mu_singles, &
                  Ka,lambda_weight,lam_length,"S")

             !// ADD Gij TO DIAGONAL ELEMENTS OF Ka
             call add_in_Gij_vec(Ka,lambda_dim,Gij_vector,case,lam_length)

             !// DO THE MATMUL OF COUPLING COEFFICIENTS WITH KA TO OBTAIN MATRIX ELEMENTS
             do a=1,lam_length
                call dgemm("N","N",lambda_dim,mu_dim,lambda_dim,1.0,Ka(1,1,a), &
                     tmp_size2,Px,tmp_size,0.0,tmp_storage,lambda_dim)
                Ka(1:lambda_dim,1:mu_dim,a) = tmp_storage(:,:)
             enddo
             
             deallocate(Gij_vector,tmp_storage,stat=deallocatestatus)
             call deallocatecheck(deallocatestatus,"one_intern_seg_compl_vec_lmo II")

             !// HIT THE CI VECTOR WITH THESE MATRIX ELEMENTS 
             !// (THESE INTERACTIONS ARE SIMILIAR TO (IJ|KA) ROUTINES)

             do a_spin =1,lambda_dim

                a_start = a_address + lam_length*(a_spin - 1)
                a_end   = a_start   + lam_length- 1

                do v_spin =1,mu_dim 

                   v_start = v_address + v_spin -1
                   v_end   = v_start 

                   !// NOW WE HAVE THE SIGMAVECTOR OVER PAO BASIS (VALENCE)

                   sigmavec(v_start) = sigmavec(v_start) +&
                        ddot(lam_length,Ka(a_spin,v_spin,1),tmp_size2*tmp_size2,civec(a_start),1)

                   !// NOW WE HAVE THE SIGMAVECTOR OVER PAO BASIS (N-1)

                   sigmavec(a_start:a_end) = sigmavec(a_start:a_end) + &
                        Ka(a_spin,v_spin,1:lam_length)*civec(v_start)

                enddo
             enddo
             
             deallocate(Ka,stat=deallocatestatus)
             call deallocatecheck(deallocatestatus,"one_intern_seg_compl_vec_lmo III")
             
             endif

          endif

       endif

    elseif (start_elec == num_elec-2) then   
       !  goto 200 

       !// HERE WE HANDLE THE FOLLOWING TERMS (a|G|i)A(Px) + Ka  DESCRIBING INTERACTIONS 
       !/ /BETWEEN N-2 AND N-1 STATES AND VICE VERSA. PLEASE SEE THE CPR PAPER BY DUCH AND 
       !// KARWOWSKI FOR MORE DETAILS.
       !// THE COUPLING COEFFICIENTS FOR LOOP_TYPE 1 IS A1 AND FOR LOOP_TYPE 2 IS A2. 
       !// THE VARIABLE CASE IS USED FOR COMPUTING Ka,ADD_GIJ_VEC         
       !//
       !// Daa and Sa (a|G|i)A(Px) + Ka  (The term (ia|aa)*A(Px) is incorporated in the
       !//                                (IA|BC) routine) 
       !// Dba and Sb (a|G|i)A(Px) + Ka  (The terms (ai|bb)A(Px) + (ab|ib)A(Py(s)) are
       !//                                incorporated in (IA|BC) routine) 

       if (loop_type == 1) then 
          case = 1
       elseif (loop_type == 2) then 
          case = 3
       endif

       ibar = internal_singles + 2 

       call one_info_simplifier(loop_type,start_elec,internal_singles,ibar,jbar, & 
            lambda_singles,mu_singles, lambda_dim, mu_dim, & 
            lambda_path%weight+1, lambda_path%rt_loop_weight+1, & 
            v_address,a_address,ab_address,aa_address,Px)

!!!!!!!!!!!!!!!!!!!!!
       !//                 
       !// BUILD {13},{35} 
       !//       {53},{75} LOOPS  
       !//                 
!!!!!!!!!!!!!!!!!!!!!

       if (lambda_dim /= 0 .and. mu_dim /= 0) then 

          if ((ab_address > 0 .or. aa_address > 0) .and. a_address > 0) then
          
             if (internal_index_vector2(lambda_path%weight+1) >= 0 .and. internal_index_vector1(lambda_path%rt_loop_weight+1) >= 0) then

             lambda_weight = lambda_path%weight+1 
             lam_length    = num_allowed_virtuals(lambda_weight,"D")
             lam_lengthC2 = lam_length*(lam_length-1)/2

             mu_weight = lambda_path%rt_loop_weight+1 
             mu_length = num_allowed_virtuals(mu_weight,"S")
             mu_lengthC2 = mu_length*(mu_length-1)/2


             ! Don't change i_dum,j_dum to i,j. The latter are used earlier!
             call get_virtuals(lambda_weight,"D",loc_scr%virt_lam_allow)
             call get_virtuals(mu_weight,"S",loc_scr%virt_mu_allow)

             !// MAKE GIJ FOR ALL A:
             tmp_size2 = max(lambda_dim,mu_dim)
             allocate(Ka(tmp_size2,tmp_size2,lam_length),Gij_vector(lam_length),tmp_storage(lambda_dim,mu_dim),stat=allocatestatus)
             call allocatecheck(allocatestatus,"one_intern_seg_compl_vec_lmo II")
             call make_gij_a_pao(lambda_path,loc_scr%virt_lam_allow,i,Gij_vector,loop_type,lambda_weight,lam_length,"D")

             !// MAKE Ka . THIS DEPENDS ON THE VALUE OF CASE WHICH DEPENDS ON THE LOOP_TYPE 
             call compute_Ka_vec_pao(lambda_path,loc_scr%virt_lam_allow,i,case,&
                  ibar, jbar, &
                  lambda_singles, mu_singles, &
                  Ka,lambda_weight,lam_length,"D")

             !// ADD Gij TO DIAGONAL ELEMENTS OF Ka
             call add_in_Gij_vec(Ka, lambda_dim, Gij_vector,case,lam_length)

             !// DO THE MATMUL OF COUPLING COEFFICIENTS WITH KA TO OBTAIN MATRIX ELEMENTS
             do a=1,lam_length
                call dgemm("N","N",lambda_dim,mu_dim,lambda_dim,1.0,Ka(1,1,a), &
                     tmp_size2,Px,tmp_size,0.0,tmp_storage,lambda_dim)
                Ka(1:lambda_dim,1:mu_dim,a) = tmp_storage(:,:)
             enddo
             
             deallocate(Gij_vector,tmp_storage,stat=deallocatestatus)
             call deallocatecheck(deallocatestatus,"one_intern_seg_compl_vec_lmo II")

             !// CALL A ROUTINE TO FETCH THE PAOVLP MATRIX 

             ! Making Changes for the lmo version

             common_virt = 0
             virt_lam_long = 0
             virt_lam_pos = 0
             virt_mu_long = 0
             virt_mu_pos = 0

             ! Keep the virt_allowed in long format and keep track of position
             do i_dum = 1, lam_length
                j_dum = loc_scr%virt_lam_allow(i_dum)
                virt_lam_long(j_dum) = j_dum
                virt_lam_pos(j_dum) = i_dum
             enddo

             do i_dum = 1, mu_length
                j_dum = loc_scr%virt_mu_allow(i_dum)
                virt_mu_long(j_dum) = j_dum
                virt_mu_pos(j_dum) = i_dum
             enddo

             ! Get the common virtuals between lam,mu and also the reverse tracking position

             common_virt_length = 0
             common_virt_lam_pos = 0
             common_virt_mu_pos = 0
             do i_dum = 1, lam_length
                j_dum = loc_scr%virt_lam_allow(i_dum)
                if ( virt_lam_long(j_dum) == virt_mu_long(j_dum) ) then
                   common_virt_length = common_virt_length + 1
                   common_virt(common_virt_length) = j_dum
                   common_virt_lam_pos(common_virt_length) = virt_lam_pos(j_dum)
                   common_virt_mu_pos(common_virt_length) = virt_mu_pos(j_dum)
                endif
             enddo
             
             !if (internal_index_vector2(lambda_path%weight+1) >= 0 .and. internal_index_vector1(lambda_path%rt_loop_weight+1) >= 0) then
             
             ! allocate some scratch arrays
             allocate(tmp_ten1(lam_length,lam_length,numthreads),tmp_ten2(lam_length,lam_length,numthreads),stat=allocatestatus)
             call allocatecheck(allocatestatus,"tmp_ten12")
             
             !// HIT THE CI VECTOR WITH THESE MATRIX ELEMENTS 
             !// (THESE INTERACTIONS ARE SIMILIAR TO (IJ|KA) ROUTINES)
             !b$bomp parallel do &
             !b$bomp schedule(static) &
             !b$bomp default(none) &
             !b$bomp private(ab_start,ab_end,aa_start,aa_end,threadID,a_start,a_end,i_dum,l_dum,rdum) &
             !b$bomp shared(lambda_dim,fsn,lam_lengthC2,lam_length,tmp_ten1,civec,internal_index_vector2, &
             !b$bomp lambda_path,mu_dim,a_address,mu_length,internal_index_vector1,common_virt_length, &
             !b$bomp ab_address,aa_address,lambda_singles,common_virt_mu_pos,common_virt_lam_pos,sigmavec, &
             !b$bomp Ka, zero,tmp_ten2,mode_1)
             do ab_spin = 1,lambda_dim

#ifdef TIGER_USE_OMP
                threadID = OMP_get_thread_num()+1
#endif
             
                ab_start = ab_address + lam_lengthC2*(ab_spin-1)
                aa_start = aa_address + lam_length*(ab_spin-1)

                mode_1 = -2
                if (ab_spin <= fsn(lambda_singles-2)) mode_1 = 1

                if (mode_1 .eq. 1) then
                    call scepper_diag_p1(tmp_ten1(:,:,threadID),civec,ab_start, &
                     lam_length,civec,aa_start)
                elseif (mode_1 .eq. -2) then
                    call scepper_diag_m2(tmp_ten1(:,:,threadID),civec,ab_start, &
                     lam_length)
                endif
                     
                tmp_ten2(:,:,threadID)= zero
                do a_spin = 1,mu_dim 

                   a_start = a_address + mu_length*(a_spin - 1)


                   !if (internal_index_vector1(lambda_path%rt_loop_weight+1) < 0) cycle

                   !// THE CIVECTOR IS ALREADY IN PAO BASIS 
                   !///////////////////////////////// For the LMO version //////////////////////////////////////////////////////////////////

                   do k_dum = 1, common_virt_length
                      i_dum = common_virt_mu_pos(k_dum)
                      l_dum = common_virt_lam_pos(k_dum)
                      sigmavec(a_start+i_dum-1) = sigmavec(a_start+i_dum-1) + &
                           ddot(lam_length,Ka(ab_spin,a_spin,1),tmp_size2*tmp_size2,tmp_ten1(1,l_dum,threadID),1)
                   enddo

                   !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                   !///////////////////////////////// For the LMO version //////////////////////////////////////////////////////////////////


                   do j_dum = 1,lam_length 
                      rdum = Ka(ab_spin,a_spin,j_dum)
                      do i_dum = 1,common_virt_length
                         k_dum = common_virt_lam_pos(i_dum)  
                         l_dum = common_virt_mu_pos(i_dum)
                         tmp_ten2(k_dum,j_dum,threadID) = tmp_ten2(k_dum,j_dum,threadID) + rdum*civec(a_start+l_dum-1)
                      enddo
                   enddo

                   !////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                enddo

                if (mode_1 .eq. 1) then
                    call unscepper_diag_p1(tmp_ten2(:,:,threadID),sigmavec,ab_start,lam_length,sigmavec,aa_start)
                elseif (mode_1 .eq. -2) then
                    call unscepper_diag_m12(tmp_ten2(:,:,threadID),sigmavec,ab_start,lam_length)
                endif    
            
             enddo
             !b$bomp end parallel do
             
             ! deallocate our scratch arrays
             deallocate(tmp_ten1,tmp_ten2,Ka,stat=deallocatestatus)
             call deallocatecheck(deallocatestatus,"one_intern_seg_compl_vec_lmo III")
             
             endif

          endif

       endif

       !  200 continue
    endif

    !// DEALLOCATE MATRIX ELEMENT STORAGE ARRAY

    deallocate(Px,stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"one_intern_seg_compl_vec_lmo")

  end subroutine one_intern_seg_compl_vec_lmo

  !*****************************************************************
  subroutine make_gij_a_pao(lambda_path,virt_lam_allow,level1,Gij_internal,loop_type,weight,virt_length,state)

    !// IN THIS SUBROUTINE WE MAKE THE CONTRIBUTION TO Gij FROM THE
    !// PORTION OF A PATH WHICH RESIDES ENTIRELY IN THE INTERNAL SPACE
    !// OR IS BEFORE THE FIRST LOOP SEGMENT.
    !// THE FORMULA IS:
    !//     Gij_internal = sum' { Np(ij|pp) - delta(Np,2)(ip|jp) }
    !// WHERE SUM IS OVER ORBITALS p AND THE PRIME
    !// MEANS DON'T INCLUDE ORBITALS i AND j IN THE SUM.  Np IS THE 
    !// OCCUPATION OF ORBITAL p.  UNLIKE MAKE_GIJ_INTERNAL, WE FORM
    !// THE MATRIX FOR ALL POSSIBLE A.  THE RESULT GETS STORED
    !// IN GIJ_MATRIX, AND WE DO INCLUDE THE ONE ELECTRON INTEGRAL  

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use locist_var_mod,only:locist_scratch
    
    implicit none

    type(orbital_path),intent(in)::lambda_path
    integer,dimension(:),intent(inout)::virt_lam_allow
    integer,intent(in)::level1
    integer,intent(in)::loop_type
    integer,intent(in)::weight
    integer,intent(in)::virt_length                !// LENGTH OF VIRTUALS 
    character,intent(in)::state                    !// "S" or "s" for N-1 and "D" or "d" for N-2
    real(real8),dimension(:),intent(inout)::Gij_internal
    
    integer::path_level                 !// FOR LOOPING OVER LEVELS
    integer::integral_offset  
    integer::a1,a

    real(real8),parameter::two = real(2.0, real8)                    !// NUMBER 2

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// GET THE allowed number of virtuals 
    call get_virtuals(weight,state,virt_lam_allow)

    !// FIRST DO THE TWO ELECTRON INTEGRALS
    Gij_internal = real(0.0,real8)

    do path_level = 1, num_internal

       if (path_level == level1) cycle 

       !      if (sep_orb(path_level,i_ind) == 0) cycle

       if (lambda_path%occupations(path_level) == 1) then

          do a = 1,virt_length 
             a1 = virt_lam_allow(a)

             integral_offset = (a1-num_internal-1)*num_orbitals
             Gij_internal(a) = Gij_internal(a) +&
                  integral_buffer(1,path_level+integral_offset)
          enddo

       elseif (lambda_path%occupations(path_level) == 2) then

          do a = 1,virt_length 
             a1 = virt_lam_allow(a)

             integral_offset = (a1-num_internal-1)*num_orbitals
             Gij_internal(a) = Gij_internal(a) +&
                  two*integral_buffer(1,path_level+integral_offset) &
                  -integral_buffer(2,path_level+integral_offset)
          enddo

       endif

    enddo

    !// ADD IN THE ONE ELECTRON INTEGRAL
    if (loop_type == 1) then 

       do a = 1,virt_length 
          a1 = virt_lam_allow(a)
          Gij_internal(a) = Gij_internal(a) + pp(index2(level1,a1))
       enddo

    else 

       do a = 1,virt_length 
          a1 = virt_lam_allow(a)
          integral_offset = (a1-num_internal-1)*num_orbitals
          Gij_internal(a) = Gij_internal(a)  + & 
               integral_buffer(1,integral_offset+level1) + & 
               pp(index2(level1,a1))   
       enddo

    endif

  end subroutine make_gij_a_pao

  !*****************************************************************
  subroutine compute_Ka_vec_pao(lambda_path, virt_lam_allow,level1, loop_type,&
       ibar, jbar, lambda_singles, mu_singles,&
       Ka,weight,virt_length,state)

    !// FOR THE {13}, {35}, {53}, AND {57} LOOPS THE MATRIX ELEMENTS ALL
    !// INVOLVE A SUM OVER SINGLY OCCUPIED ORBITALS WHERE YOU
    !// MULTIPLY AN INTEGRAL BY A TRANSPOSITION.  HERE WE TREAT THIS SUM.
    !// THE FORMULAS ARE AS FOLLOWS:
    !// 
    !//    {13} = SUM { DELTA(Np,1) (ibar, pbar) (IP|JP) }
    !//    {35} = SUM { DELTA(Np,1) (jbar, pbar) (IP|JP) }
    !//    {53} = SUM { DELTA(Np,1) (ibar, pbar) (IP|JP) }
    !//    {57} = SUM { DELTA(Np,1) ((ibar, pbar) + 1) (IP|JP) }

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use locist_var_mod,only:locist_scratch

    implicit none

    type(orbital_path),intent(in)::lambda_path
    integer,dimension(:),intent(inout)::virt_lam_allow

    integer,intent(in)::level1                         !// STARTING LEVEL FOR SUM
    integer,intent(in)::loop_type                      !// TYPE OF THE LOOP
    integer,intent(in)::ibar,jbar                      !// POSITION INDICES
    integer,intent(in)::lambda_singles                 !// NUMBERS OF SINGLES
    integer,intent(in)::mu_singles
    integer,intent(in)::weight
    integer,intent(in)::virt_length                     !// LENGTH OF VIRTUALS 
    character,intent(in)::state                         !// "S" or "s" for N-1 and "D" or "d" for N-2
    real(real8), dimension(:,:,:),intent(inout)::Ka   !// SEE FORMULAS IN INTRO
    
    integer::j,k                          !// USED FOR LOOPING OVER LEVELS
    integer::pbar                           !// KEEPS TRACK OF POSITION INDICES
    integer::lambda_dim, mu_dim             !// DIMENSIONS OF MATRICES
    !integer::length                         !// SIZE OF TRANSPOSITION SUM
    integer::transposition_label            !// LABELS THE TRANSPOSITION SUM
    integer::row_offset
    integer::orbital
    integer::offset                         !// FOR INTEGRAL ACCESS
    integer::a     
    integer::a1

    real(real8)::the_integral               !// POINTER TO A LPVOID BUNGHOLE CLASS
    real(real8)::transposition_element      !// ELEMENT OF A TRASPOSITION MATRIX
    real(real8),parameter::zero = real(0.0,real8)

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pbar = 0
    lambda_dim = fsn(lambda_singles)
    mu_dim = fsn(mu_singles)

    Ka = zero

    !// GET THE allowed number of virtuals 
    call get_virtuals(weight,state,virt_lam_allow)
!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// LOOP OVER SINGLY OCCUPIED ORBITALS IN THE INTERNAL SPACE
    do pbar = 1, lambda_path%num_singles

       !// GET THE ORBITAL LABEL
       orbital = lambda_path%singles(pbar)

       !// CYCLE IF THIS IS ONE OF THE LOOP LEVELS
       if (orbital == level1) cycle

       !if (sep_orb(orbital,twoVars%i_ind) == 0) cycle

!!!!!!!!!!!!
       !//
       !// {35} LOOPS NEED (JBAR, PBAR)
       !//
!!!!!!!!!!!!
       if (loop_type == 2) then

          if (pbar /= jbar) then

             transposition_label = index2m1(pbar,jbar)

             do a = 1,virt_length
                a1 = virt_lam_allow(a)
                offset = (a1-num_internal-1)*num_orbitals
                the_integral = integral_buffer(2, offset+orbital)

                do j = 1, lambda_dim
                   row_offset = j*(j-1)/2
                   do k = 1, j-1
                      transposition_element = &
                           transpositions(transposition_label, row_offset + k)
                      Ka(j,k,a) = Ka(j,k,a) + &
                           the_integral*transposition_element
                      Ka(k,j,a) = Ka(k,j,a) + &
                           the_integral*transposition_element
                   enddo

                   transposition_element = &
                        transpositions(transposition_label, row_offset + j)

                   Ka(j,j,a) = Ka(j,j,a) + &
                        the_integral*transposition_element
                enddo
             enddo


          else

             !// HERE WE NEED IDENTITY MATRIX
             do a = 1,virt_length
                a1 = virt_lam_allow(a)
                offset = (a1-num_internal-1)*num_orbitals
                the_integral = integral_buffer(2, offset+orbital)
                do j = 1, lambda_dim
                   Ka(j,j,a) = Ka(j,j,a) + the_integral
                enddo
             enddo

          endif

!!!!!!!!!!!!!!!!!!!! 
          !//
          !// HERE WE NEED (ibar, pbar) TO TREAT EVERYTHING BUT
          !// THE {35} LOOPS.   IN THE CASE OF THE {57} LOOPS
          !// WE WILL ADD IN THE IDENTITY AFTER 
          !//
!!!!!!!!!!!!!!!!!!!!
       else

          if (pbar /= ibar) then

             transposition_label = index2m1(pbar, ibar)

             do a = 1,virt_length
                a1 = virt_lam_allow(a)
                offset = (a1-num_internal-1)*num_orbitals
                the_integral = integral_buffer(2, offset+orbital)
                do j = 1, lambda_dim
                   row_offset = j*(j-1)/2
                   do k = 1, j-1
                      transposition_element = &
                           transpositions(transposition_label, row_offset + k)
                      Ka(j,k,a) = Ka(j,k,a) + &
                           the_integral*transposition_element
                      Ka(k,j,a) = Ka(k,j,a) + &
                           the_integral*transposition_element
                   enddo

                   transposition_element = &
                        transpositions(transposition_label, row_offset + j)

                   Ka(j,j,a) = Ka(j,j,a) + &
                        the_integral*transposition_element
                enddo
             enddo


          else

             !// HERE WE NEED IDENTITY MATRIX
             do a = 1,virt_length
                a1 = virt_lam_allow(a)
                offset = (a1-num_internal-1)*num_orbitals
                the_integral = integral_buffer(2, offset+orbital)
                do j = 1, lambda_dim
                   Ka(j,j,a) = Ka(j,j,a)+the_integral
                enddo
             enddo

          endif

       endif

!!!!!!!!!!!!!!!!
       !//
       !// ADD IN THE IDENTITY FOR LOOP {57}
       !//
!!!!!!!!!!!!!!!!
       if (loop_type == 4) then

          do a = 1,virt_length
             a1 = virt_lam_allow(a)
             offset = (a1-num_internal-1)*num_orbitals
             the_integral = integral_buffer(2, offset+orbital)
             do j = 1, lambda_dim
                Ka(j,j,a) = Ka(j,j,a)+the_integral
             enddo
          enddo

       endif

    enddo

  end subroutine compute_Ka_vec_pao

end module two_seg_mod_3

