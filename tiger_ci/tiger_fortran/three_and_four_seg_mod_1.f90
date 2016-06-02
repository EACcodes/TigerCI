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
!>      WE TREAT THE FOLLOWING CONTRIBUTIONS HERE:
!>        - PURELY INTERNAL THREE SEGMENT LOOPS
!>        - THREE SEGMENT LOOPS WITH J AND K SEGMENTS IN INTERNAL SPACE
!>        - THREE AND FOUR SEGMENT LOOPS WITH ONE SEGMENT IN THE INTERNAL SPACE
!>        - PURELY EXTERNAL THREE AND FOUR SEGMENT LOOPS
!>
!> \author Derek Walter
!> \author Arun Venkatnathan
!*****************************************************************

! To reduce symbol size the following abbreviations are made
! David Krisiloff 12/22/2010
!
! coefficients --> coeff
! coupling     --> coupl
! external     --> extern
! internal     --> intern
! complement   --> compl 

module three_and_four_seg_mod_1

  use graph_var_mod
  use global_var_mod
  use molecule_var_mod
  use integral_storage_mod
  use get_integrals_mod
  use ci_utilities_mod
  use decider_symbols
  use new_tree_search_mod
  use tree_search_mod
  use locist_mod
  use scepmaker_mod
  use three_four_seg_var_mod
  use abcd_seg_info_replacement
  use iabc_seg_info_replacement
  use spin_mod
  use spin_var_mod
  use pseudo_seg_info_mod
  use sgga_helpers
  use blocked_locks_mod
  use fortran_timing
  use IOBuffer
  use make_all_integrals
  
  implicit none
  
  private
  
#ifdef TIGER_USE_OMP
  public :: three_and_four_seg_driver_1, iabc_sig_cho_big,extern_one_seg_compl_vec_lmo
#else
  public :: three_and_four_seg_driver_1, iabc_sig_cho_big,extern_one_seg_compl_vec_lmo
#endif
  
#ifdef TIGER_USE_OMP
  !// THESE ARRAYS ARE *ONLY* USED TO COMMUNICATE INFORMATION (OFFSETS ETC)
  !// IN THE OPENMP CASE WITHIN THIS MODULE. DO NOT USE FOR ANYTHING ELSE!
  !// ALSO, DO NOT PUT ANY OTHER GLOBAL MODULE VARIABLES HERE WITHOUT 
  !// PROPER JUSTIFICATION!
  integer,dimension(:),allocatable,private::omp_offsets_abcd_rec
  integer,dimension(:),allocatable,private::omp_offsets_iabc_rec
#endif
  
  !// WE HAVE TO ALLOCATE AN ARRAY THAT CONTAINS THE LOOP INFORMATION.  
  !// I HAVE WRITTEN IT IN A NICE EASY TO UNDERSTAND FORMAT.  THE FIRST
  !// PATH IS THE LAMBDA PATH AND THE SECOND PATH IS THE MU PATH.  
  !// BELOW WE HAVE ONLY LISTED THE THREE SEGMENT LOOPS.  THE INFORMATION
  !// FOR THE FOUR SEGMENT LOOPS AND FRAGMENTS OF FOUR SEGMENT LOOPS IS
  !// CONTAINED IN A SEPARATE MODULE.  THIS HAS TO DO WITH THE DEVELOPMENT
  !// PATH OF THIS PROGRAM.  
  integer, private, parameter, dimension(12,2,3)::loops=&
       reshape((/0,1,2,&    !// B1 CHAIN STARTS HERE
       1,2,0,&    !// {156}
       
       2,1,0,&
       1,0,2,&    !// {732}
       
       0,2,1,&
       1,0,2,&    !// {165}
       
       2,0,1,&
       1,2,0,&    !// {723}
       
       2,0,1,&
       0,1,2,&    !// {615}
       
       0,2,1,&
       2,1,0,&    !// {273}
       
       1,1,0,&    !// B2 CHAIN STARTS HERE
       0,0,2,&    !// {332}
       
       1,1,2,&
       2,2,0,&    !// {556}
       
       1,0,1,&
       0,2,0,&    !// {323}
       
       1,2,1,&
       2,0,2,&    !// {565}
       
       0,1,1,&
       2,0,0,&    !// {233}
       
       2,1,1,&
       0,2,2/),&  !// {655}
       shape = (/12,2,3/),&
       order = (/3,2,1/))
       
  type threeExtScr
    integer,dimension(:),allocatable::virt_lam_allow,virt_mu_allow
    real(real8),dimension(:,:),allocatable::Px
    real(real8),dimension(:,:,:),allocatable::iabc_buffer
    real(real8),dimension(:,:),allocatable::ia_mo_vecs,nonred_iabc_buffer
  end type threeExtScr
       
  type fourExtResizeScr
      real(real8),dimension(:),allocatable::integrals_abcd_mod2,integrals_abcd_mod3,integrals_abcd_mod4
      real(real8),dimension(:),allocatable::comprSigma
      integer::length2,length2C2,lengthSigma
  end type fourExtResizeScr
       
  type fourExtScr
    real(real8),dimension(:,:),allocatable::integrals_abcd
    integer,dimension(:),allocatable::the_virtuals
    type(fourExtResizeScr)::resizeScr
  end type fourExtScr

  contains


  !*****************************************************************
  subroutine three_and_four_seg_driver_1(civec, sigmavec, iteration,cho_data,loc_scr,mod1vars)

    !// THIS IS THE MAIN DRIVER ROUTINE

    use time_var_mod
    use utilities_mod
    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use two_seg_var_mod

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::iteration
    integer::allocatestatus    !// FOR DYNAMIC MEMORY ALLOCATION
    integer::deallocatestatus

    type(orbital_path)::lambda_path
    type(cholesky_data)::cho_data
    type(locist_scratch)::loc_scr
    type(threefourmod1vars)::mod1vars
    type(clock) :: timer
    type(blockedLockVectorType) :: civec, sigmavec

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(lambda_path%arc_weights(0:num_orbitals), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%2sg_wei")

    allocate(lambda_path%occupations(0:num_orbitals), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%lam_occ")

    allocate(lambda_path%singles(0:num_orbitals), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%lam_sin")

    allocate(lambda_path%constraints(3,6), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%lam_con")

    allocate(lambda_path%encountered(3), stat=allocatestatus)
    call allocatecheck(allocatestatus,"%lam_enc")
    lambda_path%arc_weights = 0
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%encountered = .false.
    lambda_path%constraints = 0

    allocate(mod1vars%scep_ci_1(num_external,num_external),&
         mod1vars%scep_sigma_1(num_external,num_external),&
         stat=allocatestatus)
    call allocatecheck(allocatestatus,"scepper ")
    mod1vars%scep_ci_1 = 0
    mod1vars%scep_sigma_1 = 0

    allocate(mod1vars%integral_buf(num_external,num_external),stat=allocatestatus)
    call allocatecheck(allocatestatus,"int     ")
    mod1vars%integral_buf = 0

    allocate(mod1vars%Px(fsn(open_shells),fsn(open_shells)),&
         stat=allocatestatus)
    call allocatecheck(allocatestatus,"PxPy2   ")
    mod1vars%Px = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!
    !//
    !// PURELY INTERNAL 3 SEGMENT
    !// LOOPS
    !//
!!!!!!!!!!!!!!!!!!!!!!!!

    mod1vars%iter_count = iteration

    if (num_internal >= 3) then
       call start_clock(timer)
       call all_i_3_seg_loopdrv(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)
       ikjk_time = ikjk_time + get_clock_wall_time(timer)
       call print_clock(timer, "  + three internal")
       call flush(ioOutput)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!
    !//
    !// THREE AND FOUR SEGMENT LOOPS
    !// WITH 1 SEGMENT IN THE 
    !// INTERNAL SPACE
    !//
!!!!!!!!!!!!!!!!!!!!!!!!  
    call start_clock(timer)

    if (valence_ci_flag .ne. 1 .and. reference_ci_flag .ne. 1) then
       call one_internal_seg_vec_2_cho(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)
    endif

    if (valence_ci_flag .ne. 1 .and. reference_ci_flag .ne. 1) then
       call one_internal_seg_vec_1_cho(cho_data, civec, sigmavec, lambda_path,mod1vars) ! faster
    endif
    
    call print_clock(timer, "  + three external")
    iabc_time = iabc_time + get_clock_wall_time(timer)
    call flush(ioOutput)

!!!!!!!!!!!!!!!!!!!!!!!!
    !//
    !// THREE AND FOUR SEGMENT LOOPS
    !// WHICH RESIDE ENTIRELY IN
    !// THE EXTERNAL SPACE
    !//
!!!!!!!!!!!!!!!!!!!!!!!!      

    call start_clock(timer)

    if (valence_ci_flag .ne. 1 .and. reference_ci_flag .ne. 1) then
       call all_external_3_4_loops_cho_2(cho_data,civec, sigmavec, lambda_path,iteration)
    endif

    call print_clock(timer, "  + four external")
    abcd_time = abcd_time + get_clock_wall_time(timer)
    call flush(ioOutput)

    !// CLEAN UP BEFORE QUITTING.
    deallocate(lambda_path%arc_weights, stat = deallocatestatus)
    call deallocatecheck(deallocatestatus, "%2la_wei")

    deallocate(lambda_path%occupations, stat = deallocatestatus)
    call deallocatecheck(deallocatestatus, "%2la_occ")

    deallocate(lambda_path%singles, stat = deallocatestatus)
    call deallocatecheck(deallocatestatus, "%2la_sin")

    deallocate(lambda_path%constraints, stat = deallocatestatus)
    call deallocatecheck(deallocatestatus, "%2la_con")

    deallocate(lambda_path%encountered, stat = deallocatestatus)
    call deallocatecheck(deallocatestatus, "encounte")

    if (allocated(integral_buffer)) then
       deallocate(integral_buffer, stat = deallocatestatus)
       call deallocatecheck(deallocatestatus, "integral")
    endif

    if (allocated(iabc_buffer)) then
       deallocate(iabc_buffer, stat = deallocatestatus)
       call deallocatecheck(deallocatestatus, "iabc_buffer")
    endif
    
    call deallocMod1Vars(mod1vars)

  end subroutine three_and_four_seg_driver_1

  !*****************************************************************
  !*****************************************************************
  subroutine all_external_3_4_loops_cho_2(cho_data,civec, sigmavec, lambda_path,iteration)

    !// THIS ROUTINE TREATS THE LOOPS WHICH RESIDE ENTIRELY IN THE EXTERNAL
    !// SPACE.  THE TREATMENT HERE WILL INCLUDE FULL SEPARATION OF THE 
    !// EXTERNAL SPACE.  BECAUSE OF THE WAY WE DO THINGS, WE FIRST JUST NEED
    !// TO CALL A ROUTINE TO DO THE TREE SEARCH.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use utilities_mod
    use global_var_mod
    use io_unit_numbers
    use cholesky_structs
    use locist_var_mod
    use three_four_seg_var_mod
    use two_seg_var_mod

    implicit none
    
    integer,intent(in)::iteration
    type(orbital_path),intent(inout)::lambda_path
    type(blockedLockVectorType) :: civec, sigmavec
    type(cholesky_data),intent(in)::cho_data

    integer::allocatestatus
    integer::a,c                   !// FOR LOADING UP (AC|BC) IN PSEUDOSPECTRAL MODE
    integer::icount
  
    integer::abcd_rec
    type(graph_search_state) :: graph
    type(fourExtScr),dimension(:),allocatable::scr
#if defined TIGER_USE_OMP
    integer::numthreads,threadID
    numthreads = numberOfThreads
    ! also initialize some omp local communication if needed
    if(.not.allocated(omp_offsets_abcd_rec)) then
       allocate(omp_offsets_abcd_rec(num_internal+1:num_orbitals),stat=allocatestatus)
       call allocatecheck(allocatestatus,"omp_offsets_abcd_rec")
    endif
#else
    integer,parameter::numthreads = 1
    integer,parameter::threadID = 1
#endif


    !// INITIALIZE DERIVED TYPE ARRAYS
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%encountered = .false.
    lambda_path%num_singles = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (num_external < 4) then

       write(ioOutput,*) "**********************************************************"
       write(ioOutput,*) "***  THIS CODE DOES NOT HANDLE SYSTEMS WITH LESS       ***"
       write(ioOutput,*) "***  THAN FOUR VIRTUALS.  I SUGGEST YOU TRY PAPER      ***"
       write(ioOutput,*) "***  AND PENCIL FOR THIS CALCULATION.                  ***"
       write(ioOutput,*) "**********************************************************"
       call flush(ioOutput)
       stop

    endif
    
    icount = 0
    
    ! just allocate the required scratch space once. Also we kind of build our own thread-local heap space here
    ! in order to avoid both false cache sharing and avoid allocating within loops since ptmalloc2 seems to be build
    ! from chicken bones and does not satisfy such allocations from a thread-local cache...
    ! ideally, we would want to lazily allocate this by the thread (assuming pinned threads and different sockets) but
    ! this is low priority.
    allocate(scr(numthreads),stat=allocatestatus)
    call allocatecheck(allocatestatus,"all_external_3_4_loops_cho_2_scr")
    do a = 1, numthreads
      if(integralDirect .and. .not.fullyIntegralDirect .and. .not.sphere_based_integral_truncations .and. .not. directLowMem) then
       allocate(scr(a)%the_virtuals(num_orbitals),scr(a)%integrals_abcd(num_external,num_external*(num_external+1)/2),stat=allocatestatus)
      else
       allocate(scr(a)%the_virtuals(num_orbitals),scr(a)%integrals_abcd(num_external,num_external),stat=allocatestatus)
      endif
      ! is this too long? Yes. Does it avoid reallocations? Yes.
      scr(a)%resizeScr%length2 = num_external
      scr(a)%resizeScr%length2C2 = num_external*(num_external-1)/2
      allocate(scr(a)%resizeScr%integrals_abcd_mod2(scr(a)%resizeScr%length2C2), &
        scr(a)%resizeScr%integrals_abcd_mod3(scr(a)%resizeScr%length2C2), &
        scr(a)%resizeScr%integrals_abcd_mod4(scr(a)%resizeScr%length2),stat=allocatestatus)
      call allocatecheck(allocatestatus,"all_external_3_4_loops_cho_2")
      ! no need to allocate sigma here
      scr(a)%resizeScr%lengthSigma = 0
    enddo
        
   if(iteration == 1) then
    abcd_rec = 0
    
    !// FIRST THING WE NEED TO DO IS SEARCH TO THE Dab VERTEX.
    !// IF WE ARE OPERATING IN VECTORIZED MODE, WE SORT THE CI AND
    !// SIGMA VECTORS AND PASS A NEW COMMAND INTO THE DECIDER ROUTINE.

    ! first our dry run to get our data setup (should be cheap and is required to be serial)
    do a = num_internal+1, num_orbitals
#if defined TIGER_USE_OMP
       ! we create the offsets on the fly here
       omp_offsets_abcd_rec(a) = abcd_rec
#endif

       do c = num_internal+1, a 
       
          if (ignorable_pair(c,a) ) cycle
      
          call init_tree_search(graph, 0, 0, num_internal, num_elec-2)
          do while ( get_internal_path(lambda_path, graph))
             call build_ext_34seg_dryrun(lambda_path,scr(1)%the_virtuals,abcd_rec,a,c)
          end do
          
       enddo
    enddo
    
   endif ! just the first iteration

    abcd_rec = 0
    
    !$omp parallel &
    !$omp default(none) &
    !$omp private(icount,threadID,abcd_rec) &
    !$omp shared(civec,sigmavec,num_orbitals,num_internal,omp_offsets_abcd_icount, &
    !$omp omp_offsets_abcd_rec,scr,cho_data,integralDirect,sphere_based_integral_truncations, &
    !$omp fullyIntegralDirect)
    !$omp do &
    !$omp schedule(static)
    do a = num_internal+1, num_orbitals
#ifdef TIGER_USE_OMP
       if(.not. integralDirect .or. (sphere_based_integral_truncations .and. .not. fullyIntegralDirect)) then
         icount = omp_offsets_abcd_icount(a)
       endif
       threadID = OMP_get_thread_num()+1
       abcd_rec = omp_offsets_abcd_rec(a)
#endif

       call ext34_parallel_load(cho_data, civec, sigmavec, icount,scr(threadID)%integrals_abcd,scr(threadID)%the_virtuals, &
            scr(threadID)%resizeScr,a,abcd_rec,threadID)
       !call ext34_parallel_load(icount,a,abcd_rec,threadID)

    enddo
    !$omp end do nowait
    !$omp end parallel
    
    ! fortran standard specifies automatic deallocation of non-save allocatables
        
  end subroutine all_external_3_4_loops_cho_2
  !****************************************************************

  subroutine ext34_parallel_load(cho_data, civec, sigmavec, icount,integrals_abcd,the_virtuals,resizeScr,a,abcd_rec,threadID)
  implicit none
  integer,intent(in)::threadID,a
  integer,intent(inout)::abcd_rec
  integer,intent(inout)::icount
  real(real8),dimension(:,:),intent(inout)::integrals_abcd
  integer,dimension(:),intent(inout)::the_virtuals
  type(fourExtResizeScr),intent(inout):: resizeScr
  type(blockedLockVectorType),intent(inout):: civec, sigmavec
  type(cholesky_data),intent(in)::cho_data
  real(real8),dimension(:,:), allocatable::ab_mo_vecs,cd_mo_vecs
  
  integer::b,c,d,ab_ind,ab_1,stat,x,ind1,ind2,counter
  integer::nonredlength
  
   if(integralDirect .and. .not.sphere_based_integral_truncations .and. .not. fullyIntegralDirect) then
   
     nonredlength = num_external*(num_external+1)/2
     
     if(.not.directLowMem) then
        allocate(ab_mo_vecs(numcho,num_external),cd_mo_vecs(numcho,nonredlength),stat=stat)
     else
        allocate(ab_mo_vecs(numcho,num_external),cd_mo_vecs(numcho,num_external),stat=stat)
     endif
      
     call allocatecheck(stat,"ab_mo_vecs and cd_mo_vecs")
     
     ! read in this (ab| block, we will need it a couple of times
     do x = num_internal+1,num_orbitals
        ! the index of (ab|
        ind1 = max(a,x)
        ind1 = ind1*(ind1-1)/2+min(a,x)
        ind1 = cho_data%mo_ind_inv(ind1)
          
        ! read the stuff in
        if (ind1 .ne. 0 .and. .not. cdVecsInMemory) then
           call for_double_buf_readblock(mo_int_no, ind1, ab_mo_vecs(:,x-num_internal), threadID)
        else if(ind1 .ne. 0 .and. cdVecsInMemory) then
           ab_mo_vecs(:,x-num_internal) = cho_data%cho_vectors(:,ind1)
        else
           ab_mo_vecs(:,x-num_internal) = 0.0d0
        endif
     enddo
     
     
     if(.not.directLowMem) then
       ! read in all (non-redundant) |cd) in
       counter = 1
       do c = num_internal+1, num_orbitals
          do d = num_internal+1,c
                    
            ! the index of |cd)
            ind1 = c !max(c,x)
            ind1 = ind1*(ind1-1)/2+d!min(c,x)
            ind2 = cho_data%mo_ind_inv(ind1)
          
            ! read the stuff in
        
            if(counter > nonredlength) then
              write(*,*) "WTF!"
              flush(6)
              stop
            endif
          
            if (ind2 .ne. 0 .and. .not. cdVecsInMemory) then
               call for_double_buf_readblock(mo_int_no, ind2, cd_mo_vecs(:,counter), threadID)
            else if(ind2 .ne. 0 .and. cdVecsInMemory) then
               cd_mo_vecs(:,counter) = cho_data%cho_vectors(:,ind2)
            else
               cd_mo_vecs(:,counter) = 0.0d0
            endif
          
            counter = counter + 1
          
         enddo
       enddo
            
       if (counter /= nonredlength+1) then
         write(*,*) "WTF 2!!!!"
         flush(6)
         stop
       endif
          
       ! and here it comes, one single big dgemm
       call dgemm("T","N",num_external,nonredlength,numcho,1.0d0,ab_mo_vecs,numcho,cd_mo_vecs,numcho,0.0d0,integrals_abcd,num_external)
     endif
     
   endif

   loop1: do c = num_internal+1, a

     if(integralDirect .and. .not.sphere_based_integral_truncations .and. .not.fullyIntegralDirect .and. directLowMem) then
   
       ! the problem at hand is easier than one might think: for a given a,c we want all (ab|cd)
       ! this matrix is (albeit cut form a highly symmetric hypercube) not necessarily FULLY symmetric
       ! itself
       ! IMPORTANT: THIS DOES NOT KNOW ABOUT ANY TRUNCATIONS! IT WILL EVENTUALLY THEREFORE ONLY BE USED FOR
       !            CANONICAL/NON-LOCAL CALCULATIONS
       do x = num_internal+1,num_orbitals
                    
          ! the index of |cd)
          ind2 = max(c,x)
          ind2 = ind2*(ind2-1)/2+min(c,x)
          ind2 = cho_data%mo_ind_inv(ind2)
          
          ! read the stuff in
          if (ind2 .ne. 0 .and. .not. cdVecsInMemory) then
             call for_double_buf_readblock(mo_int_no, ind2, cd_mo_vecs(:,x-num_internal), threadID)
          else if(ind2 .ne. 0 .and. cdVecsInMemory) then
             cd_mo_vecs(:,x-num_internal) = cho_data%cho_vectors(:,ind2)
          else
             cd_mo_vecs(:,x-num_internal) = 0.0d0
          endif
          
       enddo
       
       ! just run one simple dgemm
       call dgemm("T","N",num_external,num_external,numcho,1.0d0,ab_mo_vecs,numcho,cd_mo_vecs,numcho,0.0d0,integrals_abcd,num_external)
         
     else if(integralDirect .and. fullyIntegralDirect) then
         
      ! slightly more involved: set everything to our flag, 42
      ! YES, I do realize this may complicate debugging in some weird ass case but since the flag is ONLY used to
      ! recompute the integral, if said integral should ACTUALLY be 42, we will recompute it and get 42 back out
      integrals_abcd = magic_number
      
     else if(integralDirect .and. .not. sphere_based_integral_truncations .and. .not.fullyIntegralDirect .and. .not. directLowMem) then
     
        ! done above w/ O(N^3) memory
        
     else
   
       if (ignorable_pair(c,a)) cycle
          
       integrals_abcd = 0.0D0
       do d = num_internal+1,num_orbitals
          if (ignorable_pair(d,a)) cycle
          if (ignorable_pair(d,c)) cycle 

          do b = num_internal+1,num_orbitals
             if (ignorable_pair(b,a) ) cycle
             if (ignorable_pair(b,c) ) cycle
             if (ignorable_pair(b,d) ) cycle 

             icount = icount + 1
             call for_double_buf_readElement(cd_abcd_no,icount,integrals_abcd(b-num_internal,d-num_internal),threadID)
          enddo
       enddo
     
     endif
     
     ab_ind = a*(a-1)/2+c
     ab_1 = ab_ind

     do while (ab_1 .eq. ab_ind)
        call ext_3_4_seg_loops_vec_lmo_res_2(a,c,civec, sigmavec, ab_1,integrals_abcd,abcd_rec,the_virtuals,cho_data,resizeScr)
     enddo

  enddo loop1
  
  end subroutine ext34_parallel_load

  !*****************************************************************
  !****************************************************************
  subroutine one_internal_seg_vec_1_cho(cho_data, civec, sigmavec, lambda_path,mod1vars)

    !// IN THIS ROUTINE WE CALL THE SEARCH ALGORITHM TO FIND
    !// THE 1 SEGMENT OPEN LOOPS IN THE INTERNAL SPACE.  AFTER
    !// THIS ROUTINE WE CALL A ROUTINE TO BUILD THE EXTERNAL SPACE
    !// CONTRIBUTIONS.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use two_seg_var_mod

    implicit none

    type(orbital_path)::lambda_path
    type(threefourmod1vars)::mod1vars
    type(blockedLockVectorType) :: civec, sigmavec
    type(cholesky_data),intent(in)::cho_data

    integer::i,a,b,c          !// THE INTERNAL LEVEL
    integer::i_1,icount
    integer::allocatestatus
    integer::iabc_rec,iabc_c
    type(threeExtScr),dimension(:),allocatable::scr
    
    integer::stat,idum,idum2,nonredlength
    real(real8),dimension(:,:),allocatable::bc_mo_vecs
        
#ifdef TIGER_USE_OMP
    integer::numthreads,ten_pointer
    numthreads = numberOfThreads
    
    ! also initialize some omp local communication if needed
    if(.not.allocated(omp_offsets_iabc_rec)) then
       allocate(omp_offsets_iabc_rec(num_internal),stat=allocatestatus)
       call allocatecheck(allocatestatus,"omp_offsets_iabc_rec")
    endif
#else
    integer,parameter::numthreads = 1
    integer,parameter::ten_pointer = 1
#endif

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// INITIALIZE DERIVED TYPE ARRAYS
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%encountered = .false.
    lambda_path%num_singles = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// REWIND THE INTEGRAL FILES FOR THE THREE SEGMENT LOOPS.
    !// THE INTEGRALS FOR THE FOUR SEGMENT LOOPS ARE HANDLED
    !// USING DIRECT ACCESS IN THE SUBSEQUENT ROUTINE.

    mod1vars%iabc_rec = 0
        
    if(mod1vars%iter_count .eq. 1) then
       
       icount = 0
       do i = 1, num_internal
        
#ifdef TIGER_USE_OMP
          ! we create the offsets on the fly here
          omp_offsets_iabc_rec(i) = mod1vars%iabc_rec
          ten_pointer = OMP_get_thread_num()+1
#endif

          mod1vars%i_ind = i
          
          !// REINITIALIZE FOR SUPER SAFETY!
          lambda_path%occupations = 0
          lambda_path%constraints = 0
          lambda_path%arc_weights = 0
          lambda_path%singles = 0
          lambda_path%num_singles = 0
          !// SET UP THE LIST OF CONSTRAINED LEVELS
          lambda_path%constraints(2,1) = i
          lambda_path%level1 = i
          lambda_path%constraints(2,2) = num_internal+1

          !// THREE SEGMENTS TO HANDLE HERE: {1},{5},{2}

!!!!!!!!!!!!!!!
          !//
          !//  {1} SEGMENT
          !//
!!!!!!!!!!!!!!!

          !// SET UP MU AND LAMBDA OCCUPATIONS
          lambda_path%constraints(1,1) = 0
          lambda_path%constraints(3,1) = 1

          !// INITIALIZE CONSTRAINT COUNT
          lambda_path%constraints(2,6) = 1

          !// NEXT TASK
          lambda_path%constraints(1,6) = CHO_IABC 

          !// STORE A LABEL TO REMIND US THAT THIS IS A THREE SEGMENT
          lambda_path%loop_type = 1

          !// DO THE SEARCH
          lambda_path%num_singles = 0
          call broken_constrained_search(lambda_path,0,0,mod1vars=mod1vars)


!!!!!!!!!!!!!!!
          !//
          !//  {5} SEGMENT
          !//
!!!!!!!!!!!!!!!

          !// REINITIALIZE FOR SUPER SAFETY!
          lambda_path%occupations = 0
          lambda_path%arc_weights = 0
          lambda_path%singles = 0
          lambda_path%num_singles = 0

          !// SET UP MU AND LAMBDA OCCUPATIONS
          lambda_path%constraints(1,1) = 1
          lambda_path%constraints(3,1) = 2

          !// INITIALIZE CONSTRAINT COUNT
          lambda_path%constraints(2,6) = 1

          !// NEXT TASK
          lambda_path%constraints(1,6) = CHO_IABC

          !// STORE A LABEL TO REMIND US THAT THIS IS A THREE SEGMENT
          lambda_path%loop_type = 2

          !// DO THE SEARCH
          lambda_path%num_singles = 0
          call broken_constrained_search(lambda_path,0,0,mod1vars=mod1vars)
          
       enddo
       
    endif ! second or later iteration
    
        allocate(scr(numthreads),stat=allocatestatus)
        call allocatecheck(allocatestatus,"one_internal_seg_vec_1_cho_scr")
    
    
        nonredlength = num_external*(num_external-1)+num_external
        do i = 1,numthreads
           allocate(scr(i)%Px(fsn(open_shells),fsn(open_shells)), &
           scr(i)%virt_lam_allow(num_external),scr(i)%virt_mu_allow(num_external),stat=allocatestatus)
           call allocatecheck(allocatestatus,"one_internal_seg_vec_1_cho")
           if(fullyIntegralDirect .and. .not. directSuperLowMem) then
              allocate(scr(i)%nonred_iabc_buffer(num_external,nonredlength),stat=stat)
              call allocatecheck(stat,"one_internal_seg_vec_1_cho: ia_mo_vecs")
           else if(integralDirect .and..not.sphere_based_integral_truncations .and. .not. fullyIntegralDirect .and. .not. directSuperLowMem) then
              allocate(scr(i)%ia_mo_vecs(numcho,num_external), &
                       scr(i)%nonred_iabc_buffer(num_external,nonredlength),stat=stat)
              call allocatecheck(stat,"one_internal_seg_vec_1_cho: ia_mo_vecs")
           else if(integralDirect .and. directSuperLowMem) then
              ! nothing, no memory
           else
              allocate(scr(i)%iabc_buffer(num_external,num_external,num_external), &
                       stat = stat)
              call allocatecheck(stat,"one_internal_seg_vec_1_cho: iabc_buffer")
           endif
       enddo
       
       if(integralDirect .and. .not.sphere_based_integral_truncations .and. .not. fullyIntegralDirect .and. .not. directSuperLowMem) then
       
           ! allocate and populate all non-redundant |bc)'s
           allocate(bc_mo_vecs(numcho,nonredlength),stat=stat)
           call allocatecheck(stat,"one_internal_seg_vec_1_cho: bc_mo_vecs")
           
           ! now read in all the non-redundant |bc) in
           do b = num_internal+1,num_orbitals
              do c = num_internal+1,b
             
                 idum = b*(b-1)/2+c
                 idum2 = cho_data%mo_ind_inv(idum)
                 if (idum2 .ne. 0) then
                    idum = (b-num_internal)*((b-num_internal)-1)/2+(c-num_internal)
                    if(cdVecsInMemory) then
                      bc_mo_vecs(:,idum) = cho_data%cho_vectors(:,idum2)
                    else
                      call for_double_buf_readblock(mo_int_no, idum2, bc_mo_vecs(:,idum), 1)
                    endif
                 else
                    bc_mo_vecs(:,idum) = 0.0d0
                 endif
                
              enddo
           enddo         
       endif

       
       icount = 0
       iabc_rec = 0
       
       !$omp parallel &
       !$omp default(none) &
       !$omp shared(civec,sigmavec,num_internal,omp_offsets_iabc_icount,omp_offsets_iabc_rec,ignorable_pair, &
       !$omp num_orbitals,scr,nonredlength,integralDirect,sphere_based_integral_truncations,numcho, &
       !$omp cho_data,num_external,bc_mo_vecs,fullyIntegralDirect,cdVecsInMemory,directSuperLowMem) &
       !$omp private(icount,ten_pointer,iabc_rec,iabc_c,i_1,idum,idum2)
       !$omp do &
       !$omp schedule(static)
       do i = 1, num_internal
       
#ifdef TIGER_USE_OMP
        ten_pointer = OMP_get_thread_num()+1
        iabc_rec = omp_offsets_iabc_rec(i)
#endif
       
        if(integralDirect .and. .not.sphere_based_integral_truncations .and. .not. fullyIntegralDirect .and. .not. directSuperLowMem) then
                 
          ! read in all (ia| for the i in question
          do a = num_internal+1,num_orbitals
             idum = a*(a-1)/2+i
             idum = cho_data%mo_ind_inv(idum)
             
             if (idum .ne. 0) then
                if(cdVecsInMemory) then
                  scr(ten_pointer)%ia_mo_vecs(:,a-num_internal) = cho_data%cho_vectors(:,idum)
                else
                  call for_double_buf_readblock(mo_int_no, idum, scr(ten_pointer)%ia_mo_vecs(:,a-num_internal), ten_pointer)
                endif
             else
                scr(ten_pointer)%ia_mo_vecs(:,a-num_internal) = 0.0d0
             endif
             
          enddo
          
          ! do one big ass matrix multiplication
          call dgemm("T","N",num_external,nonredlength,numcho,1.0d0,scr(ten_pointer)%ia_mo_vecs,numcho,bc_mo_vecs, &
               numcho,0.0d0,scr(ten_pointer)%nonred_iabc_buffer,num_external)
               
        else if(integralDirect .and. fullyIntegralDirect .and. .not. directSuperLowMem) then
        
          ! set it to our magic number flag... :-)
          scr(ten_pointer)%nonred_iabc_buffer = magic_number

        else if(integralDirect .and. directSuperLowMem) then
          ! nothing, we do not need memory :-)
        else
       
#ifdef TIGER_USE_OMP
          icount = omp_offsets_iabc_icount(i)
#endif

          scr(ten_pointer)%iabc_buffer = 0.0D0

          do c = num_internal+1, num_orbitals

             if (ignorable_pair(c,i) ) cycle

             iabc_c = c-num_internal  ! Careful, not c

             do b = num_internal+1,num_orbitals
                if (ignorable_pair(b,i) ) cycle
                if (ignorable_pair(b,c) ) cycle
                do a = num_internal+1,num_orbitals
                   if (ignorable_pair(a,c) ) cycle
                   if (ignorable_pair(a,b) ) cycle
                   if (ignorable_pair(a,i) ) cycle
                   icount = icount + 1
                   call for_double_buf_readElement(cd_iabc_no,icount,scr(ten_pointer)%iabc_buffer(a-num_internal,b-num_internal,iabc_c),ten_pointer)
                enddo
             enddo

          enddo
          
        endif

          i_1 = i
          do while (i_1 .eq. i)
             call iabc_sig_cho_big_res(civec, sigmavec, i_1,scr(ten_pointer)%virt_lam_allow,scr(ten_pointer)%virt_mu_allow, &
                                       iabc_rec,scr(ten_pointer)%Px, &
                                       scr(ten_pointer)%iabc_buffer, scr(ten_pointer)%nonred_iabc_buffer,cho_data,i)
          enddo
          
       enddo
       !$omp end do nowait
       !$omp end parallel
       
       ! fortran standard specifies automatic deallocation of non-save variables.

  end subroutine one_internal_seg_vec_1_cho
  !*****************************************************************
  !*****************************************************************

  !******************************************************************
  subroutine iabc_sig_cho_big(lambda_path,mod1vars)

    !// IN THIS SUBROUTINE WE TREAT THE CONTRIBUTION OF THREE
    !// AND FOUR SEGMENT LOOPS BY BUILDING THE EXTERNAL SPACE
    !// COMPLEMENT TO A SINGLE SEGMENT IN THE INTERNAL SPACE.
    !// THE CALCULATION IS PERFORMED IN VECTORIZED MODE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use locist_var_mod,only:locist_scratch

    implicit none

    type(orbital_path),intent(inout)::lambda_path
    type(threefourmod1vars),intent(inout)::mod1vars

    integer::start_elec                !// NUMBER OF ELECTRONS WHERE WE HAVE LEFT OFF
    integer::lam_or_mu                 !// SET TO 2 TO INDICATE WE ARE BUILDING A MU PATHS
    integer::constraint_count          !// KEEPS TRACK OF CONSTRAINED LEVELS
    integer::top_level                !// LABELS THE TOP AND BOTTOM LEVELS OF THE LOOP
    integer::mu_levels                 !// LABELS THE LEVELS IN THE MU PATH
    integer::path_elecs                !// NUMBER OF ELECTRONS IN THE MU PATH    
    integer::current_vertex            !// VERTICES IN THE MU PATHS
    integer::step_type                 !// STEPS IN THE MU PATH
    integer::loop_levels               !// LABELS LEVELS IN THE LOOP
    integer::head_levels               !// LABELS LEVELS IN THE HEAD OF THE LOOP
    integer::ibar                      !// POSITIONS OF IMPORTANT ORBITALS
    integer::i_segment                !// ORBITAL INDICES OF IMPORTANT ORBITALS
    integer::internal_singles          !// NUMBER OF SINGLY OCCUPIED ORBITALS IN INTERNAL PART OF
    !// CONFIGURATION
    integer::loop_type                 !// STORES THE LOOP TYPE                                   
    integer::lambda_weight,mu_weight   !// WEIGHTS AND ADDRESS OF LAMBDA AND MU PATHS
    integer::nm1_dim,nm2_dim

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    start_elec = sum(lambda_path%occupations(0:num_internal))
    if (start_elec /= num_elec-2) return

    !// STORE THE WEIGHT OF THE LAMBDA_PATH
    lambda_path%weight = sum(lambda_path%arc_weights(0:num_internal))
    if (skip_this_internal(lambda_path%weight,start_elec)) return

    !// INITIALIZE VARIABLES FOR BUILDING THE MU PATH
    loop_type = lambda_path%loop_type
    lam_or_mu = 2
    constraint_count = 1
    lambda_path%rt_loop_weight = 0
    top_level = lambda_path%level1
    i_segment = lambda_path%constraints(2,1)

    !// NUMBER OF ELECTRONS IN PATH BEFORE LOOP HEAD
    path_elecs = sum(lambda_path%occupations(0:top_level-1))

    !// NUMBER OF SINGLES IN THE INTERNAL SPACE PART OF PATH
    internal_singles = lambda_path%num_singles
    ibar = 0

    nm2_dim = 0
    nm1_dim = 0
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


!!!!!!!!!!!!!!!!!!!!
    !//
    !// GET IBAR
    !//
!!!!!!!!!!!!!!!!!!!!

    !// NOW GET IBAR
    do loop_levels = top_level, num_internal
       if (lambda_path%occupations(loop_levels)==1) ibar = ibar + 1
       if (loop_levels==i_segment) exit
    enddo

    !// ADJUST IBAR TO TAKE INTO ACCOUNT THE FACT THAT IT IS IN THE 
    !// MU PATH   
    if (loop_type==1) ibar = ibar + 1

    !// ADJUST  IBAR FOR THE HEAD SINGLES
    do head_levels = 1, top_level-1
       if (lambda_path%occupations(head_levels)==1) then
          ibar = ibar + 1
       endif
    enddo

!!!!!!!!!!!!!!!!!!!!!
    !//
    !//  DO MULTIPLICATIONS
    !//
!!!!!!!!!!!!!!!!!!!!!
    lambda_weight  = lambda_path%weight+1
    mu_weight      = lambda_path%rt_loop_weight+1

    if (loop_type == 1) then
        nm2_dim = fsn(internal_singles+2)
        nm1_dim = nm2_dim
    elseif (loop_type == 2) then
        nm2_dim = fsn(internal_singles+2)
        nm1_dim = fsn(internal_singles)
    endif

    dim1: if (nm2_dim /=0 .and. nm1_dim /=0) then

       if ((internal_index_vector2(lambda_weight) > 0 .or. internal_index_vector3(lambda_weight) > 0) .and. internal_index_vector1(mu_weight) > 0) then

          ! record for later
          mod1vars%iabc_rec = mod1vars%iabc_rec + 1
          call iabc_seg_writeData(mod1vars%iabc_rec,pseudo_iabc_info,1,mod1vars%i_ind,lambda_weight,mu_weight,internal_singles,ibar,loop_type)
          
       endif
    endif dim1

  end subroutine iabc_sig_cho_big
  
  !******************************************************************
  subroutine iabc_sig_cho_big_res(civec, sigmavec, i_1,virt_lam_allow,virt_mu_allow,iabc_rec,Px,iabc_buffer,iabc_symbuffer, &
    cho_data,i)

    !// IN THIS SUBROUTINE WE TREAT THE CONTRIBUTION OF THREE
    !// AND FOUR SEGMENT LOOPS BY BUILDING THE EXTERNAL SPACE
    !// COMPLEMENT TO A SINGLE SEGMENT IN THE INTERNAL SPACE.
    !// THE CALCULATION IS PERFORMED IN VECTORIZED MODE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer,dimension(:),intent(inout)::virt_lam_allow,virt_mu_allow
    integer,intent(inout)::iabc_rec
    real(real8),dimension(:,:),intent(inout)::Px
    real(real8),dimension(:,:,:),intent(in)::iabc_buffer
    real(real8),dimension(:,:),intent(inout)::iabc_symbuffer
    integer,intent(inout)::i_1
    integer,intent(in)::i
    type(blockedLockVectorType),intent(inout):: civec, sigmavec
    type(cholesky_data),intent(in)::cho_data

    integer::a,b                       !// INTEGERS TO LABEL EXTERNAL ORBITALS
    integer::c                         !// Label the N-1 hole 
    integer::ibar                      !// POSITIONS OF IMPORTANT ORBITALS
    integer::lambda_singles,mu_singles !// SINGLES IN LAMBDA AND MU PATHS
    integer::internal_singles          !// NUMBER OF SINGLY OCCUPIED ORBITALS IN INTERNAL PART OF
    !// CONFIGURATION
    integer::loop_type                 !// STORES THE LOOP TYPE                                   
    integer::lambda_weight,mu_weight   !// WEIGHTS AND ADDRESS OF LAMBDA AND MU PATHS
    integer::allocatestatus, deallocatestatus      !// FOR DYNAMIC MEMORY
    integer::nm2_dim,nm1_dim                       !// NM2 STANDS FOR N-2, NM1 = N-1
    integer::ab_spin,a_spin                        !// TO INDEX SPIN FUNCTIONS
    integer::mode_1                                !// SYMMETRIC OR ANTISYMMETRIC SCEPPER
    integer::a_address,ab_address,aa_address       !// A IS FOR N-1, AB AND AA FOR N-2, V FOR VALENCE
    integer::a_start,a_end
    integer::ab_start
    integer::aa_start

    integer::lam_length,lam_lengthC2
    integer::mu_length,mu_lengthC2
    integer::a1,b1,c1,i_2 
    integer::iost,x1,x2

    real(real8),parameter::zero = real(0.0,real8)
    real(real8)::cc,tmp

    real(real8),dimension(:,:,:),allocatable::iabc_small
    real(real8),dimension(:,:),allocatable::scep_ci_1,scep_sigma_1
    real(real8),dimension(:),allocatable::tmp_arr
    
    real(real8),external::ddot
#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num()+1
#else
    integer,parameter::threadID = 1
#endif

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    lambda_singles = 0
    mu_singles = 0
    a_address = 0
    ab_address = 0
    nm2_dim = 0
    aa_address = 0
    nm1_dim = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    i_2 = 0                   ! added so that in case of I/O failure i_2 is not uninitialized  
    iabc_rec = iabc_rec + 1
    call iabc_seg_readData(iabc_rec,pseudo_iabc_info,threadID,i_2,lambda_weight,mu_weight,internal_singles,ibar,loop_type,iost)
    if (i_1 .ne. i_2) then
       iabc_rec = iabc_rec-1
       i_1 = i_2
       return
    endif

    if (iost .ne. 0) then
       i_1 = 0
       return
    endif

    call iabc_info_simplifier (loop_type,ibar,& 
         internal_singles,lambda_singles, & 
         mu_singles,nm2_dim,nm1_dim,lambda_weight, & 
         mu_weight,a_address,ab_address,aa_address, & 
         Px) 


    dim1: if (nm2_dim /=0 .and. nm1_dim /=0) then 

       if ((ab_address > 0 .or. aa_address > 0) .and. a_address > 0) then 

          lam_length = num_allowed_virtuals(lambda_weight,"D") 
          lam_lengthC2 = lam_length*(lam_length-1)/2
          call get_virtuals(lambda_weight,"D", virt_lam_allow)

          mu_length = num_allowed_virtuals(mu_weight,"S")
          mu_lengthC2 = mu_length*(mu_length-1) 
          call get_virtuals(mu_weight,"S", virt_mu_allow)

          allocate(iabc_small(lam_length,lam_length,mu_length),tmp_arr(mu_length), &
                   scep_ci_1(lam_length,lam_length),scep_sigma_1(lam_length,lam_length),stat=allocatestatus)
          call allocatecheck(allocatestatus,"iabc_sig_cho_big_res")

          if(integralDirect .and. .not.sphere_based_integral_truncations .and. .not. fullyIntegralDirect .and. .not. directSuperLowMem) then
             do c=1,mu_length
                c1 = virt_mu_allow(c)
                do a=1,lam_length
                   a1 = virt_lam_allow(a)
                   do b=1,lam_length
                      b1 = virt_lam_allow(b)
                      x1 = max(a1,c1)-num_internal
                      x2 = min(a1,c1)-num_internal

                      iabc_small(b,a,c) = iabc_symbuffer(b1-num_internal,x1*(x1-1)/2+x2)
                   enddo
                enddo
             enddo
          
          else if(integralDirect .and. fullyIntegralDirect .and. .not. directSuperLowMem) then
          
            do c=1,mu_length
                c1 = virt_mu_allow(c)
                do a=1,lam_length
                   a1 = virt_lam_allow(a)
                   do b=1,lam_length
                      b1 = virt_lam_allow(b)
                      x1 = max(a1,c1)-num_internal
                      x2 = min(a1,c1)-num_internal

                      tmp = iabc_symbuffer(b1-num_internal,x1*(x1-1)/2+x2)
                      
                      if(tmp == magic_number) then
                        ! recalculate
                        tmp = makeOneIABC(i,b1,a1,c1,cho_data)
                        
                        ! store
                        iabc_symbuffer(b1-num_internal,x1*(x1-1)/2+x2) = tmp
                      endif
                      iabc_small(b,a,c) = tmp
                      
                   enddo
                enddo
             enddo

          else if(integralDirect .and. directSuperLowMem) then
             do c=1,mu_length
                c1 = virt_mu_allow(c)
                do a=1,lam_length
                   a1 = virt_lam_allow(a)
                   do b=1,lam_length
                      b1 = virt_lam_allow(b)
                      x1 = max(a1,c1)-num_internal
                      x2 = min(a1,c1)-num_internal

                      ! calculate
                      if(cdVecsInMemory) then
                        ! this also directly does CS prescreeing (therefore: if you keep the CD vectors in memory whilst running in low memory mode
                        ! and you DO NOT want any CS based truncation, you WILL need to specify an integral threshold of 0.0.
                        ! as all my tests indicate that running with rather high CS thresholds ( up to 1E-4) does not impact the total energy beyond
                        ! 1 kj/mol, I consider this to be an unimportant edge case (and actually, you WANT to run with CS prescreeing on)
                        tmp = makeOneIABC(i,b1,a1,c1,cho_data)
                      else
                        ! this would only do CS in the fullyIntegralDirect case (which actually never end up here)
                        tmp = makeOneIntegral(i,b1,a1,c1,cho_data,threadID)
                      endif

                      iabc_small(b,a,c) = tmp

                   enddo
                enddo
             enddo

          else
         
             do c=1,mu_length
                c1 = virt_mu_allow(c)
                do a=1,lam_length
                   a1 = virt_lam_allow(a)
                   do b=1,lam_length
                      b1 = virt_lam_allow(b)
                      iabc_small(b,a,c) = iabc_buffer(b1-num_internal,a1-num_internal,c1-num_internal)
                   enddo
                enddo
             enddo
          endif

          do ab_spin  = 1,nm2_dim

             if (internal_index_vector2(lambda_weight) < 0) cycle

             mode_1 = -1 
             if (ab_spin <= fsn(lambda_singles-2)) mode_1 = 1 
             if (loop_type == 2 .and. mode_1 == -1) mode_1 = -2 ! little... klotzy

             ab_start = ab_address + lam_lengthC2*(ab_spin-1)
             aa_start = aa_address + lam_length*(ab_spin-1)

             if (mode_1 .eq. 1) then
                 call scepper_diag_p1(scep_ci_1,civec%v,ab_start,&
                      lam_length,civec%v,aa_start)
             elseif (mode_1 .eq. -1) then
                 call scepper_diag_m1(scep_ci_1,civec%v,ab_start,&
                      lam_length)
             elseif (mode_1 .eq. -2) then
                 call scepper_diag_m2(scep_ci_1,civec%v,ab_start,&
                      lam_length)
             endif

             !// DOUBLET/N-2 INTERACTIONS ..
             ! jmd: this is kind of opaque, as we internally treat our 3D tensor as a 2D matrix
             ! also, this used to be computed within the following loop, which obviously is a bad idea
             ! I will keep the original code (well, the original didn't even use BLAS1) as reference therefore
             call dgemv('T',lam_length*lam_length,mu_length,1.d0,iabc_small,lam_length*lam_length,scep_ci_1,1,zero,tmp_arr,1)
             !do c = 1,mu_length
             !   tmp_arr(c) = cc*ddot(lam_length*lam_length,scep_ci_1,1,iabc_small(1,1,c),1)
             !enddo

             scep_sigma_1 = zero
             do a_spin = 1,nm1_dim 

                if (internal_index_vector1(mu_weight) < 0) cycle

                cc =  Px(ab_spin,a_spin)
                if (abs(cc) > integral_threshold) then 

                   a_start = a_address + mu_length*(a_spin-1)
                   a_end =   a_start   + mu_length- 1

                   !// N-2/DOUBLET INTERACTIONS ..
                   ! again, we formulate this in terms of a BLAS2 operation, treating the tensor as a matrix (thanks to the memory layout)
                   ! this call proudly presented by: WCW and JMD
                   ! find the original code below
                   call dgemv('N',lam_length*lam_length,mu_length,cc,iabc_small,lam_length*lam_length,civec%v(a_start),1,1.d0,scep_sigma_1,1)
                   !do a=1,mu_length
                   !   call daxpy(lam_length*lam_length,cc*civec%v(a_start+a-1),iabc_small(1,1,a),1,scep_sigma_1,1)
                   !enddo

#ifdef TIGER_USE_OMP
                   !call setSigmaLocks(a_start,a_end)
                   call setLocks(sigmavec%l,a_start,a_end) 
#endif
                   call daxpy(mu_length,cc,tmp_arr,1,sigmavec%v(a_start),1)
                   !sigmavec%v(a_start:a_end) = sigmavec%v(a_start:a_end)+ cc*tmp_arr
#ifdef TIGER_USE_OMP
                   !call unsetSigmaLocks(a_start,a_end)
                   call unsetLocks(sigmavec%l,a_start,a_end)
#endif

                endif

             enddo
             
             if (loop_type == 2 .and. mode_1 /=1 ) scep_sigma_1 = transpose(scep_sigma_1)
            
             if (mode_1 .eq. 1) then 
                call unscepper_diag_sigma_atomic_p1(scep_sigma_1,sigmavec,ab_start,lam_length,&
                        aa_start)
             elseif ((mode_1 .eq. -1) .or. (mode_1 .eq. -2)) then
                call unscepper_diag_sigma_atomic_m12(scep_sigma_1,sigmavec,ab_start,lam_length)
             endif
             
             
          enddo

          deallocate(iabc_small,tmp_arr,scep_ci_1,scep_sigma_1,stat=deallocatestatus)
          call deallocatecheck(deallocatestatus,"iabc_sig_cho_big_res")

       endif

    endif dim1


  end subroutine iabc_sig_cho_big_res
  !*******************************************************************
  !****************************************************************


  !*****************************************************************
  subroutine iabc_info_simplifier (loop_type,ibar,& 
       internal_singles,lambda_singles, & 
       mu_singles,lambda_dim,mu_dim,lambda_weight, & 
       mu_weight,a_address,ab_address,aa_address, & 
       Px) 

    !// THIS ROUTINE GETS THE COUPLING COEFFICIENTS AND PACKAGES
    !// THEM IN A NICE AND TRANSPARENT WAY FOR THE VECTORIZED
    !// MULTIPLICATION                                   

    !use spin_var_mod
    
    implicit none 


    integer,intent(in):: loop_type
    integer,intent(in):: ibar                                
    integer,intent(inout):: internal_singles,lambda_singles,mu_singles
    integer,intent(inout):: lambda_dim,mu_dim
    integer,intent(in):: lambda_weight,mu_weight
    integer,intent(out):: a_address,ab_address,aa_address

    real(real8),parameter:: sqrt2 = real(sqrt(real(2.0,real8)),real8)
    real(real8),dimension(:,:),intent(inout):: Px
    integer::maxD
    
    if (loop_type  == 1) then 

       !// CHECK ALL THIS STUFF . 

       lambda_singles = internal_singles + 2
       mu_singles     = lambda_singles

       lambda_dim = fsn(lambda_singles) 
       mu_dim     = fsn(mu_singles) 

       call contracted_cycle_stateless(ibar,lambda_singles-1,0,0,0,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)

       a_address = internal_index_vector1(mu_weight) 
       ab_address= internal_index_vector2(lambda_weight) 
       aa_address= internal_index_vector3(lambda_weight) 

    elseif (loop_type == 2) then 


       lambda_singles = internal_singles + 2 
       mu_singles     = internal_singles 

       lambda_dim = fsn(lambda_singles) 
       mu_dim     = fsn(mu_singles) 

       call contracted_cycle_stateless(lambda_singles-1,ibar,0,0,0,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)
       maxD = max(lambda_dim,mu_dim)
       Px(1:maxD,1:maxD) = Px(1:maxD,1:maxD) *sqrt2

       a_address = internal_index_vector1(mu_weight) 
       ab_address= internal_index_vector2(lambda_weight) 
       aa_address= internal_index_vector3(lambda_weight) 

    endif

  end subroutine iabc_info_simplifier

  !*****************************************************************
  !*****************************************************************
  subroutine one_internal_seg_vec_2_cho(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)

    !// IN THIS ROUTINE WE CALL THE SEARCH ALGORITHM TO FIND
    !// THE 1 SEGMENT OPEN LOOPS IN THE INTERNAL SPACE.  AFTER
    !// THIS ROUTINE WE CALL A ROUTINE TO BUILD THE EXTERNAL SPACE
    !// CONTRIBUTIONS.  

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use two_seg_var_mod

    implicit none

    type(cholesky_data),intent(in)::cho_data
    type(orbital_path),intent(inout)::lambda_path
    type(locist_scratch),intent(inout)::loc_scr
    type(threefourmod1vars),intent(inout)::mod1vars
    type(blockedLockVectorType),intent(inout):: civec, sigmavec

    integer::i          !// THE INTERNAL LEVEL
    integer::a,b
    integer::status
    integer::ia_count,counter

    real(real8),dimension(:,:),allocatable::ia_cho
    real(real8)::int1,int2

#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num()+1
#else
    integer,parameter::threadID = 1
#endif

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// INITIALIZE DERIVED TYPE ARRAYS
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%encountered = .false.

    lambda_path%num_singles = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(iaib)) then
       allocate(iaib(num_external*(num_external-1)/2),stat=status)
       call allocatecheck(status,"iaib    ")
    endif
    if (.not.allocated(iaia)) then
       allocate(iaia(num_external),stat=status)
       call allocatecheck(status,"iaia")
    endif
    if(.not.fullyIntegralDirect) then
      allocate(ia_cho(numcho,num_external),stat=status)
      call allocatecheck(status,"ia_cho")
      ia_cho = 0               ! THIS ARRAY MUST BE INITIALIZED TO ZERO OR YOU ARE GOING TO GET 
      ! THE WRONG ANSWER ... IM GOING TO NEED TO TAKE A LOOK AT THIS 
      ! LATER.  
      ! JMD: No idea who wrote this when. The current logic is pretty simple: if the local approximation
      !      dictates it, we skip reading in -> the vector should be zero, so we must initialize if we
      !      do not set this vector zero then. Easy as 3.14..... ;-)
    endif



    if (mod1vars%iter_count .gt. 1 .and. .not.(integralDirect .or. fullyIntegralDirect)) then
       rewind(unit=rediaia_no)
       rewind(unit=rediaib_no)
    endif

    if (mod1vars%iter_count .eq. 1) then

       do i = 1, num_internal


          ! read in all (ia| which are also obviously |ib)
          iaia = 0.0D0
          do a = num_internal+1,num_orbitals

             if (ignorable_pair(i,a) ) cycle

             ia_count = a*(a-1)/2+i
             ia_count = cho_data%mo_ind_inv(ia_count)

             if (ia_count .eq. 0) cycle
             
             if(.not. fullyIntegralDirect) then
               if(cdVecsInMemory) then
                 ia_cho(:,a-num_internal) = cho_data%cho_vectors(:,ia_count)
               else
                 call for_double_buf_readblock(mo_int_no, ia_count, ia_cho(:,a-num_internal), threadID)
               endif
               int1 = dot_product(ia_cho(:,a-num_internal),ia_cho(:,a-num_internal))
               iaia(a-num_internal) = int1
               if(.not.integralDirect) then
                 write(unit=rediaia_no) int1
               endif
             else
               iaia(a-num_internal) = makeOneIntegral(i,a,i,a,cho_data,threadID)
             endif
          enddo
          
          iaib = 0.0D0
          counter = 0
          do a = num_internal+1,num_orbitals-1
             do b = a+1,num_orbitals

                counter = counter + 1

                if (ignorable_pair(i,a) .or. ignorable_pair(i,b) ) cycle

                if(.not. fullyIntegralDirect) then
                  int2 = dot_product(ia_cho(:,a-num_internal),ia_cho(:,b-num_internal))
                  iaib(counter) = int2
                  if(.not.integralDirect) then
                    write(unit=rediaib_no) int2
                  endif
                else
                  iaib(counter) = makeOneIntegral(i,a,i,b,cho_data,threadID)
                endif
             enddo
          enddo

          !// REINITIALIZE FOR SUPER SAFETY!
          lambda_path%occupations = 0
          lambda_path%constraints = 0
          lambda_path%arc_weights = 0
          lambda_path%singles = 0
          lambda_path%num_singles = 0

!!!!!!!!!!!!!!!
          !// 
          !//  {2} SEGMENT
          !// 
!!!!!!!!!!!!!!!

          !// SET UP THE LIST OF CONSTRAINED LEVELS
          lambda_path%constraints(2,1) = i
          lambda_path%level1 = i
          lambda_path%constraints(2,2) = num_internal+1

          !// SET UP MU AND LAMBDA OCCUPATIONS
          lambda_path%constraints(1,1) = 0
          lambda_path%constraints(3,1) = 2

          !// INITIALIZE CONSTRAINT COUNT
          lambda_path%constraints(2,6) = 1

          !// NEXT TASK
!          if (virtual_truncation_flag .eq. 1 .and. local_ortho_mos .eq. 0) then
!             lambda_path%constraints(1,6) = 51
!          else if (virtual_truncation_flag .eq. 1 .and. local_ortho_mos .eq. 1) then
             lambda_path%constraints(1,6) = LMO_1_EXT_COMP
!          endif

          !// STORE A LABEL TO REMIND US THAT THIS IS A THREE SEGMENT
          lambda_path%loop_type = 5

          !// DO THE SEARCH
          lambda_path%num_singles = 0
          call broken_constrained_search(lambda_path,0,0,loc_scr=loc_scr,mod1vars=mod1vars,a=civec,b=sigmavec)
       enddo

    else ! second or later iteration, potential to read from storage

       do i = 1, num_internal


        if(integralDirect .and. .not. fullyIntegralDirect) then
          ! compute on the fly

          ! read in all (ia| which are also obviously |ib)
          iaia = 0.0D0
          do a = num_internal+1,num_orbitals

             if (ignorable_pair(i,a) ) cycle

             ia_count = a*(a-1)/2+i
             ia_count = cho_data%mo_ind_inv(ia_count)

             if (ia_count .eq. 0) cycle
             if(cdVecsInMemory) then
               ia_cho(:,a-num_internal) = cho_data%cho_vectors(:,ia_count)
             else
               call for_double_buf_readblock(mo_int_no, ia_count, ia_cho(:,a-num_internal), threadID)
             endif
             int1 = dot_product(ia_cho(:,a-num_internal),ia_cho(:,a-num_internal))
             iaia(a-num_internal) = int1
          enddo
          
          iaib = 0.0D0
          counter = 0
          do a=num_internal+1,num_orbitals-1
             do b = a+1,num_orbitals

                counter = counter + 1

                if (ignorable_pair(i,a) .or. ignorable_pair(i,b) ) cycle

                int2 = dot_product(ia_cho(:,a-num_internal),ia_cho(:,b-num_internal))
                iaib(counter) = int2
             enddo
          enddo
          
        else if(fullyIntegralDirect) then
        
        
          iaia = 0.0D0
          do a = num_internal+1,num_orbitals

             if (ignorable_pair(i,a) ) cycle

             iaia(a-num_internal) = makeOneIntegral(i,a,i,a,cho_data,threadID)
               
          enddo
          
          iaib = 0.0D0
          counter = 0
          do a=num_internal+1,num_orbitals-1
          
             do b = a+1,num_orbitals

                counter = counter + 1
                
                if (ignorable_pair(i,a) .or. ignorable_pair(i,b) ) cycle
                
                iaib(counter) = makeOneIntegral(i,a,i,b,cho_data,threadID)
             enddo
          enddo
          
        else
          ! read from storage
          iaia = 0.0D0
          iaib = 0.0D0
          counter = 0

          do a = num_internal+1,num_orbitals

             if (ignorable_pair(i,a) ) cycle

             read(unit=rediaia_no) iaia(a-num_internal)
          enddo

          do a = num_internal+1,num_orbitals-1
             do b = a+1,num_orbitals
                counter = counter + 1

                if (ignorable_pair(i,a) .or. ignorable_pair(i,b) ) cycle

                read(unit=rediaib_no) iaib(counter) 
             enddo
          enddo
          
        endif

          !// REINITIALIZE FOR SUPER SAFETY!
          lambda_path%occupations = 0
          lambda_path%constraints = 0
          lambda_path%arc_weights = 0
          lambda_path%singles = 0
          lambda_path%num_singles = 0

!!!!!!!!!!!!!!!
          !// 
          !//  {2} SEGMENT
          !// 
!!!!!!!!!!!!!!!

          !// SET UP THE LIST OF CONSTRAINED LEVELS
          lambda_path%constraints(2,1) = i
          lambda_path%level1 = i
          lambda_path%constraints(2,2) = num_internal+1

          !// SET UP MU AND LAMBDA OCCUPATIONS
          lambda_path%constraints(1,1) = 0
          lambda_path%constraints(3,1) = 2

          !// INITIALIZE CONSTRAINT COUNT
          lambda_path%constraints(2,6) = 1

          !// NEXT TASK
!          if (virtual_truncation_flag .eq. 1 .and. local_ortho_mos .eq. 0) then
!             lambda_path%constraints(1,6) = 51
!          else if (virtual_truncation_flag .eq. 1 .and. local_ortho_mos .eq. 1) then
             lambda_path%constraints(1,6) = LMO_1_EXT_COMP 
!          endif

          !// STORE A LABEL TO REMIND US THAT THIS IS A THREE SEGMENT
          lambda_path%loop_type = 5

          !// DO THE SEARCH
          lambda_path%num_singles = 0
          call broken_constrained_search(lambda_path,0,0,loc_scr=loc_scr,mod1vars=mod1vars,a=civec,b=sigmavec)
       enddo

    endif

    !// DEALLOCATE ARRAYS
    if(.not. fullyIntegralDirect) then
      deallocate(iaib,iaia,ia_cho,stat=status)
      call deallocatecheck(status,"one_internal_seg_vec_2_cho")
    endif


  end subroutine one_internal_seg_vec_2_cho
  !*****************************************************************
  !*****************************************************************
  subroutine all_i_3_seg_loopdrv(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)

    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use two_seg_var_mod
  
    implicit none

    type(cholesky_data)::cho_data
    type(orbital_path)::lambda_path
    type(locist_scratch)::loc_scr
    type(threefourmod1vars)::mod1vars
    type(graph_search_state)::graph
    type(blockedLockVectorType) :: civec, sigmavec

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// INITIALIZE DERIVED TYPES
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%encountered = .false.
    lambda_path%num_singles = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (mod1vars%iter_count .ne. 1) then
       rewind (unit=307)
    endif
    
    call init_tree_search(graph, 0, 0, num_internal, num_elec-2)
    do while ( get_internal_path(lambda_path, graph))
       call all_i_3_seg_loopdrv_csfs(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)
    end do

    !// INITIALIZE DERIVED TYPES
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%encountered = .false.
    lambda_path%num_singles = 0

    call init_tree_search(graph, 0, 0, num_internal, num_elec-1)
    do while ( get_internal_path(lambda_path, graph))
       call all_i_3_seg_loopdrv_csfs(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)
    end do

    !// INITIALIZE DERIVED TYPES
    lambda_path%occupations = 0
    lambda_path%singles = 0
    lambda_path%arc_weights = 0
    lambda_path%constraints = -1
    lambda_path%encountered = .false.
    lambda_path%num_singles = 0

    call init_tree_search(graph, 0, 0, num_internal, num_elec)
    do while ( get_internal_path(lambda_path, graph))
       call all_i_3_seg_loopdrv_csfs(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)
    end do

  end subroutine all_i_3_seg_loopdrv
  !*****************************************************************
  subroutine all_i_3_seg_loopdrv_csfs(civec, sigmavec, cho_data,lambda_path,loc_scr,mod1vars)

    use cholesky_structs
    use locist_var_mod,only:locist_scratch

    implicit none

    type(orbital_path)::lambda_path
    type(cholesky_data)::cho_data
    type(locist_scratch)::loc_scr
    type(threefourmod1vars)::mod1vars
    type(blockedLockVectorType) :: civec, sigmavec

    integer::loop_type
    integer::current_vertex         !// USED IN BUILDING MU PATH
    integer::step_type              !// KEEPS TRACK OF STEP IN MU PATH CONSTRUCTION
    integer::path_elecs             !// CUMULATIVE OCCUPATION IN MU PATH
    integer::top_level              !// LEVEL OF FIRST LOOP SEGMENT
    integer::mu_levels              !// LABELS LEVELS IN MU PATH
    integer::constraint_count       !// KEEP TRACK OF CONSTRAINTS IN BUILDING MU PATH
    integer::start_elec             !// OCCUPATION WHERE INTERNAL CSF MEETS UP WITH EXTERNAL SPACE
    integer::lam_step_1,lam_step_2,lam_step_3

    integer:: j,k,p,q,i_orb,j_orb                         !// LOOP VARIABLES
    integer::j_occ,k_occ   
    integer::level_1,level_2,level_3
    integer:: ij_count,ik_count,jk_count  

    real(real8):: ijik,ijjk,ikjk    !// INTEGRALS
    real(real8),dimension(:),allocatable:: ij_cho,jk_cho,ik_cho !// Cholesky vectors

    real(real8),parameter::zero=real(0.0,real8)
    
    real(real8), external::ddot ! BLAS1 (dot product)
    integer::status

#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num()+1
#else
    integer,parameter::threadID = 1
#endif
    
    !********************************************************


    lambda_path%weight = sum(lambda_path%arc_weights(0:num_internal)) 
    start_elec = sum(lambda_path%occupations(0:num_internal))
    if (skip_this_internal(lambda_path%weight,start_elec)) return

    if (mod1vars%iter_count .eq. 1) then
        allocate(ij_cho(numcho),jk_cho(numcho),ik_cho(numcho),stat=status)
        call allocatecheck(status,"all_i_3_seg_loopdrv_csfs")
    end if
        
    do p = 1, lambda_path%num_singles 
       i_orb = lambda_path%singles(p) 

       do j = 1,num_internal 
          if ( j .eq. i_orb ) cycle  
          j_occ = lambda_path%occupations(j)
          if (j_occ .eq. 1) cycle

          do k = 1,j-1

             if ( k .eq. i_orb  ) cycle  
             if ( ignorable_pair(i_orb,k)  .and. ignorable_pair(i_orb,j)   .and. ignorable_pair(j,k) ) cycle

             k_occ = lambda_path%occupations(k)
             if (k_occ .eq. 1) cycle
             if (j_occ .eq. k_occ) cycle

             level_1 = min(i_orb,k)
             level_3 = max(i_orb,j)
             level_2 = min(j,max(i_orb,k))

             !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             ! For Cholesky implementation, build the integrals here                                        ^
             !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

             if (mod1vars%iter_count .eq. 1) then

                ij_count = level_3*(level_3-1)/2+level_2
                ij_count = cho_data%mo_ind_inv(ij_count)

                ij_cho = 0.0D0
                if (ij_count .ne. 0) then
                   call for_double_buf_readblock(mo_int_no, ij_count, ij_cho(1:numcho), threadID)
                endif

                ik_count = level_3*(level_3-1)/2+level_1
                ik_count = cho_data%mo_ind_inv(ik_count) 

                ik_cho = 0.0D0 
                if (ik_count .ne. 0) then
                   call for_double_buf_readblock(mo_int_no, ik_count, ik_cho(1:numcho), threadID)
                endif

                ijik= ddot(numcho,ij_cho,1,ik_cho,1)

                jk_count = level_2*(level_2-1)/2+level_1
                jk_count = cho_data%mo_ind_inv(jk_count)

                jk_cho = 0.0D0
                if (jk_count .ne. 0) then
                   call for_double_buf_readblock(mo_int_no, jk_count, jk_cho(1:numcho), threadID)
                endif

                ijjk = ddot(numcho,ij_cho,1,jk_cho,1)

                ikjk = ddot(numcho,ik_cho,1,jk_cho,1)

                lambda_path%element1 = ijik
                lambda_path%element2 = ijjk
                lambda_path%element3 = ikjk


                write(unit=307) ijik,ijjk,ikjk

             else

                read(unit=307) ijik,ijjk,ikjk
                lambda_path%element1 = ijik
                lambda_path%element2 = ijjk
                lambda_path%element3 = ikjk
                !
             endif

             !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

             top_level = level_1
             lambda_path%constraints(2,1)=level_1
             lambda_path%constraints(2,2)=level_2
             lambda_path%constraints(2,3)=level_3
             lambda_path%constraints(2,4)=num_internal+1


             lam_step_1 = lambda_path%occupations(level_1)
             lam_step_2 = lambda_path%occupations(level_2)
             lam_step_3 = lambda_path%occupations(level_3)


             ! if *all* the integrals are below our threshold, we cycle as the kernel (all_i_seg_loopdrv_sigma)
             ! simply uses one of them as a multiplier
             if ( abs(ijik) .lt. integral_threshold .and. abs(ijjk) .lt. integral_threshold .and. &
                  abs(ikjk) .lt. integral_threshold) cycle

             loop_type_l: do loop_type = 1,6

                if (lam_step_1 .ne. loops(loop_type,1,1) ) cycle  
                if (lam_step_2 .ne. loops(loop_type,1,2) ) cycle 
                if (lam_step_3 .ne. loops(loop_type,1,3) ) cycle  

                lambda_path%rt_loop_weight = 0                 

                !// SET CONSTRAINTS FOR MU PATHS
                lambda_path%constraints(3,1:3) = &
                     loops(loop_type,2,1:3)

                constraint_count = 1
                path_elecs = sum(lambda_path%occupations(0:top_level-1))

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
                         cycle loop_type_l
                      endif

                   elseif(step_type == 1) then

                      if (add1(mu_levels-1,path_elecs)) then
                         path_elecs = path_elecs + 1
                         lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                              abs(y1(current_vertex))
                      else
                         cycle loop_type_l
                      endif

                   elseif(step_type == 0) then

                      if (.not.add0(mu_levels-1,path_elecs)) cycle loop_type_l

                   endif

                enddo

                !// FOR THE RT_LOOP_WEIGHT ADD IN THE WEIGHTS OF THE RELEVANT PARTS
                !// OF THE HEAD AND TAIL PATHS.  THIS WILL END UP BEING THE MU PATH WEIGHT
                lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                     sum(lambda_path%arc_weights(0:top_level-1))

                if (skip_this_internal(lambda_path%rt_loop_weight,path_elecs)) cycle loop_type_l


                !if (local_ortho_mos == 1) then
                call all_i_seg_loopdrv_sigma(civec, sigmavec, lambda_path,loc_scr,loop_type,mod1vars)
                !else 
                !   call all_i_seg_loopdrv_sigma_pao(lambda_path,loop_type)
                !endif

             enddo loop_type_l ! loop_type

          enddo !  k

       enddo !  j

    enddo !  num singles


    do p = 1, lambda_path%num_singles 
       i_orb = lambda_path%singles(p) 
       do q = p+1,lambda_path%num_singles
          j_orb = lambda_path%singles(q)  
          if ( ignorable_pair(i_orb,j_orb) ) cycle

          do k = 1,num_internal 
             if ( k .eq. i_orb ) cycle  
             if ( k .eq. j_orb ) cycle  
             if (lambda_path%occupations(k) .eq. 1) cycle

             if (ignorable_pair(i_orb,k) ) cycle
             if (ignorable_pair(j_orb,k) ) cycle

             level_1 = min(i_orb,k)
             level_3 = max(j_orb,k)
             level_2 = min(j_orb,max(i_orb,k))

             !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             ! For Cholesky implementation, build the integrals here                                        ^
             !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

             if (mod1vars%iter_count .eq. 1) then

                ij_count = level_3*(level_3-1)/2+level_2
                ij_count = cho_data%mo_ind_inv(ij_count)

                ij_cho = 0.0D0
                if (ij_count .ne. 0) then
                   call for_double_buf_readblock(mo_int_no, ij_count, ij_cho(1:numcho), threadID)
                endif

                ik_count = level_3*(level_3-1)/2+level_1
                ik_count = cho_data%mo_ind_inv(ik_count) 

                ik_cho = 0.0D0
                if (ik_count .ne. 0) then
                   call for_double_buf_readblock(mo_int_no, ik_count, ik_cho(1:numcho), threadID)
                endif

                ijik= ddot(numcho,ij_cho,1,ik_cho,1)

                jk_count = level_2*(level_2-1)/2+level_1
                jk_count = cho_data%mo_ind_inv(jk_count)

                jk_cho= 0.0D0
                if (jk_count .ne. 0) then
                   call for_double_buf_readblock(mo_int_no, jk_count, jk_cho(1:numcho), threadID)
                endif

                ijjk = ddot(numcho,ij_cho,1,jk_cho,1)

                ikjk = ddot(numcho,ik_cho,1,jk_cho,1)

                lambda_path%element1 = ijik
                lambda_path%element2 = ijjk
                lambda_path%element3 = ikjk

                write(unit=307) ijik,ijjk,ikjk

             else

                read(unit=307) ijik,ijjk,ikjk
                lambda_path%element1 = ijik
                lambda_path%element2 = ijjk
                lambda_path%element3 = ikjk
                !
             endif

             !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

             top_level = level_1
             lambda_path%constraints(2,1)=level_1
             lambda_path%constraints(2,2)=level_2
             lambda_path%constraints(2,3)=level_3
             lambda_path%constraints(2,4)=num_internal+1


             lam_step_1 = lambda_path%occupations(level_1)
             lam_step_2 = lambda_path%occupations(level_2)
             lam_step_3 = lambda_path%occupations(level_3)


             ! if *all* the integrals are below our threshold, we cycle as the kernel (all_i_seg_loopdrv_sigma)
             ! simply uses one of them as a multiplier
             if ( abs(ijik) .lt. integral_threshold .and. abs(ijjk) .lt. integral_threshold .and. &
                  abs(ikjk) .lt. integral_threshold) cycle

             loop_type_lX: do loop_type = 7,12

                if (lam_step_1 .ne. loops(loop_type,1,1) ) cycle  
                if (lam_step_2 .ne. loops(loop_type,1,2) ) cycle 
                if (lam_step_3 .ne. loops(loop_type,1,3) ) cycle  

                lambda_path%rt_loop_weight = 0                  

                !// SET CONSTRAINTS FOR MU PATHS
                lambda_path%constraints(3,1:3) = &
                     loops(loop_type,2,1:3)

                constraint_count = 1
                path_elecs = sum(lambda_path%occupations(0:top_level-1))

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
                         cycle loop_type_lX
                      endif

                   elseif(step_type == 1) then

                      if (add1(mu_levels-1,path_elecs)) then
                         path_elecs = path_elecs + 1
                         lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                              abs(y1(current_vertex))
                      else
                         cycle loop_type_lX
                      endif

                   elseif(step_type == 0) then

                      if (.not.add0(mu_levels-1,path_elecs)) cycle loop_type_lX

                   endif

                enddo

                !// FOR THE RT_LOOP_WEIGHT ADD IN THE WEIGHTS OF THE RELEVANT PARTS
                !// OF THE HEAD AND TAIL PATHS.  THIS WILL END UP BEING THE MU PATH WEIGHT
                lambda_path%rt_loop_weight = lambda_path%rt_loop_weight +&
                     sum(lambda_path%arc_weights(0:top_level-1))

                if (skip_this_internal(lambda_path%rt_loop_weight,path_elecs)) cycle loop_type_lX

                !if ( local_ortho_mos == 1) then
                call all_i_seg_loopdrv_sigma(civec, sigmavec, lambda_path,loc_scr,loop_type,mod1vars)
                !else
                !   call all_i_seg_loopdrv_sigma_pao(lambda_path,loop_type)
                !endif

             enddo loop_type_lX ! loop_type

          enddo !  k

       enddo !  num_singles

    enddo !  num singles
    
    if (mod1vars%iter_count .eq. 1) then
       deallocate(ij_cho,jk_cho,ik_cho,stat=status)
       call deallocatecheck(status,"all_i_3_seg_loopdrv_csfs")
    end if

  end subroutine all_i_3_seg_loopdrv_csfs
  !*****************************************************************
  subroutine all_i_seg_loopdrv_sigma(civec, sigmavec, lambda_path,loc_scr,loop_type,mod1vars)

    !// THIS SUBROUTINE BUILDS UP THE COMPLEMENTS TO THE THREE SEGMENT LOOPS
    !// WHICH RESIDE ENTIRELY IN THE INTERNAL SPACE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use locist_var_mod,only:locist_scratch

    implicit none

    type(orbital_path)::lambda_path
    type(locist_scratch)::loc_scr
    type(threefourmod1vars)::mod1vars
    type(blockedLockVectorType) :: civec, sigmavec

    integer::ibar,jbar              !// ORBITAL POSITION INDICES
    integer::internal_singles       !// NUMBER OF SINGLES IN INTERNAL PATH OF PATH
    integer::lambda_singles         !// SINGLES IN LAMBDA AND MU PATH
    integer::mu_singles
    integer::lambda_dim,mu_dim      !// DIMENSIONS FOR LAMBDA AND MU CONFIGURATIONS
    integer::lambda_weight          !// WEIGHTS OF LAMBDA AND MU PATH
    integer::mu_weight             
    integer::lambda_address         !// ADDRESS OF MU AND LAMBDA PATHS
    integer::mu_address
    integer::start_elec             !// OCCUPATION WHERE INTERNAL CSF MEETS UP WITH EXTERNAL SPACE
    integer::s1,s2,s3,s4            !// NUMBERS OF SINGLES IN LOOP SEGMENTS
    integer::loop_type              !// LABELS THE LOOP TYPE

    real(real8)::the_integral       !// THE INTEGRAL WE WILL NEED

    integer::lambda_start,lambda_start2,mu_start,mu_start2
    integer::lambda_end,lambda_end2,mu_end,mu_end2
    integer::aa_lam_address,aa_mu_address      

    integer::mu_length1
    integer::mu_length2,mu_length2C2
    integer::lam_length1
    integer::lam_length2,lam_length2C2

    integer::mode_1,mode_2

    integer::i,j 
    real(real8)::cc
    integer::status


    !//////////// EXTRA VARIABLES

    integer::common_virt_length,lam_offset,mu_offset,i_dum,j_dum,i_x,j_x,i_y,j_y
    !>todo REMOVE ULTIMATELY!
    integer,dimension(:),allocatable::common_virt,virt_lam_long,virt_lam_pos,virt_mu_long,virt_mu_pos,&
         common_virt_lam_pos,common_virt_mu_pos
    real(real8),dimension(:),allocatable::sig_nm2_lam,sig_nm2_mu
    real(real8),dimension(:),allocatable::sig_nm2_lam_d,sig_nm2_mu_d
    
    allocate(common_virt(num_orbitals),virt_lam_long(num_orbitals),virt_lam_pos(num_orbitals), &
         virt_mu_long(num_orbitals),virt_mu_pos(num_orbitals),common_virt_lam_pos(num_orbitals), &
         common_virt_mu_pos(num_orbitals),sig_nm2_lam(num_external*(num_external-1)/2), &
         sig_nm2_mu(num_external*(num_external-1)/2),sig_nm2_lam_d(num_external), &
         sig_nm2_mu_d(num_external),stat=status)
    call allocatecheck(status,"all_i_seg_loopdrv_sigma I")
    
    
    ! START AUTOGENERATED INITIALIZATION 
    mu_length2c2 = 0
    s1 = 0
    s2 = 0
    s3 = 0
    common_virt_mu_pos = 0
    lambda_start = 0
    mu_dim = 0
    i = 0
    i_x = 0
    i_y = 0
    mu_address = 0
    aa_lam_address = 0
    mu_length2 = 0
    mu_length1 = 0
    mu_weight = 0
    aa_mu_address = 0
    mu_offset = 0
    j_y = 0
    lambda_singles = 0
    j_x = 0
    the_integral = 0.0
    mu_singles = 0
    jbar = 0
    lambda_weight = 0
    mode_1 = 0
    mode_2 = 0
    lam_length1 = 0
    lam_length2 = 0
    j_dum = 0
    lambda_address = 0
    ibar = 0
    mu_start2 = 0
    s4 = 0
    mu_end = 0
    mu_start = 0
    cc = 0.0
    lam_length2c2 = 0
    mu_end2 = 0
    virt_mu_pos = 0
    virt_lam_pos = 0
    lambda_start2 = 0
    lambda_dim = 0
    common_virt_lam_pos = 0
    j = 0
    i_dum = 0
    lambda_end = 0
    lambda_end2 = 0
    virt_mu_long = 0
    common_virt_length = 0
    lam_offset = 0
    ! END AUTOGENERATED INITIALIZATION 

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    start_elec = sum(lambda_path%occupations(0:num_internal))
    internal_singles = lambda_path%num_singles

    sig_nm2_mu_D = 0.0
    sig_nm2_lam_d = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !//
    !// GET NUMBERS OF SINGLES
    !//
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call get_numbers_of_singles(s1,s2,s3,s4,lambda_path)                            

    call three_seg_indices_info(loop_type,ibar,jbar,s1,s2,s3,lambda_path,&  
         the_integral)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// THREE CASES HERE: ONE FOR THE D,S,AND V VERTICES
    if (start_elec == num_elec-2) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !// 
       !//  TWO VIRTUALS SINGLY OCCUPIED AND ONE VIRTUAL DOUBLY OCCUPIED
       !// 
!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !// GET DIMENSIONS
       if (loop_type <= 6) then
          lambda_singles = internal_singles+2
          mu_singles = lambda_singles
       else
          lambda_singles = internal_singles+2
          mu_singles = lambda_singles-2
       endif

       lambda_dim = fsn(lambda_singles)
       mu_dim     = fsn(mu_singles)

       lambda_address = internal_index_vector2(lambda_path%weight+1) 
       mu_address     = internal_index_vector2(lambda_path%rt_loop_weight +1) 

       aa_lam_address = internal_index_vector3(lambda_path%weight+1) 
       aa_mu_address  = internal_index_vector3(lambda_path%rt_loop_weight +1) 




       lambda_weight  = lambda_path%weight+1
       mu_weight      = lambda_path%rt_loop_weight+1

       !DEBUG expanded if statement !
       if (lambda_address > 0 .and. mu_address > 0) then  

          if (lambda_dim /= 0 .and. mu_dim /= 0) then

             !// GET COUPLING COEFFICIENTS
             call compute_3seg_coupl_coeff_vec(lambda_singles,mu_singles,&
                  loop_type,ibar,jbar,mod1vars%Px)

             lam_length2 = num_allowed_virtuals(lambda_weight,"D") 
             lam_length2C2 = lam_length2*(lam_length2-1)/2
             mu_length2  = num_allowed_virtuals(mu_weight,"D") 
             mu_length2C2 = mu_length2*(mu_length2-1)/2

             call get_virtuals(lambda_weight,"D", loc_scr%virt_lam_allow)
             call get_virtuals(mu_weight,"D", loc_scr%virt_mu_allow)


             common_virt = 0
             virt_lam_long = 0
             virt_lam_pos = 0
             virt_mu_long = 0
             virt_mu_pos = 0

             do i = 1, lam_length2
                i_dum = loc_scr%virt_lam_allow(i)
                virt_lam_long(i_dum) = i_dum
                virt_lam_pos(i_dum) = i
             enddo

             do i = 1, mu_length2
                i_dum = loc_scr%virt_mu_allow(i)
                virt_mu_long(i_dum) = i_dum
                virt_mu_pos(i_dum) = i
             enddo


             common_virt_length = 0
             common_virt_lam_pos = 0
             common_virt_mu_pos = 0
             do i = 1, lam_length2
                i_dum = loc_scr%virt_lam_allow(i)
                if ( virt_lam_long(i_dum) == virt_mu_long(i_dum) ) then
                   common_virt_length = common_virt_length + 1
                   common_virt(common_virt_length) = i_dum 
                   common_virt_lam_pos(common_virt_length) = virt_lam_pos(i_dum)
                   common_virt_mu_pos(common_virt_length) = virt_mu_pos(i_dum)
                endif
             enddo


             do i=1,lambda_dim

                if (internal_index_vector2(lambda_path%weight +1) < 0) cycle

                mode_1 = -1 
                if (i <= fsn(lambda_singles-2)) mode_1 = 1 

                lambda_start = lambda_address + lam_length2C2*(i-1)
                lambda_end   = lambda_start + lam_length2C2 - 1

                lambda_start2 = aa_lam_address + lam_length2*(i-1)
                lambda_end2   = lambda_start2 + lam_length2 - 1


                do j=1,mu_dim


                   if (internal_index_vector2(lambda_path%rt_loop_weight+1) < 0) cycle

                   mode_2 = -1 
                   if (j <= fsn(mu_singles-2)) mode_2 = 1 


                   mu_start = mu_address + mu_length2C2*(j-1)
                   mu_end   = mu_start + mu_length2C2- 1

                   mu_start2 = aa_mu_address +  mu_length2*(j-1)
                   mu_end2   = mu_start2 +  mu_length2 - 1


                   cc = mod1vars%Px(i,j)
                   if (abs(cc) > integral_threshold) then 

                      !////////////////////////////////////////////////////////////////////////////////!
                      !// NOW FOR THE SIGMA_{LAMBDA} !SEE THE DIMENSION OF SIGMA_NM2_LAMBDA


                      sig_nm2_lam = 0.0D0
                      sig_nm2_lam_d = 0.0D0

                      if (mode_2 == 1) then ! DEBUG: Changed mode_1 into mode_2 because aa_mu_start is the critical vairable here. CMK

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2
                            mu_offset = (i_y-1)*(i_y-2)/2 + mu_start - 1 ! Include the (mu_start-1) here
                            do j_dum = 1,i_dum-1
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_lam(lam_offset+j_x) = cc*the_integral*civec%v(mu_offset+j_y)
                            enddo
                            sig_nm2_lam_d(i_x) = cc*the_integral*civec%v(mu_start2+i_y-1)
                         enddo

                      else

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2
                            mu_offset = (i_y-1)*(i_y-2)/2 + mu_start - 1 ! Include the (mu_start-1) here
                            do j_dum = 1,i_dum-1
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_lam(lam_offset+j_x) = cc*the_integral*civec%v(mu_offset+j_y)
                            enddo
                         enddo

                      endif

                      !/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                      sig_nm2_mu= 0.0D0
                      sig_nm2_mu_d = 0.0D0

                      if (mode_1 == 1) then ! DEBUG: Changed mode_2 into mode_1 because aa_lambda_start is the critical vairable here. CMK

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2 + lambda_start-1 ! Include lambda_start-1 here
                            mu_offset = (i_y-1)*(i_y-2)/2 
                            do j_dum = 1,i_dum-1
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_mu(mu_offset+j_y) = cc*the_integral*civec%v(lam_offset+j_x)
                            enddo
                            sig_nm2_mu_d(i_y) = cc*the_integral*civec%v(lambda_start2+i_x-1)
                         enddo

                      else

                         do i_dum = 1, common_virt_length
                            i_x = common_virt_lam_pos(i_dum)
                            i_y = common_virt_mu_pos(i_dum)
                            lam_offset = (i_x-1)*(i_x-2)/2 + lambda_start-1 ! Include lambda_start-1 here
                            mu_offset = (i_y-1)*(i_y-2)/2
                            do j_dum = 1,i_dum-1
                               j_x = common_virt_lam_pos(j_dum)
                               j_y = common_virt_mu_pos(j_dum)
                               sig_nm2_mu(mu_offset+j_y) = cc*the_integral*civec%v(lam_offset+j_x)
                            enddo
                         enddo

                      endif
                      !////////////////////////////////////////////////////////////////////////////////////////////////////////

                      sigmavec%v(lambda_start:lambda_end) = sigmavec%v(lambda_start:lambda_end) + sig_nm2_lam(1:lam_length2C2)
                      if (aa_lam_address > 0) then
                            sigmavec%v(lambda_start2:lambda_end2) = sigmavec%v(lambda_start2:lambda_end2) &
                                   + sig_nm2_lam_d(1:lam_length2)
                      endif

                      sigmavec%v(mu_start:mu_end) = sigmavec%v(mu_start:mu_end) + sig_nm2_mu(1:mu_length2C2)
                      if (aa_mu_address > 0) then
                              sigmavec%v(mu_start2:mu_end2) = sigmavec%v(mu_start2:mu_end2) + sig_nm2_mu_d(1:mu_length2)
                      endif

                   endif

                enddo

             enddo

          endif

       endif

    elseif (start_elec == num_elec -1) then 

!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !// 
       !//  ONE VIRTUAL SINGLY OCCUPIED
       !// 
!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !// GET DIMENSIONS
       if (loop_type <= 6) then
          lambda_singles = internal_singles+1
          mu_singles = lambda_singles
       else
          lambda_singles = internal_singles+1
          mu_singles = lambda_singles-2
       endif

       lambda_dim = fsn(lambda_singles)
       mu_dim = fsn(mu_singles)

       lambda_address = internal_index_vector1(lambda_path%weight+1) 
       mu_address     = internal_index_vector1(lambda_path%rt_loop_weight +1) 

       lambda_weight  = lambda_path%weight+1
       mu_weight      = lambda_path%rt_loop_weight+1

       if (lambda_address > 0 .and. mu_address > 0) then  

          if (lambda_dim /= 0 .and. mu_dim /= 0) then

             !// GET COUPLING COEFFICIENTS
             call compute_3seg_coupl_coeff_vec(lambda_singles,mu_singles,&
                  loop_type,ibar,jbar,mod1vars%Px)

             lam_length1 = num_allowed_virtuals(lambda_weight,"S") 
             mu_length1 = num_allowed_virtuals(mu_weight,"S") 
             call get_virtuals(lambda_weight,"S", loc_scr%virt_lam_allow)
             call get_virtuals(mu_weight,"S", loc_scr%virt_mu_allow)

             ! Making Changes for the lmo version

             common_virt = 0
             virt_lam_long = 0
             virt_lam_pos = 0
             virt_mu_long = 0
             virt_mu_pos = 0

             do i = 1, lam_length1
                i_dum = loc_scr%virt_lam_allow(i)
                virt_lam_long(i_dum) = i_dum 
                virt_lam_pos(i_dum) = i
             enddo

             do i = 1, mu_length1
                i_dum = loc_scr%virt_mu_allow(i)
                virt_mu_long(i_dum) = i_dum
                virt_mu_pos(i_dum) = i
             enddo

             common_virt_length = 0
             common_virt_lam_pos = 0
             common_virt_mu_pos = 0

             do i = 1, lam_length1
                i_dum = loc_scr%virt_lam_allow(i)
                if ( virt_lam_long(i_dum) == virt_mu_long(i_dum) ) then
                   common_virt_length = common_virt_length + 1
                   common_virt(common_virt_length) = i_dum
                   common_virt_lam_pos(common_virt_length) = virt_lam_pos(i_dum)
                   common_virt_mu_pos(common_virt_length) = virt_mu_pos(i_dum)
                endif
             enddo

             do i=1,lambda_dim

                if (internal_index_vector1(lambda_path%weight +1) < 0) cycle

                lambda_start = lambda_address + lam_length1*(i-1)
                lambda_end   = lambda_start + lam_length1- 1

                do j=1,mu_dim

                   if (internal_index_vector1(lambda_path%rt_loop_weight +1) < 0) cycle

                   mu_start = mu_address + mu_length1*(j-1)
                   mu_end   = mu_start + mu_length1- 1

                   cc = mod1vars%Px(i,j)
                   if (abs(cc) > integral_threshold) then 

                      !////////////////////////////////////////////////////////////////////////////////!

                      !// EVALUATE SIGMA_LAMBDA IN PAO BASIS !SEE THE DIMENSION OF SIGMA_NM2_LAMBDA

                      do i_dum = 1, common_virt_length
                         i_x = common_virt_lam_pos(i_dum)
                         i_y = common_virt_mu_pos(i_dum)
                         sigmavec%v(lambda_start+i_x-1) = sigmavec%v(lambda_start+i_x-1)+ &
                              cc*the_integral*civec%v(mu_start+i_y-1)
                      enddo

                      !/////////////////////////////////////////////////////////////////////////////////////////////

                      do i_dum = 1, common_virt_length
                         i_x = common_virt_lam_pos(i_dum)
                         i_y = common_virt_mu_pos(i_dum)
                         sigmavec%v(mu_start+i_y-1) = sigmavec%v(mu_start+i_y-1) + &
                              cc*the_integral*civec%v(lambda_start+i_x-1)
                      enddo
                      !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                   endif

                enddo
             enddo

          endif

       endif

    elseif (start_elec == num_elec) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !// 
       !//  VALENCE STATE; NO EXTERNALS OCCUPIED
       !// 
!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !// GET DIMENSIONS
       if (loop_type <= 6) then
          lambda_singles = internal_singles
          mu_singles = lambda_singles
       else
          lambda_singles = internal_singles
          mu_singles = lambda_singles-2
       endif
       lambda_dim = fsn(lambda_singles)
       mu_dim = fsn(mu_singles)

       lambda_address = internal_index_vector0(lambda_path%weight+1) 
       mu_address     = internal_index_vector0(lambda_path%rt_loop_weight +1) 

       if (lambda_address > 0 .and. mu_address > 0) then  

          if (lambda_dim /= 0 .and. mu_dim /= 0) then

             !// GET COUPLING COEFFICIENTS
             call compute_3seg_coupl_coeff_vec(lambda_singles,mu_singles,&
                  loop_type,ibar,jbar,mod1vars%Px)

             do i=1,lambda_dim

                lambda_start = lambda_address + i-1 
                lambda_end   = lambda_start 

                do j=1,mu_dim

                   mu_start = mu_address + j-1
                   mu_end   = mu_start 

                   cc = mod1vars%Px(i,j)
                   if (abs(cc) > integral_threshold) then 

                      sigmavec%v(lambda_start:lambda_end) = sigmavec%v(lambda_start:lambda_end)+ & 
                           cc*the_integral*civec%v(mu_start:mu_end)

                      sigmavec%v(mu_start:mu_end)         = sigmavec%v(mu_start:mu_end) + & 
                           cc*the_integral*civec%v(lambda_start:lambda_end)

                   endif

                enddo
             enddo

          endif

       endif

    endif
    
    deallocate(common_virt,virt_lam_long,virt_lam_pos, &
         virt_mu_long,virt_mu_pos,common_virt_lam_pos, &
         common_virt_mu_pos,sig_nm2_lam, &
         sig_nm2_mu,sig_nm2_lam_d, &
         sig_nm2_mu_d,stat=status)
    call deallocatecheck(status,"all_i_seg_loopdrv_sigma I")

  end subroutine all_i_seg_loopdrv_sigma
  !*****************************************************************
  !*****************************************************************
  !*****************************************************************
  subroutine three_seg_indices_info(loop_type,ibar,jbar,s1,s2,s3,lambda_path,& 
       the_integral)

    !// THIS ROUTINE GETS THE COUPLING COEFFICIENT INFORMATION
    !// FOR THE 3 SEGMENT LOOPS WHICH RESIDE ENTIRELY IN 
    !// THE INTERNAL SPACE AND PACKAGES THEM NICELY
    !// FOR THE VECTORIZED MULTIPLICATION                                   

    implicit none
    type (orbital_path),intent(in)::lambda_path 

    integer,intent(in):: loop_type 
    integer,intent(out):: ibar,jbar
    integer,intent(in):: s1,s2,s3 

    real(real8),intent(out)::the_integral       !// THE INTEGRAL WE WILL NEED

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !//
    !// GET IBAR, JBAR, THE INTEGRAL
    !//
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (loop_type <= 2.or.&
         loop_type == 7.or.&
         loop_type == 8) then

       jbar = s1 
       ibar = s1 + s2
       the_integral = -lambda_path%element1

    elseif (loop_type <= 4.or.&
         loop_type == 9.or.&
         loop_type == 10) then

       jbar = s1
       ibar = s1 + s2 + s3
       the_integral = -lambda_path%element2

    else

       jbar = s1 + s2
       ibar = s1 + s2 + s3
       the_integral = -lambda_path%element3

    endif

    if (loop_type <= 6) jbar = jbar + 1

    if (loop_type == 7.or.&
         loop_type == 9.or.&
         loop_type == 11) the_integral = - the_integral


  end subroutine three_seg_indices_info

  !*****************************************************************
  subroutine compute_3seg_coupl_coeff_vec(lambda_singles,mu_singles,&
       loop_type, ibar, jbar,Px)

    !// THIS ROUTINE IS A GENERIC ROUTINE THAT COMPUTES COUPLING COEFFICIENTS
    !// FOR THE THREE SEGMENT CASE.  THE COUPLING COEFFICIENTS ARE STORED IN CONTRACTED

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer::lambda_singles        !// NUMBERS OF SINGLES
    integer::mu_singles            
    integer::loop_type             !// LABELS THE LOOPS
    integer::ibar                  !// NEED THESE INDICES FOR COUPLING COEFFICIENTS
    integer::jbar

    real(real8),dimension(:,:)::Px 
    real(real8),parameter::sqrt2=real(sqrt(real(2.0,real8)),real8)
    integer::maxD

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    maxD = max(fsn(lambda_singles),fsn(mu_singles))
    if (ibar == jbar) then

       !// JUST NEED AN IDENTITY MATRIX
       !contracted = identity_matrix(fsn(open_shells))
       Px(1:maxD,1:maxD) = identity_matrix(maxD)

       !if (loop_type > 6) contracted = sqrt2*contracted
       if (loop_type > 6) Px(1:maxD,1:maxD) = sqrt2*Px(1:maxD,1:maxD)

    elseif (loop_type <= 6) then
    
       !// HERE WE NEED CYCLES
       !contracted = cycles(index2m1(ibar,jbar),1:fsn(open_shells):1,&
       Px(1:maxD,1:maxD) = transpose(cycles(1:maxD,1:maxD,index2m1(ibar,jbar)))

    else

       !// HERE WE NEED CONTRACTED CYCLES
       !call contracted_cycle(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles)
       call contracted_cycle_stateless(jbar,ibar,0,0,jbar,0,lambda_singles,mu_singles,Px,cycles,sft,spin_matrix)

    endif

    !Px = contracted 

  end subroutine compute_3seg_coupl_coeff_vec

  !*****************************************************************
  !****************************************************************
  !****************************************************************
  Subroutine extern_one_seg_compl_vec_lmo(civec, sigmavec, lambda_path,loc_scr,mod1vars)

    !// IN THIS SUBROUTINE WE TREAT THE CONTRIBUTION OF THREE
    !// AND FOUR SEGMENT LOOPS BY BUILDING THE EXTERNAL SPACE
    !// COMPLEMENT TO A SINGLE SEGMENT IN THE INTERNAL SPACE.
    !// THE CALCULATION IS PERFORMED IN VECTORIZED MODE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!

    use locist_var_mod,only:locist_scratch

    implicit none

    type(orbital_path)::lambda_path
    type(locist_scratch)::loc_scr
    type(threefourmod1vars)::mod1vars
    type(blockedLockVectorType) :: civec, sigmavec

    integer::start_elec                !// NUMBER OF ELECTRONS WHERE WE HAVE LEFT OFF
    integer::a,b                      !// INTEGERS TO LABEL EXTERNAL ORBITALS
    integer::buffer_meter              !// USED TO COUNT INTEGRAL
    integer::lam_or_mu                 !// SET TO 2 TO INDICATE WE ARE BUILDING A MU PATHS
    integer::constraint_count          !// KEEPS TRACK OF CONSTRAINED LEVELS
    integer::top_level                !// LABELS THE TOP AND BOTTOM LEVELS OF THE LOOP
    integer::mu_levels                 !// LABELS THE LEVELS IN THE MU PATH
    integer::path_elecs                !// NUMBER OF ELECTRONS IN THE MU PATH    
    integer::current_vertex            !// VERTICES IN THE MU PATHS
    integer::step_type                 !// STEPS IN THE MU PATH
    integer::ibar                      !// POSITIONS OF IMPORTANT ORBITALS
    integer::i_segment                !// ORBITAL INDICES OF IMPORTANT ORBITALS
    integer::lambda_singles,mu_singles !// SINGLES IN LAMBDA AND MU PATHS
    integer::internal_singles          !// NUMBER OF SINGLY OCCUPIED ORBITALS IN INTERNAL PART OF
    !// CONFIGURATION
    integer::loop_type                 !// STORES THE LOOP TYPE                                   
    integer::lambda_weight,mu_weight   !// WEIGHTS AND ADDRESS OF LAMBDA AND MU PATHS
    integer::nm2_dim                             !// NM2 STANDS FOR N-2, NM1 = N-1
    integer::ab_spin                              !// TO INDEX SPIN FUNCTIONS
    integer::mode_1                                !// SYMMETRIC OR ANTISYMMETRIC SCEPPER
    integer::ab_address,aa_address       !// A IS FOR N-1, AB AND AA FOR N-2, V FOR VALENCE
    integer::ab_start,ab_end
    integer::aa_start,aa_end
    integer::v_address,v_dim
    integer::v_start,v_end

    integer::lam_length,lam_lengthC2        
    integer::a1,b1
    integer::status

    real(real8),parameter::zero=real(0.0,real8),one=real(1.0,real8),sqrt2=real(sqrt(real(2.0,real8)),real8)
    real(real8)::cc

    real(real8), dimension(:,:),allocatable::temp_buf_small  !FOR (IA|IB) AND (IA|IA) INTEGRALS
    ! START AUTOGENERATED INITIALIZATION 
    v_address = 0
    lambda_singles = 0
    v_end = 0
    lam_length = 0
    b = 0
    cc = 0.0
    aa_start = 0
    ab_start = 0
    aa_address = 0
    mode_1 = 0
    lam_lengthc2 = 0
    ab_end = 0
    b1 = 0
    ab_spin = 0
    nm2_dim = 0
    a1 = 0
    ab_address = 0
    a = 0
    v_start = 0
    aa_end = 0
    ! END AUTOGENERATED INITIALIZATION 

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// INITIALIZE VARIABLES FOR BUILDING THE MU PATH
    start_elec = sum(lambda_path%occupations(0:num_internal))
    loop_type = lambda_path%loop_type
    lam_or_mu = 2
    constraint_count = 1
    lambda_path%rt_loop_weight = 0
    top_level = lambda_path%level1
    i_segment = lambda_path%constraints(2,1)

    !// NUMBER OF ELECTRONS IN PATH BEFORE LOOP HEAD
    path_elecs = sum(lambda_path%occupations(0:top_level-1))

    !// NUMBER OF SINGLES IN THE INTERNAL SPACE PART OF PATH
    internal_singles = lambda_path%num_singles
    ibar = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (start_elec /= num_elec-2) return

    !// STORE THE WEIGHT OF THE LAMBDA_PATH
    lambda_path%weight = sum(lambda_path%arc_weights(0:num_internal))

    if (skip_this_internal(lambda_path%weight,start_elec)) return

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


!!!!!!!!!!!!!!!!!!!!!
    !//
    !//  DO MULTIPLICATIONS
    !//
!!!!!!!!!!!!!!!!!!!!!

    if (loop_type ==  5) then 

       !//  LETS SCEP THE (IA|IB) INTEGRAL TO A MATRIX   


       buffer_meter = 0
       do a = 1,num_external-1
          do b = a+1,num_external
             buffer_meter = buffer_meter + 1
             mod1vars%integral_buf(a,b) = iaib(buffer_meter)
             mod1vars%integral_buf(b,a) = mod1vars%integral_buf(a,b)
          enddo
       enddo


       do a = num_internal+1,num_orbitals  
          mod1vars%integral_buf(a-num_internal,a-num_internal) = & 
               iaia(a-num_internal) 
       enddo

   
       lambda_weight  = lambda_path%weight+1
       mu_weight      = lambda_path%rt_loop_weight+1

       lambda_singles = internal_singles + 2
       mu_singles     = internal_singles

       nm2_dim = fsn(lambda_singles) 
       v_dim   = fsn(mu_singles) 

       !// COUPLING COEFFICIENT IS IDENTITY 
       v_address = internal_index_vector0(mu_weight) 
       ab_address= internal_index_vector2(lambda_weight) 
       aa_address= internal_index_vector3(lambda_weight) 


       dim2: if (nm2_dim /=0 .and. v_dim /=0 ) then 

          if ((ab_address > 0 .or. aa_address > 0) .and. v_address > 0) then 

             cc = one/sqrt2   !Factor need to be included (sqrt2*1/2) 

             lam_length = num_allowed_virtuals(lambda_weight,"D") 
             lam_lengthC2 = lam_length*(lam_length-1)/2
             call get_virtuals(lambda_weight,"D", loc_scr%virt_lam_allow)
             
             allocate(temp_buf_small(lam_length,lam_length),stat=status)
             call allocatecheck(status,"temp_buf_small")

             do b=1,lam_length
                b1 = loc_scr%virt_lam_allow(b)
                do a=1,lam_length 
                   a1 = loc_scr%virt_lam_allow(a)
                   temp_buf_small(a,b) = mod1vars%integral_buf(a1-num_internal,b1-num_internal)
                enddo
             enddo

             do ab_spin  = 1,fsn(lambda_singles-2)


                ab_start = ab_address + lam_lengthC2*(ab_spin-1)
                aa_start = aa_address + lam_length*(ab_spin-1)

                v_start = v_address + ab_spin - 1
                v_end =   v_start


                mode_1 = -1
                if (ab_spin <= fsn(lambda_singles-2)) mode_1 = 1


                if (internal_index_vector2(lambda_path%weight +1) < 0) cycle

                if (mode_1 .eq. 1) then
                    call scepper_diag_p1(mod1vars%scep_ci_1,civec%v,ab_start, &
                        lam_length, civec%v,aa_start)
                elseif (mode_1 .eq. -1) then
                    call scepper_diag_m1(mod1vars%scep_ci_1,civec%v,ab_start, &
                        lam_length)
                endif

                !// VALENCE/N-2 INTERACTIONS ..
                sigmavec%v(v_start) = sigmavec%v(v_start)+ cc*sum(mod1vars%scep_ci_1(1:lam_length,1:lam_length)* &
                     temp_buf_small(1:lam_length,1:lam_length))

                !                 !// N-2/VALENCE INTERACTIONS ..

                mod1vars%scep_sigma_1(1:lam_length,1:lam_length) = cc*civec%v(v_start)*temp_buf_small(1:lam_length,1:lam_length)

                if (mode_1 .eq. 1) then
                    call unscepper_diag_p1(mod1vars%scep_sigma_1,sigmavec%v,ab_start,lam_length,&
                        sigmavec%v,aa_start)
                elseif (mode_1 .eq. -1) then
                    call unscepper_diag_m12(mod1vars%scep_sigma_1,sigmavec%v,ab_start,lam_length)
                endif


             enddo
             
             deallocate(temp_buf_small,stat=status)
             call deallocatecheck(status,"temp_buf_small")

          endif

       endif dim2

    endif

  end subroutine extern_one_seg_compl_vec_lmo
  !*****************************************************************
  
  subroutine build_ext_34seg_dryrun(the_path,the_virtuals,abcd_rec,a_ic,b_ic)
  
    implicit none
  
    type(orbital_path),intent(in)::the_path
    integer,dimension(:),intent(inout)::the_virtuals
    integer,intent(inout)::abcd_rec
    integer,intent(in)::a_ic,b_ic
  
    integer::head_weight              !// INTERNAL PATH WEIGHT
    integer::ab_address, aa_address   !// ADDRESS OF THE CONFIGURATIONS
    integer::total_spin               !// TOTAL NUMBER OF SPIN FUNCTIONS
    integer::single_spin              !// NUMBER OF SINGLET COUPLED SPIN FUNCTIONS
    integer::single_spinm2            !// NUMBER OF SINGLEST COUPLED SPIN FUNCTIONS FOR
    !integer::mode_1                   !// MODE FOR DISTINQUISHING SYMMETRIC / ANTISYMMETRIC
    !//   SCEP MATRIX
    integer::num_ab                   !// TOTAL NUMBER OF [AB] = Nv*(Nv+1)/2

    integer::a_check2,b_check2

    integer::icount
    integer::length2,length2C2
    logical::a_check,b_check
    
#ifdef TIGER_USE_OMP
     integer::threadID
     threadID = OMP_get_thread_num()+1
#else
     integer,parameter::threadID = 1
#endif

    !//////////////////////////////////////////////////////////////////////////////////////////////

    head_weight=sum(the_path%arc_weights)

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// CHECK TO SEE IF THE CALCULATION SHOULD BE CARRIED OUT.

    if (skip_this_internal(head_weight, num_elec-2)) return

!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    single_spin         = fsn(the_path%num_singles)
    !single_spinm2       = fsn(the_path%num_singles-2)    ! see below
    total_spin          = fsn(the_path%num_singles+2)
    num_ab              = num_external*(num_external+1)/2


    !BUG... we can go off the end of the array of num_singles < 2 ... not that this value seems to do anything in the code
    !but better safe than sorry
    ! David Krisiloff 5/18/2011
    if ( the_path%num_singles < 2 ) then 
       single_spinm2 = 0 
    else
       single_spinm2 = fsn(the_path%num_singles-2)
    end if

    !// LOCATE THE STARTING SECTION OF THE CI VECTOR.  REMEMBER THE CI
    !// VECTOR SHOULD BE FORMATTED BY SPIN FUNCTIONS AT THIS POINT.
    !// AB_ADDRESS INDEXES CONFIGURATIONS HAVING TWO DIFFERENT
    !// VIRTUALS SINGLY OCCUPIED AND AA_ADDRESS INDEXES CONFIGURATIONS
    !// HAVING ONE VIRTUAL DOUBLY OCCUPIED
    ab_address = internal_index_vector2(head_weight+1)
    aa_address = internal_index_vector3(head_weight+1)


    length2 = num_allowed_virtuals(head_weight+1,"D")
    length2C2 = length2*(length2-1)/2

    the_virtuals = 0
    call get_virtuals(head_weight+1,"D",the_virtuals)

    a_check = .FALSE.
    b_check = .FALSE.

    a_check2 = 0
    b_check2 = 0

    do icount = 1,length2
       if (the_virtuals(icount) .eq. a_ic) then 
          a_check2 = icount
          a_check = .TRUE.
       endif
       if (the_virtuals(icount) .eq. b_ic) then 
          b_check2 = icount
          b_check = .TRUE.
       endif
    enddo

    if (.not.(a_check .AND. b_check) ) return

    abcd_rec = abcd_rec+1
    icount = a_ic*(a_ic-1)/2 + b_ic
    ! store some necessary data for later
    call abcd_seg_writeData(abcd_rec,pseudo_abcd_info,threadID,icount,head_weight,single_spin,single_spinm2,total_spin,a_check2,b_check2)
  
  end subroutine build_ext_34seg_dryrun
  
  subroutine ext_3_4_seg_loops_vec_lmo_res_2(a,c,civec, sigmavec, ab_1,integrals_abcd,abcd_rec,the_virtuals,cho_data,resizeScr)

    !// THIS ROUTINE IS THE MAIN DRIVER ROUTINE FOR TREATING THE THREE AND
    !// FOUR SEGMENT LOOPS WHICH RESIDE ENTIRELY IN THE EXTERNAL SPACE.  HOWEVER,
    !// RATHER THAN LOOPING OVER THE VARIOUS CONFIGURATIONS, THIS ROUTINE LOOPS
    !// OVER THE SPIN FUNCTIONS. CONSEQUENTLY, IT IS SAID TO BE IN VECTORIZED MODE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    
    integer,intent(inout)::ab_1
    real(real8),dimension(:,:),intent(inout)::integrals_abcd
    integer,intent(inout)::abcd_rec
    integer,dimension(:),intent(inout)::the_virtuals
    type(blockedLockVectorType) :: civec, sigmavec
    integer,intent(in)::a,c
    type(cholesky_data),intent(in)::cho_data
    type(fourExtResizeScr),intent(inout):: resizeScr
    
    integer::head_weight              !// INTERNAL PATH WEIGHT
    integer::ab_address, aa_address   !// ADDRESS OF THE CONFIGURATIONS
    integer::total_spin               !// TOTAL NUMBER OF SPIN FUNCTIONS
    integer::single_spin              !// NUMBER OF SINGLET COUPLED SPIN FUNCTIONS
    integer::single_spinm2            !// NUMBER OF SINGLEST COUPLED SPIN FUNCTIONS FOR
    !//   STATES HAVING ONE VIRTUAL DOUBLY OCCUPIED

    integer::a_check2,b_check2
    integer::ab_2,ab_cpd

    integer::length2,length2C2  
    real(real8),parameter::zero=real(0.0,real8),one=real(1.0,real8),minOne=real(-1.0,real8),sqrt2=real(sqrt(real(2.0,real8)),real8),invsqrt2=real(1/sqrt(real(2.0,real8)),real8)

    integer::flag
    integer::status
    
    integer::current_spin,aa_start,ab_start
    real(real8)::rdum
#ifdef TIGER_USE_SLOW_BLAS1
    real(real8)::rdum2
#endif
    real(real8),external::ddot
    
    
#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num()+1
#else
    integer,parameter::threadID = 1
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ab_2 = 0               ! added so that in case of I/O failure ab_2 is not uninitialized  
    abcd_rec = abcd_rec + 1
    call abcd_seg_readData(abcd_rec,pseudo_abcd_info,threadID,ab_2,head_weight,single_spin,single_spinm2,total_spin,a_check2,b_check2,flag)
    if (ab_1 .ne. ab_2) then
       abcd_rec = abcd_rec -1
       ab_1 = ab_2
       return
    endif

    if (flag .ne. 0) then
       ab_1 = 0
       return
    endif

    ab_cpd = (a_check2-1)*(a_check2-2)/2 + b_check2

    !// LOCATE THE STARTING SECTION OF THE CI VECTOR.  REMEMBER THE CI
    !// VECTOR SHOULD BE FORMATTED BY SPIN FUNCTIONS AT THIS POINT.
    !// AB_ADDRESS INDEXES CONFIGURATIONS HAVING TWO DIFFERENT
    !// VIRTUALS SINGLY OCCUPIED AND AA_ADDRESS INDEXES CONFIGURATIONS
    !// HAVING ONE VIRTUAL DOUBLY OCCUPIED
    ab_address = internal_index_vector2(head_weight+1)
    aa_address = internal_index_vector3(head_weight+1)

    if(ab_address <= 0) return

    length2 = num_allowed_virtuals(head_weight+1,"D")
    length2C2 = length2*(length2-1)/2
    
    ! we can save some allocation by caching the biggest one
    if(length2 > resizeScr%length2 .or. length2C2 > resizeScr%length2C2) then
        deallocate(resizeScr%integrals_abcd_mod2,resizeScr%integrals_abcd_mod3,resizeScr%integrals_abcd_mod4,stat=status) !temp
        call deallocatecheck(status,"ext_3_4_seg_loops_vec_lmo_res_2")
        allocate(resizeScr%integrals_abcd_mod2(length2C2),resizeScr%integrals_abcd_mod3(length2C2), &
            resizeScr%integrals_abcd_mod4(length2),stat=status)
        call allocatecheck(status,"ext_3_4_seg_loops_vec_lmo_res_2")
        resizeScr%length2 = length2
        resizeScr%length2C2 = length2C2
    endif
    
#ifndef TIGER_USE_SLOW_BLAS1
    if(total_spin > resizeScr%lengthSigma) then
        if(allocated(resizeScr%comprSigma)) then
            deallocate(resizeScr%comprSigma,stat=status)
            call deallocatecheck(status,"ext_3_4_seg_loops_vec_lmo_res_2_2")
        end if
        allocate(resizeScr%comprSigma(total_spin),stat=status)
        call allocatecheck(status,"ext_3_4_seg_loops_vec_lmo_res_2_2")
        resizeScr%lengthSigma = total_spin
    end if
#endif
    
    call get_virtuals(head_weight+1,"D",the_virtuals)

    call sort_ints_34_vec_lmo(a,c,integrals_abcd,resizeScr%integrals_abcd_mod2,resizeScr%integrals_abcd_mod3, &
        resizeScr%integrals_abcd_mod4, the_virtuals,num_internal,length2,cho_data)

    if (aa_address > 0) then

       if (a_check2/= b_check2) then

#ifndef TIGER_USE_SLOW_BLAS1
          call dgemv('t', length2, single_spin, &
               sqrt2, civec%v(aa_address), length2, &
               resizeScr%integrals_abcd_mod4, 1, &
               zero, resizeScr%comprSigma, 1)

          call dgemv('t', length2C2, single_spin, &
               one, civec%v(ab_address), length2C2, &
               resizeScr%integrals_abcd_mod2, 1, &
               one, resizeScr%comprSigma, 1)
               
          ab_start = ab_address
          do current_spin = 1,single_spin
              rdum = resizeScr%comprSigma(current_spin)
              !$omp atomic
              sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) + rdum
              ab_start = ab_start+length2C2
          enddo
               
          if(total_spin - single_spin .ne. 0) then
            call dgemv('t', length2C2, total_spin-single_spin, &
               minOne, civec%v(ab_address+single_spin*length2C2), length2C2, &
               resizeScr%integrals_abcd_mod3, 1, &
               zero, resizeScr%comprSigma, 1)
               
            do current_spin = single_spin+1,total_spin
                rdum = resizeScr%comprSigma(current_spin-single_spin)
                !$omp atomic
                sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) + rdum
                ab_start = ab_start+length2C2
            enddo
          endif
#else
          aa_start = aa_address
          ab_start = ab_address
          do current_spin = 1,single_spin
               rdum = ddot(length2C2,civec%v(ab_start),1,resizeScr%integrals_abcd_mod2,1)
               rdum2 = ddot(length2,civec%v(aa_start),1,resizeScr%integrals_abcd_mod4,1)
               rdum = rdum + sqrt2*rdum2
               !$omp atomic
               sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) + rdum
               aa_start = aa_start+length2
               ab_start = ab_start+length2C2
          enddo
          
          do current_spin = single_spin+1,total_spin
               rdum = ddot(length2C2,civec%v(ab_start),1,resizeScr%integrals_abcd_mod3,1)
               !$omp atomic
               sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) - rdum
               ab_start = ab_start+length2C2
          enddo
#endif

       else ! (a_check2 == b_check2)
#ifndef TIGER_USE_SLOW_BLAS1
          call dgemv('t',length2, single_spin, &
               one, civec%v(aa_address), length2, &
               resizeScr%integrals_abcd_mod4, 1, &
               zero, resizeScr%comprSigma, 1)

          call dgemv('t',length2C2, single_spin, &
               invsqrt2, civec%v(ab_address), length2C2, &
               resizeScr%integrals_abcd_mod2, 1, &
               one, resizeScr%comprSigma, 1)
               
          aa_start = aa_address
          do current_spin = 1,single_spin
              rdum = resizeScr%comprSigma(current_spin)
              !$omp atomic
              sigmavec%v(aa_start+b_check2-1) =  sigmavec%v(aa_start+b_check2-1) + rdum
              aa_start = aa_start+length2
          enddo
                  
          if(total_spin - single_spin .ne. 0) then
            call dgemv('t',length2C2,total_spin-single_spin, &
               invsqrt2, civec%v(ab_address+single_spin*length2C2), length2C2, &
               resizeScr%integrals_abcd_mod3, 1, &
               zero, resizeScr%comprSigma,1)
               
            do current_spin = single_spin+1,total_spin
                rdum = resizeScr%comprSigma(current_spin-single_spin)
                !$omp atomic
                sigmavec%v(aa_start+b_check2-1) = sigmavec%v(aa_start+b_check2-1) + rdum
                aa_start = aa_address+length2
            enddo               
          endif
#else
          aa_start = aa_address
          ab_start = ab_address
          do current_spin = 1, single_spin
               rdum = ddot(length2C2,civec%v(ab_start),1,resizeScr%integrals_abcd_mod2,1)
               rdum2 = ddot(length2,civec%v(aa_start),1,resizeScr%integrals_abcd_mod4,1)
               rdum = invsqrt2*rdum + rdum2
               !$omp atomic               
               sigmavec%v(aa_start+b_check2-1) =  sigmavec%v(aa_start+b_check2-1) + rdum
               aa_start = aa_start+length2
               ab_start = ab_start+length2C2
          enddo
          
          do current_spin = single_spin+1,total_spin
               rdum = invsqrt2*ddot(length2C2,civec%v(ab_start),1,resizeScr%integrals_abcd_mod3,1)
               !$omp atomic
               sigmavec%v(aa_start+b_check2-1) = sigmavec%v(aa_start+b_check2-1) + rdum
               aa_start = aa_start+length2
               ab_start = ab_start+length2C2
          enddo
#endif
       endif ! (a_check2/= b_check2)

    elseif (aa_address < 0) then

       if (a_check2 /= b_check2) then

#ifndef TIGER_USE_SLOW_BLAS1
          total_spin = fsn(nm2_singles(head_weight+1)+2)
          single_spin = fsn(nm2_singles(head_weight+1))

          call dgemv('t', length2C2, single_spin, &
               one, civec%v(ab_address), length2C2, &
               resizeScr%integrals_abcd_mod2, 1, &
               zero, resizeScr%comprSigma, 1)
               
          ab_start = ab_address
          do current_spin = 1,single_spin
              rdum = resizeScr%comprSigma(current_spin)
              !$omp atomic
              sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) + rdum
              ab_start = ab_start+length2C2
          enddo
               
          if(total_spin - single_spin .ne. 0) then
            call dgemv('t', length2C2, total_spin-single_spin, &
               minOne, civec%v(ab_address+single_spin*length2C2), length2C2, &
               resizeScr%integrals_abcd_mod3, 1, &
               zero, resizeScr%comprSigma, 1)
               
            do current_spin = single_spin+1,total_spin
                rdum = resizeScr%comprSigma(current_spin-single_spin)
                !$omp atomic
                sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) + rdum
                ab_start = ab_start+length2C2
            enddo
          endif
#else
          aa_start = aa_address
          ab_start = ab_address
          do current_spin = 1, single_spinm2
               rdum = ddot(length2C2,civec%v(ab_start),1,resizeScr%integrals_abcd_mod2,1)
               !$omp atomic
               sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) + rdum
               aa_start = aa_start+length2
               ab_start = ab_start+length2C2
          enddo
          
          do current_spin = single_spin+1,total_spin
               rdum = ddot(length2C2,civec%v(ab_start),1,resizeScr%integrals_abcd_mod3,1)
               !$omp atomic
               sigmavec%v(ab_start+ab_cpd-1) = sigmavec%v(ab_start+ab_cpd-1) - rdum
               ab_start = ab_start+length2C2
          enddo
#endif
       endif

    endif ! a_check2 /= b_check2

  end subroutine ext_3_4_seg_loops_vec_lmo_res_2
  
  subroutine sort_ints_34_vec_lmo(a,c,integrals_abcd,integrals_abcd_mod2,integrals_abcd_mod3,integrals_abcd_mod4, &
    the_virtuals,num_internal,length2,cho_data)
        
  implicit none
  
  integer,intent(in)::a,c
  real(real8),dimension(:,:),intent(inout)::integrals_abcd
  real(real8),dimension(:),intent(inout)::integrals_abcd_mod2,integrals_abcd_mod3,integrals_abcd_mod4
  integer,intent(in)::num_internal,length2
  integer,dimension(:),intent(in)::the_virtuals
  type(cholesky_data),intent(in)::cho_data
  
  integer::i,j,counter,idum,idum2
  real(real8)::tmp1,tmp2
  integer::tmp3,tmp4
 
  !real(real8),dimension(num_external,num_external)::integrals_abcd
  
  if(integralDirect .and. fullyIntegralDirect) then
  
  ! perhaps we need to calculate part of the integrals
  counter = 0
  do i=1,length2
    idum = the_virtuals(i)-num_internal
    do j = 1,i-1
      idum2 = the_virtuals(j)-num_internal
      counter = counter + 1
      tmp1 = integrals_abcd(idum2,idum)
      if(tmp1 == magic_number) then
        ! here is our 42: calculate it
        tmp1 = makeOneABCD(a,idum2+num_internal,c,idum+num_internal,cho_data)
        
        ! store it back in our "cache"
        integrals_abcd(idum2,idum) = tmp1
        
      endif
      
      tmp2 = integrals_abcd(idum,idum2)
      if(tmp2 == magic_number) then
        ! here is our 42: calculate it
        tmp2 = makeOneABCD(a,idum+num_internal,c,idum2+num_internal,cho_data)
        
        ! store it back in our "cache"
        integrals_abcd(idum,idum2) = tmp2
        
      endif
      
      ! compare...
      !write(*,*) "DEBUG: idum,idum2,integral1, integral2,ref1,ref2"
      !write(*,*) idum,idum2,integrals_abcd(idum2,idum),integrals_abcd(idum,idum2),integrals_abcd_real(idum2,idum),integrals_abcd_real(idum,idum2)
      
      integrals_abcd_mod2(counter) = tmp1 + tmp2
      integrals_abcd_mod3(counter) = tmp1 - tmp2
    enddo
    
    ! and one last time for the diagonal integral
    tmp1 = integrals_abcd(idum,idum)
    if(tmp1 == magic_number) then
      ! 42: calculate
      tmp1 = makeOneABCD(a,idum+num_internal,c,idum+num_internal,cho_data)
      ! set back in
       integrals_abcd(idum,idum) = tmp1
    endif
    
    integrals_abcd_mod4(i) = tmp1
  enddo
  
  else if(integralDirect .and. .not. sphere_based_integral_truncations .and. .not. fullyIntegralDirect .and. .not. directLowMem) then
  
  ! all integrals filled in already
  counter = 0
  do i=1,length2
    idum = the_virtuals(i)-num_internal
    tmp3 = max(idum,c-num_internal)
    tmp3 = tmp3*(tmp3-1)/2+min(idum,c-num_internal)
    do j = 1,i-1
      idum2 = the_virtuals(j)-num_internal
      tmp4 = max(idum2,c-num_internal)
      tmp4 = tmp4*(tmp4-1)/2+min(idum2,c-num_internal)
      counter = counter + 1
      tmp1 = integrals_abcd(idum2,tmp3)
      tmp2 = integrals_abcd(idum,tmp4)
      integrals_abcd_mod2(counter) = tmp1 + tmp2
      integrals_abcd_mod3(counter) = tmp1 - tmp2
    enddo
    integrals_abcd_mod4(i) = integrals_abcd(idum,tmp3)
  enddo

  else
  
  ! all integrals filled in already (slightly different order)
  counter = 0
  do i=1,length2
    idum = the_virtuals(i)-num_internal
    do j = 1,i-1
      idum2 = the_virtuals(j)-num_internal
      counter = counter + 1
      tmp1 = integrals_abcd(idum2,idum)
      tmp2 = integrals_abcd(idum,idum2)
      integrals_abcd_mod2(counter) = tmp1 + tmp2
      integrals_abcd_mod3(counter) = tmp1 - tmp2
    enddo
    integrals_abcd_mod4(i) = integrals_abcd(idum,idum)
  enddo
  
  endif
  
  end subroutine sort_ints_34_vec_lmo

end module three_and_four_seg_mod_1
