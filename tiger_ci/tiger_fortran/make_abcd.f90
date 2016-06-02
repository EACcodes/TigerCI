! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module abcd_mod
  use global_var_mod
  use cholesky_structs
  use io_unit_numbers
  use molecule_var_mod
  use fortran_timing
  use molecule_var_mod
  use wp_tov_mod
  use utilities_mod
  use IOBuffer
#ifdef TIGER_USE_OMP
  use three_four_seg_var_mod
  use omp_lib
#endif

contains
  subroutine make_abcd(cho_data)
    ! A subroutine to construct the (ab|cd) integrals for the sigmavector routines
    implicit none

    integer::a,b,c,d ! orbital indices
    integer::max_ac
    integer::c_lead,a_label
    integer::ab_count,ab_count2
    integer::icount
    integer::idum,idum2,idum3,idum4,idum5,idum6,idum7,idum8
    integer::idummy,idummy2,idummy3
    integer::ivar,ivar2,ivar3,ivar4,ivar5,max_abcd
    integer::c_filter,bc_filter
    integer::ac_pairs
    integer::rec_count,rec_count2
    integer::num_d_s,max_a_dom,active_a_counter
    integer::num_a_blocks,num_a_active,unfilled_a_s
    integer::allocatestatus,deallocatestatus
    integer,dimension(num_orbitals):: a_list,a_dom_size
    integer,dimension(num_orbitals,numberOfThreads)::b_mat,d_mat,d_mat2
    integer,dimension(num_external,num_external)::ac_block
    integer,dimension(:),allocatable::a_ind,b_ind

    real(real8)::zero,one
    real(real8),allocatable,dimension(:,:)::ab_cho
    real(real8),allocatable,dimension(:,:,:)::cd_cho_ten,small_ab_cho_ten,small_cd_cho_ten,abcd_block_ten
    integer::numthreads,ten_pointer

    logical::transform_int  ! Decides if we need to transform integrals
    type(cholesky_data)::cho_data

#ifdef TIGER_FINE_TIMES
    type(clock) :: timer
#endif

#ifdef TIGER_USE_OMP
    numthreads = numberOfThreads
#else
    numthreads = 1
#endif

    !*******************************************************************************

    ac_pairs = 0  ! Number of kept (ac) pairs
    ac_block = 0  ! Number of (ab|cd) integrals in each (ac) block
    ab_count = 0 
    idum6 = 0
    idum7= 0
    idum8 = 0

    ! setup the buffered abcd file
#ifdef DEBUG_TIGER
    write(ioOutput,*) "DEBUG: Opening file cd_abcd in buffer."
    flush(ioOutput)
#endif

    call for_double_buf_openfile(for_buf_int_poolID,cd_abcd_no,scratch_directory//'cd_abcd.dat',len(scratch_directory)+11)

#if defined TIGER_USE_OMP && !defined TIGER_IOWA_OMP
                                ! I am not entirely happy with this solution as we violate the "no public state" rule
                                ! but we anyways need to revamp the integrals and a clean transport of this info will then be
                                ! one of the things to take care of
     if(.not.allocated(omp_offsets_abcd_icount)) then
    allocate(omp_offsets_abcd_icount(num_internal+1:num_orbitals), stat=allocatestatus)
    call allocatecheck(allocatestatus,"omp_offsets_abcd_icount")
 endif
#endif

 do a = num_internal+1,num_orbitals
    do c = num_internal+1,num_orbitals
       if (ignorable_pair(c,a) ) cycle
       ac_pairs = ac_pairs + 1 
    enddo
 enddo

 do a = num_internal+1,num_orbitals
#if defined TIGER_USE_OMP && !defined TIGER_IOWA_OMP
     omp_offsets_abcd_icount(a) = idum6
#endif
    do c = num_internal+1,a

       if (ignorable_pair(c,a) ) cycle

       ab_count = 0 

       do b = num_internal+1,num_orbitals
          if (ignorable_pair(b,a) ) cycle
          if (ignorable_pair(b,c) ) cycle

          do d = num_internal+1,num_orbitals

             if (ignorable_pair(d,a) ) cycle
             if (ignorable_pair(d,c) ) cycle
             if (ignorable_pair(d,b) ) cycle

             ab_count = ab_count + 1
             idum6 = idum6 + 1

          enddo
       enddo
       !        ac_block(c-num_internal,a-num_internal) = ab_count
       ac_block(a-num_internal,c-num_internal) = ab_count ! # of (ab|cd) integrals in this (ac) block
       !        idum7 = idum7 + ab_count
    enddo
 enddo

 ! Stepped back counting 

 icount = 0 
 do a = num_internal+1,num_orbitals
    do c = num_internal+1,a
       if (ac_block(a-num_internal,c-num_internal) .eq. 0) cycle
       idum = ac_block(a-num_internal,c-num_internal)
       ac_block(a-num_internal,c-num_internal) = icount 
       icount = icount + idum 
    enddo
 enddo
 !  write(6,*) "Expt2 icount",icount

 !*******************************************

 ! Find the size of each domain corresponding to each a
 ! Find the maximum size of domain amongst all {a}s
 a_dom_size = 0
 max_a_dom = 0

 do a = num_internal+1,num_orbitals
    idum = 0
    do b = num_internal+1,num_orbitals
       if (.not. ignorable_pair(b,a) ) then
          idum = idum + 1
       endif
    enddo

    a_dom_size(a) = idum
    if (max_a_dom .lt. idum) max_a_dom = idum

    if (idum == 0) then
       write(6,*) "Domain size for orbital",a,"is zero.Check it out!"
       call flush(6)
       stop
    endif

 enddo


 !******************************************************************************
 ! Set a maximum buffer size
 max_ac = int(max_mem_ints/numcho)
 max_ac = max_ac - 1

 ! Test if memory assigned is big enough
 if (2*max_a_dom .gt. max_ac) then
    write(6,*) "There is not eneough memory for assembling the (ab|cd) integrals"
    call flush(6)
    stop
 endif

 max_abcd = 0

 !******************************************************************************

 allocate(ab_cho(numcho,max_ac), cd_cho_ten(numcho,max_a_dom,numthreads), small_cd_cho_ten(numcho,max_a_dom,numthreads) ,&
      small_ab_cho_ten(numcho,max_a_dom,numthreads), abcd_block_ten(max_a_dom,max_a_dom,numthreads), stat=allocatestatus)
 call allocatecheck(allocatestatus, "allocations in make_abcd ... if this is failing check the CD MEMORY keyword")
 !  ab_cho = 0.0D0
 !  cd_cho_ten = 0.0D0
 !  small_cd_cho_ten = 0.0D0
 !  small_ab_cho_ten = 0.0D0
 !  abcd_block_ten = 0.0D0

 allocate(a_ind(max_ac),b_ind(max_ac),stat=allocatestatus)
 call allocatecheck(allocatestatus, "a_ind and b_ind in make_abcd.f90")
 !  a_ind = 0
 !  b_ind = 0

 !******************************************************************************
 ! Construct the (ab|cd) integrals. Buffered version.
 ! Step 1; Read max block of (ab) buffer
 ! Step 2: Go to work construct all (ab|cd) integrals with (ab) block. 
 ! Step 3: Get new (ab) block and repeat till all (ab) blocks are exhausted

 transform_int = .false.

 ab_count = 0           
 ab_count2 = 0
 rec_count = 0
 rec_count2 = 0

 zero = 0.0D0
 one = 1.0D0

 ivar5 = 0

 do a = num_internal+1,num_orbitals
    do b = num_internal+1,num_orbitals

       if (ignorable_pair(b,a) ) cycle

       idum = max(a,b)
       idum = idum*(idum-1)/2+min(a,b)
       idum = cho_data%mo_ind_inv(idum)

       ! Read in max {ab} block   

       ab_count = ab_count + 1   ! For counting ab_index in (ab) buffer
       ab_count2 = ab_count2 + 1 ! Actual progress of ab loop

       if (idum .ne. 0) then
          call for_double_buf_readblock(mo_int_no, idum, ab_cho(1:numcho,ab_count), 1)
       endif

       a_ind(ab_count) = a
       b_ind(ab_count) = b

       if (ab_count .eq. max_ac) then ! Buffer is now full or everything has been read into buffer
          transform_int = .true.   ! Go to transform integrals portion because buffer is full

          !       Check number of different a blocks. Be careful of the exception where a has 0 size block.. Anomaly.
          a_list = 0
          num_a_blocks = 1
          idum = a_ind(1) 
          a_list(1)  = idum   ! Keeps a list of i_s
          do idum2 = 1,ab_count
             if (idum /= a_ind(idum2)) then
                num_a_blocks = num_a_blocks + 1
                idum = a_ind(idum2)
                a_list(num_a_blocks) =  idum   
             endif
          enddo

          num_a_active = num_a_blocks ! Update the number of useful a_s 

          ! Consider the case when one cannot fill the entire a{b}
          ! Do boundary testing
          idum = a_ind(ab_count)
          idum2 = b_ind(ab_count)
          idum4 = 0
          do idum3 = idum2+1,num_orbitals
             if (ignorable_pair(idum,idum3) ) cycle
             idum4 = idum4 + 1
          enddo

          unfilled_a_s = 0
          if (idum4 .ge. 1) then    ! There are unfilled a{b}. Find how many of those
             idum2 = a_list(num_a_blocks) ! get the index of the unfilled a ! Compare against a_ind(ab_count)
             idum3 = 0
             do idum = 1, ab_count
                if (idum2 /= a_ind(idum)) cycle
                idum3 = idum3 + 1
             enddo
             unfilled_a_s = idum3      ! This is the number of unfilled a_s
             num_a_active = num_a_blocks - 1 ! Update the actual number of useful a_blocks

          endif

       endif ! endif ab_count = max_ac


       if (ab_count2 .eq. ac_pairs) then ! Everything has been read into buffer
          transform_int = .true.   ! Go to transform integrals portion because buffer is full
          !        write(6,*) "Everything in Buffer"

          !       Check number of different a blocks. Be careful of the exception where a has 0 size block.. Anomaly.
          a_list = 0
          num_a_blocks = 1
          idum = a_ind(1)
          a_list(1) = idum ! Keeps a list of a_s
          do idum2 = 1,ab_count
             if (idum /= a_ind(idum2)) then
                num_a_blocks = num_a_blocks + 1
                idum = a_ind(idum2)
                a_list(num_a_blocks) = idum  
             endif
          enddo
          num_a_active = num_a_blocks ! Update the number of useful a_s

          ! Consider the case when one cannot fill the entire a{b}
          ! For this case all the a_s must be filled. Can't have any unfilled a_s
          unfilled_a_s = 0
          idum = a_ind(ab_count)
          idum2 = b_ind(ab_count)
          idum4 = 0
          do idum3 = idum2+1,num_orbitals
             if (ignorable_pair(idum,idum3)) cycle
             idum4 = idum4 + 1
          enddo

          if (idum4 .ge. 1) then
             write(6,*) "Something seriously wrong!"
             write(6,*) "You have unfilled a_s in (ab|cd)"
             call flush(6)
             stop
          endif

       endif  ! endif ab_count2 = ac_pairs


       !*************************************************************************************** 
       if (transform_int) then ! (ab) buffer is full. Use them to construct (ab|cd) integrals

          !           write(6,*) "num_a_active",num_a_active
          !           write(6,*) "a_list",a_list
          ! Read entire c{d} block which is relevant to at least one of the a_s in the {ab} block
          c_lead = a_list(num_a_active)   ! Set limit for maximum c

#ifdef TIGER_FINE_TIMES
          call start_clock(timer)
#endif

          !$omp parallel &
          !$omp default(none) &
          !$omp private(idum,idum2,c_filter,ten_pointer,num_d_s,active_a_counter,a_label,icount,ivar4,ivar2, &
          !$omp ivar5,bc_filter,idummy3,idummy,idum3,idum5) &
          !$omp shared(a_list,ignorable_pair,numcho,num_internal,c_lead,num_a_active,num_orbitals,cho_data,a_dom_size, &
          !$omp ac_block, one, zero,b_ind,ab_cho,cd_cho_ten,abcd_block_ten,small_ab_cho_ten,small_cd_cho_ten, &
          !$omp b_mat,d_mat,d_mat2,max_a_dom)

          !$omp do &
          !$omp schedule(static)
          do c = num_internal+1,c_lead
#ifdef TIGER_USE_OMP
             ten_pointer = OMP_get_thread_num()+1
#else
             ten_pointer = 1
#endif
             c_filter = 0
             do idum = 1,num_a_active   ! Still have to check at a later stage if a particular i is useful for this b
                idum2 = a_list(idum)
                if ( .not. ignorable_pair(idum2,c)) c_filter = 1
             enddo

             if (c_filter == 0)  cycle ! This entire c block is useless

             ! Now read in the c{d} block
             cd_cho_ten(:,:,ten_pointer) = 0.0D0    ! Size of cd_cho not right
             d_mat(:,ten_pointer) = 0
             num_d_s = 0

             do d = num_internal+1,num_orbitals

                if (ignorable_pair(d,c)) cycle

                idum = max(d,c)
                idum = idum*(idum-1)/2+min(d,c)
                idum = cho_data%mo_ind_inv(idum)

                num_d_s = num_d_s + 1
                d_mat(num_d_s,ten_pointer) = d   ! Keeps a list of d_s in the c{d} block

                if (idum .ne. 0) then
                   call for_double_buf_readblock(mo_int_no, idum, cd_cho_ten(1:numcho,num_d_s,ten_pointer), ten_pointer)
                endif

             enddo  ! enddo d   


             active_a_counter = 1
             do idum = 1,num_a_active  ! Go through all the complete a blocks

                a_label  = a_list(idum)   ! a_label stores the current a value
                icount =  a_dom_size(a_label)

                if (a_label .ge. c) then ! a should not be less than c

                   if (.not. ignorable_pair(a_label,c) ) then ! Otherwise the current c is not useful for this a_block

                      bc_filter = 0                    ! Check for WP between b and c and pack small_ia_cho
                      !small_ab_cho = 0.0D0 // jmd: IMHO we do not need this and spend a lot of time in here
                      b_mat(1:num_orbitals,ten_pointer) = 0
                      do idummy2 = 1,icount
                         idummy3 = b_ind(active_a_counter+idummy2-1)
                         if (ignorable_pair(idummy3,c)) cycle
                         bc_filter = bc_filter + 1
                         small_ab_cho_ten(1:numcho,bc_filter,ten_pointer) = ab_cho(1:numcho,active_a_counter+idummy2-1)
                         b_mat(bc_filter,ten_pointer) = idummy3
                      enddo


                      if (bc_filter .gt. 0) then 
                         ! Here you want to sift through the a{b} and c{d} block for useful pieces.
                         ! Get only the d_s in c{d} which you can find in a{b}
                         idum3 = 0
                         d_mat2(:,ten_pointer) = 0
                         do idum4 = 1,num_d_s
                            idum5 = d_mat(idum4,ten_pointer) ! Get the c value
                            if (ignorable_pair(idum5,a_label) ) cycle  ! WP between a and d
                            idummy = 0            ! Check for wp between a and c
                            do idummy2 = 1,bc_filter
                               idummy3 = b_mat(idummy2,ten_pointer)
                               if (.not. ignorable_pair(idummy3,idum5) ) idummy = 1 
                            enddo
                            if (idummy == 1) idum3 = idum3 + 1     ! Count the useful d_s
                            d_mat2(idum3,ten_pointer) = idum5
                            small_cd_cho_ten(1:numcho,idum3,ten_pointer) = cd_cho_ten(1:numcho,idum4,ten_pointer)
                         enddo

                         ! Now you have the relevant small c{d} block. Construct the (ab|cd) integrals
                         call dgemm('T','N',bc_filter,idum3,numcho,one,small_ab_cho_ten(:,:,ten_pointer),numcho,&
                              small_cd_cho_ten(:,:,ten_pointer),numcho,&
                              zero,abcd_block_ten(:,:,ten_pointer),max_a_dom)
                         ! Write to file 

                         ivar5 = ac_block(a_label-num_internal,c-num_internal)
                         do ivar = 1,idum3
                            ivar2 =  d_mat2(ivar,ten_pointer)
                            do ivar3 = 1,bc_filter
                               ivar4 = b_mat(ivar3,ten_pointer)
                               if (ignorable_pair(ivar4,ivar2) ) cycle
                               ivar5 = ivar5 + 1
                               !max_abcd = max(max_abcd,ivar5)
                               call for_double_buf_writeElement(cd_abcd_no,ivar5,abcd_block_ten(ivar3,ivar,ten_pointer),ten_pointer)
                               !idum7 = idum7 + 1
                            enddo
                         enddo

                      endif ! if wp(c,b)

                   endif  ! if wp (a,c)

                endif ! if a >= c

                active_a_counter = active_a_counter + icount


             enddo ! enddo num_i_active

          enddo ! enddo c
          !$omp end do nowait
          !$omp end parallel

#ifdef TIGER_FINE_TIMES
          write(*,*) "Walltime for relevant make_abcd loop:", get_clock_wall_time(timer)
#endif

          ! There might be some orphan  i{a} who are not used because buffer does not contain the entire block
          ! Shift these blocks to the top of the buffer and adjust for ia_count

          if (unfilled_a_s /= 0) then

             !           write(6,*) "There are leftovers"
             idum2 = ab_count-unfilled_a_s
             do idum = 1,unfilled_a_s
                idum2 = idum2 + 1
                ab_cho(1:numcho,idum) = ab_cho(1:numcho,idum2)
                a_ind(idum) = a_ind(idum2)
                b_ind(idum) = b_ind(idum2)
             enddo
             ab_count = unfilled_a_s  ! Adjust for ab_count as the buffer is not empty
          else
             ab_count = 0 ! jmd: reset the ab_count as this otherwise gets us into out-of-bounds
          endif

          transform_int = .false.

       endif ! transform_int
       !********************************************************************************************

    enddo ! b
 enddo ! a

 deallocate(ab_cho,stat=deallocatestatus)
 deallocate(cd_cho_ten,stat=deallocatestatus)
 deallocate(small_cd_cho_ten,stat=deallocatestatus)
 deallocate(small_ab_cho_ten,stat=deallocatestatus) 
 deallocate(abcd_block_ten,stat=deallocatestatus)

 deallocate(a_ind,stat=deallocatestatus)
 deallocate(b_ind,stat=deallocatestatus)


 !  do a = 1, max_abcd
 !     call for_double_buf_readElement(cd_abcd_no,a,tmp,1)
 !     write(*,*) "icount and integral ",a,tmp
 !  enddo

end subroutine make_abcd

function makeOneABCD(a,b,c,d,cho_data) 
 implicit none
 integer,intent(in)::a,b,c,d
 type(cholesky_data),intent(in)::cho_data

 real(real8)::makeOneABCD

 integer::cho_point1,cho_point2
 real(real8),parameter::zero=real(0.0,real8)
 real(real8),external::ddot

 cho_point1 = max(a,b)
 cho_point1 = cho_point1*(cho_point1-1)/2+min(a,b)
 cho_point1 = cho_data%mo_ind_inv(cho_point1)

 cho_point2 = max(c,d)
 cho_point2 = cho_point2*(cho_point2-1)/2+min(c,d)
 cho_point2 = cho_data%mo_ind_inv(cho_point2)

 ! we have them in memory
 ! check an adapted CS criterion
 if(cho_data%cho_norms(cho_point1)*cho_data%cho_norms(cho_point2) >= integral_threshold) then
    makeOneABCD = ddot(numcho,cho_data%cho_vectors(:,cho_point1),1,cho_data%cho_vectors(:,cho_point2),1)
    !     makeOneABCD = dot_product(cho_data%cho_vectors(:,cho_point1),cho_data%cho_vectors(:,cho_point2))
 else
    makeOneABCD = zero
 endif

end function makeOneABCD

end module
