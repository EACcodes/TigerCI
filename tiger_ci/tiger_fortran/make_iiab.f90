! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
! Change the write file number
! Change the file to direct access 
module iiab_mod
use IOBuffer
contains
  subroutine make_iiab(cho_data)
    ! A subroutine to construct the (ii|ab) integrals for the sigmavector routines

    use global_var_mod
    use molecule_var_mod
    use cholesky_structs
    use utilities_mod
    use wp_tov_mod
#ifdef TIGER_USE_OMP
    use omp_lib
#endif

    implicit none

    integer::a,b,c  ! orbital indices
    integer::i
    integer::max_ii
    integer::i1
    integer::ii_count,ii_count2,ii_count3,ab_count
    integer::icount,idum
    integer::ii_pairs,ii_off
    integer::rec_count
    integer::block_len
    integer::allocatestatus,deallocatestatus
    integer,dimension(num_internal)::ii_block
    integer,dimension(:),allocatable::i_ind,ii_step

    real(real8),allocatable,dimension(:,:)::ii_cho
    real(real8),dimension(numcho)::ab_vec
    real(real8)::iiab_int

    real(kind=real8), external::ddot ! BLAS1 (dot product)

    logical::transform_int  ! Decides if we need to transform integrals
    type(cholesky_data)::cho_data
    integer::threadID
#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif

    !*******************************************************************************

    ! Count number of relevant ii pairs 

    ii_pairs = 0  ! Number of kept (ii) pairs
    ii_block = 0  ! Number of (ii|ab) integrals in each (ii) block
    ab_count = 0  ! Don't remove this line because of stepped back counting

    do i = 1,num_internal
       ii_pairs = ii_pairs + 1 

       ii_block(i) = ab_count ! # of (ii|ab) integrals in an (ii) block. Stepped back counting

       ab_count = 0 

       do a = num_internal+1,num_orbitals

          do b = num_internal+1,a
             if (ignorable_pair(b,a) ) cycle

             idum = a*(a-1)/2+b
             idum = cho_data%mo_ind_inv(idum)
             if (idum .eq. 0) cycle   ! no ab_cho on disk 

             ab_count = ab_count + 1

          enddo
       enddo
    enddo

    ! Now we want to make ii_block a cumulative counter
    ! Useful for writing to file as record counter

    icount = 0
    do i = 1,num_internal
       icount = icount + ii_block(i)
       ii_block(i) = icount 
    enddo

    !******************************************************************************
    ! Set a maximum buffer size, 
    max_ii = int(max_mem_ints/numcho)
    max_ii = max_ii - 1

    !******************************************************************************

    allocate(ii_cho(numcho,max_ii),stat=allocatestatus)
    call allocatecheck(allocatestatus,"ii_cho  ")
    allocate(i_ind(max_ii),stat=allocatestatus)
    call allocatecheck(allocatestatus,"i_ind   ")
    allocate(ii_step(max_ii),stat=allocatestatus)
    call allocatecheck(allocatestatus,"ii_step ")

    !******************************************************************************
    ! Construct the (ii|ab) integrals. Buffered version.
    ! Step 1; Read max block of (ii) buffer
    ! Step 2: Go to work constructall (ii|ab) integrals with (ii) block. 
    ! Step 3: Get new (ii) block and repeat 

    transform_int = .false.

    ii_count = 0            ! Offset position for (ii) buffer  block
    ii_count2 = 0           ! Uninteruppted counting in ii loop
    ii_count3 = 0
    rec_count = 0
    ii_off = 0
    block_len = 0

    do i = 1,num_internal

       ! Read in max {ii} block   

       ii_count = ii_count + 1   ! For counting ii_index in (ii) buffer
       ii_count2 = ii_count2 + 1 ! Actual progress of ii loop

       idum = (i)*(i+1)/2        ! Not (i-1)*(i-2)/2+j
       idum = cho_data%mo_ind_inv(idum)   ! ii should be present
       call for_double_buf_readblock(mo_int_no, idum, ii_cho(1:numcho,ii_count), threadID)
       !read(unit=mo_int_no,rec=idum) ii_cho(1:numcho,ii_count)
       i_ind(ii_count) = i

       if (ii_count .eq. max_ii)then ! Buffer is now full or everything has been read into buffer
          transform_int = .true.   ! Go to transform integrals portion because buffer is full
          block_len = max_ii
       endif

       if (ii_count2 .eq. ii_pairs) then ! Everything has been read into buffer
          transform_int = .true.   ! Go to transform integrals portion because buffer is full
          block_len = ii_count 
       endif




       !*************************************************************************************** 
       if (transform_int) then ! (ii) buffer is full. Use them to construct (ii|ab) integrals

          ii_step = 0             ! To keep the offset position of each ia block on disk

          do a = num_internal+1,num_orbitals   ! Note that  we will construct all the (ii|ab) integrals for the current (ii) block

             do b = num_internal+1,a

                !              if (ignorable_pair(b,a) .eq. 0) cycle  !	 NOT HERE

                ab_vec = 0.0D0
                idum = a*(a-1)/2+b
                idum = cho_data%mo_ind_inv(idum)

                if (idum .eq. 0) cycle
                call for_double_buf_readblock(mo_int_no, idum, ab_vec(1:numcho), threadID)
                !read(unit=mo_int_no,rec=idum) ab_vec(1:numcho)

                ii_count3 = 0 
                do c = 1,block_len
                   ii_count3 = ii_count3 + 1 ! Don't misplace this line
                   i1 = i_ind(ii_count3) ! Retrieves the index of i

                   !                 if (sep_orb(i1,a) .eq. 0) cycle
                   !                 if (sep_orb(i1,b) .eq. 0) cycle
                   if (ignorable_pair(b,a) ) cycle 

                   ii_step(ii_count3) = ii_step(ii_count3) + 1 ! Don't misplace this line
                   rec_count = ii_block(i1) + ii_step(ii_count3)
                   iiab_int = ddot(numcho,ii_cho(:,ii_count3),1,ab_vec(:),1)

                   call for_double_buf_writeElement(cd_iiab_no,rec_count,iiab_int,threadID)

                enddo

             enddo
          enddo

          transform_int = .false.  ! Reset to false when (ii) buffer is exhausted
          i_ind = 0 
          ii_cho = 0.0D0
          ii_count = 0

       endif ! transform_int
       !********************************************************************************************

    enddo ! i


    deallocate(ii_cho,stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"ii_cho  ")
    deallocate(i_ind,stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"i_ind   ")
    deallocate(ii_step,stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"ii_step ")

  end subroutine make_iiab
end module iiab_mod
