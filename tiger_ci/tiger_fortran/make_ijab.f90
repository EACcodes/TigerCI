! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module ijab_mod
use IOBuffer
contains
  subroutine make_ijab(cho_data)
    ! A subroutine to construct the (ij|ab) integrals for the sigmavector routines

    use global_var_mod
    use molecule_var_mod
    use cholesky_structs
    use wp_tov_mod
    use io_unit_numbers
#ifdef TIGER_USE_OMP
    use omp_lib
#endif

    implicit none

    integer::a,b,c  ! orbital indices
    integer::i,j
    integer::max_ab
    integer::a1,b1
    integer::ab_count,ab_count2,ab_count3
    integer::icount,idum,idum2,idum3
    integer::ij_pairs,ab_pairs
    integer::rec_count
    integer::block_len
    integer::allocatestatus,deallocatestatus

    integer,dimension(num_internal,num_internal)::ij_block
    integer,dimension(:),allocatable::a_ind,b_ind,ij_step,itempvec

    real(real8),allocatable,dimension(:)::tempvec
    real(real8),allocatable,dimension(:,:)::ab_cho
    real(real8),dimension(numcho)::ij_vec

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
    ! We are going to create a matrix that relates to whether an (a,b) pair shairs a common virtual domain

    ! Count number of relevant ij pairs 

    ij_pairs = 0  ! Number of kept (ij) pairs
    ij_block = 0  ! Number of (ij|ab) integrals in each (ij) block
    ab_count = 0  ! Don't remove this line because of stepped back counting

    do i = 1,num_internal
       do j = 1,i-1

          !        idum = i*(i-1)/2+j
          !        if (pprr(idum) .lt. 1.0D-8) cycle
          !        if (mo_ind_inv(idum) .eq. 0) cycle
          if (ignorable_pair(i,j) ) cycle

          !        ij_block(i,j) = ab_count ! # of (ij|ab) integrals in an (ij) block. Stepped back counting
          ! Don't mess up the position of this lin
          ij_pairs = ij_pairs + 1 

          ab_count = 0 

          do a = num_internal+1,num_orbitals
             if (ignorable_pair(i,a) ) cycle
             if (ignorable_pair(j,a) ) cycle

             do b = num_internal+1,a

                !              icount = a*(a-1)/2+b
                !              if (mo_ind_inv(icount) .eq. 0) cycle

                if (ignorable_pair(i,b) ) cycle
                if (ignorable_pair(j,b) ) cycle
                if (ignorable_pair(a,b) ) cycle ! Not this

                ab_count = ab_count + 1
             enddo
          enddo
          ij_block(i,j) = ab_count ! # of (ij|ab) integrals in an (ij) block
       enddo
    enddo

    ! Now we want to make ij_block a cumulative counter
    ! Useful for writing to file as record counter

    icount = ij_block(2,1)
    ij_block(2,1) = 0
    do i = 3,num_internal
       do j = 1,i-1
          if (ij_block(i,j) .eq. 0) cycle
          idum = ij_block(i,j)
          ij_block(i,j) = icount 
          icount = icount + idum 
       enddo
    enddo

    ab_pairs = 0
    do a = num_internal+1,num_orbitals
       do b = num_internal+1,a
          if (ignorable_pair(b,a) ) cycle
          ab_pairs = ab_pairs + 1
       enddo
    enddo

    !******************************************************************************
    ! Set a maximum buffer size
    max_ab = int(max_mem_ints/numcho)
    max_ab = max_ab - 1

    !******************************************************************************

    allocate(ab_cho(numcho,max_ab),stat=allocatestatus)
    allocate(a_ind(max_ab),stat=allocatestatus)
    a_ind = 0
    allocate(b_ind(max_ab),stat=allocatestatus)
    b_ind = 0
    allocate(ij_step(max_ab),stat=allocatestatus)
    ij_step = 0

    allocate(tempvec(max_ab),stat=allocatestatus)

    allocate(itempvec(max_ab),stat=allocatestatus)


    !******************************************************************************
    ! Construct the (ij|ab) integrals. Buffered version.
    ! Step 1; Read max block of (ij) buffer
    ! Step 2: Go to work constructall (ij|ab) integrals with (ij) block. 
    ! Step 3: Get new (ij) block and repeat 

    transform_int = .false.

    ab_count = 0            ! Offset position for (ij) buffer  block
    ab_count2 = 0           ! Uninteruppted counting in ij loop
    ab_count3 = 0
    rec_count = 0
    block_len = 0

    do a = num_internal+1,num_orbitals
       do b = num_internal+1,a           ! Not i. Think about how to handle (ii|ab)

          if (ignorable_pair(a,b) ) cycle

          icount = a*(a-1)/2+b
          !        if (pprr(idum) .lt. 1.0D-8) cycle

          idum = cho_data%mo_ind_inv(icount)
          !        if (idum .eq. 0) cycle


          ! Read in max {ab} block   

          ab_count = ab_count + 1   ! For counting ab_index in (ab) buffer
          ab_count2 = ab_count2 + 1 ! Actual progress of ab loop

          if ( idum .ne. 0) then
             call for_double_buf_readblock(mo_int_no, idum, ab_cho(1:numcho,ab_count), threadID)
             !read(unit=mo_int_no,rec=idum) ab_cho(1:numcho,ab_count)
          endif

          a_ind(ab_count) = a
          b_ind(ab_count) = b

          if (ab_count .eq. max_ab)then ! Buffer is now full or everything has been read into buffer
             transform_int = .true.   ! Go to transform integrals portion because buffer is full
             block_len = max_ab
          endif

          if (ab_count2 .eq. ab_pairs) then ! Everything has been read into buffer
             transform_int = .true.   ! Go to transform integrals portion because buffer is full
             block_len = ab_count 
          endif


          !*************************************************************************************** 
          if (transform_int) then ! (ab) buffer is full. Use them to construct (ij|ab) integrals

             idum3 = 0     
             do i= 1,num_internal   ! Note that  we will construct all the (ij|ab) integrals for the current (ab) block
                do j = 1,i-1 

                   if (ignorable_pair(i,j)) cycle

                   ij_vec = 0.0D0
                   icount = i*(i-1)/2+j
                   idum = cho_data%mo_ind_inv(icount)

                   if (idum .ne. 0) then
                      call for_double_buf_readblock(mo_int_no, idum, ij_vec(1:numcho), threadID)
                      !read(unit=mo_int_no,rec=idum) ij_vec(1:numcho)
                   endif

                   ab_count3 = 0 
                   idum2 = 0
                   itempvec = 0
                   tempvec = 0.0D0
                   rec_count = ij_block(i,j)

                   idum3 = idum3 + 1

                   do c = 1,block_len
                      ab_count3 = ab_count3 + 1 ! Don't misplace this line
                      a1 = a_ind(ab_count3) ! Retrieves the index of a
                      b1 = b_ind(ab_count3)
                      if (ignorable_pair(a1,i) ) cycle
                      if (ignorable_pair(b1,i) ) cycle
                      if (ignorable_pair(a1,j) ) cycle
                      if (ignorable_pair(b1,j) ) cycle

                      idum2 = idum2 + 1
                      ij_step(idum3) = ij_step(idum3) + 1 ! Don't misplace this line
                      itempvec(idum2) = rec_count + ij_step(idum3)

                      tempvec(idum2) = ddot(numcho,ab_cho(:,ab_count3),1,ij_vec,1)
                   enddo

                   do c = 1,idum2
                      idum = itempvec(c)
                      call for_double_buf_writeElement(cd_ijab_no,idum,tempvec(c),threadID)
                   enddo

                enddo ! enddo j
             enddo  ! enddo i 

             transform_int = .false.  ! Reset to false when (ij) buffer is exhausted
             a_ind = 0 
             b_ind = 0
             ab_cho = 0.0D0
             ab_count = 0

          endif ! transform_int
          !********************************************************************************************

       enddo ! b
    enddo ! a


    deallocate(ab_cho,stat=deallocatestatus)
    deallocate(a_ind,stat=deallocatestatus)
    deallocate(b_ind,stat=deallocatestatus) 
    deallocate(ij_step,stat=deallocatestatus)
    deallocate(tempvec,stat=deallocatestatus)
    deallocate(itempvec,stat=deallocatestatus)

  end subroutine make_ijab

  function makeOneIJAB(i,j,a,b,cho_data)

    use global_var_mod
    use cholesky_structs
    use io_unit_numbers
    use molecule_var_mod

    implicit none

    integer,intent(in)::i,j,a,b
    type(cholesky_data),intent(in)::cho_data
    real(real8)::makeOneIJAB

    integer::cho_point1,cho_point2
    real(real8),parameter::zero=real(0.0,real8)
    real(real8),external::ddot

    cho_point1 = max(i,j)
    cho_point1 = cho_point1*(cho_point1-1)/2+min(i,j)
    cho_point1 = cho_data%mo_ind_inv(cho_point1)

    cho_point2 = max(a,b)
    cho_point2 = cho_point2*(cho_point2-1)/2+min(a,b)
    cho_point2 = cho_data%mo_ind_inv(cho_point2)

    ! we have them in memory
    ! check an adapted CS criterion
    if(cho_data%cho_norms(cho_point1)*cho_data%cho_norms(cho_point2) >= integral_threshold) then
       makeOneIJAB = ddot(numcho,cho_data%cho_vectors(:,cho_point1),1,cho_data%cho_vectors(:,cho_point2),1)
       !      makeOneIJAB = dot_product(cho_data%cho_vectors(:,cho_point1),cho_data%cho_vectors(:,cho_point2))
    else
       makeOneIJAB = zero
    endif

  end function makeOneIJAB


end module ijab_mod
