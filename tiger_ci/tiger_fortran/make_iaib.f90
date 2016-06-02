! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module iaib_mod
use IOBuffer
contains
  subroutine make_iaib(cho_data)
    ! A subroutine to construct the (ia|ib) integrals for the sigmavector routines

    use global_var_mod
    use molecule_var_mod
    use cholesky_structs
    use utilities_mod
    use wp_tov_mod
#ifdef TIGER_USE_OMP
    use omp_lib
#endif

    implicit none

    type(cholesky_data),intent(in)::cho_data

    integer::i,a,b,ia_count,status,offset,ab_count,iaib_count
    real(real8),dimension(:,:),allocatable::ia_cho
    real(real8)::int2

    real(kind=real8), external::ddot ! BLAS1 (dot product)

#ifdef TIGER_USE_OMP
    integer::threadID
    threadID = OMP_get_thread_num()+1
#else
    integer,parameter::threadID = 1
#endif

    allocate(ia_cho(numcho,num_orbitals-num_internal),stat=status)
    call allocatecheck(status,"make_iaib")

    iaib_count = 0
    do i = 1, num_internal

       ia_cho = 0.0D0
       do a = 1, num_external

          if (ignorable_pair(a+num_internal,i) ) cycle

          ia_count = a+num_internal
          ia_count = ia_count*(ia_count-1)/2+i
          ia_count = cho_data%mo_ind_inv(ia_count)

          if ( ia_count .eq. 0 ) cycle
          call for_double_buf_readblock(mo_int_no, ia_count, ia_cho(1:numcho,a), threadID)
       enddo

       offset = num_internal*num_external*(num_external+1)/2 + (i-1)*num_external*(num_external+1)/2
       ab_count = 0

       do a = 1, num_external
          do b = 1,a

             ab_count = ab_count + 1

             if (ignorable_pair(a+num_internal,i) ) cycle
             if (ignorable_pair(b+num_internal,i) ) cycle

             iaib_count = iaib_count + 1
             int2 = ddot(numcho,ia_cho(1,a),1,ia_cho(1,b),1)

             call for_double_buf_writeElement(cd_iaib_no,iaib_count,int2,threadID)

          enddo
       enddo
    enddo

    deallocate(ia_cho,stat=status)
    call deallocatecheck(status,"make_iaib")

  end subroutine make_iaib
end module iaib_mod
