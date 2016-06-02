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
!//  GET_INTEGRAL_MOD - THIS SUBROUTINE CONTAINS A NUMBER OF ROUTINES
!//  RESPONSIBLE FOR GETTING INTEGRALS.  I HAVE INCLUDED THIS AS A 
!//  SEPARATE MODULE TO INSURE THAT IT IS EASY TO PORT THIS PROGRAM
!//  TO OTHER SOFTWARE PACKAGES WHICH GENERATE THE INTEGRALS DIFFERENTLY
!//  LATER ON.  ALL THE ROUTINES WILL CHECK TO SEE IF THE ARRAY
!//  IS CURRENTLY ALLOCATED BEFORE ALLOCATING SPACE TO FETCH THE INTEGRALS
!// 
!//  WRITTEN BY DEREK WALTER, 1999
!//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************

module get_integrals_mod

  use global_var_mod
  use integral_storage_mod
  use molecule_var_mod
  use utilities_mod
  use wp_tov_mod
  use IOBuffer
#ifdef TIGER_USE_OMP
  use omp_lib
#endif

  implicit none

  contains

  !*****************************************************************
  subroutine getpp

    !// THIS SUBROUTINE GETS THE ONE ELECTRON INTEGRALS.

    implicit none

    integer::num_pairs              !// NUMBER OF 1E- INTEGRALS
    integer::allocatestatus         !// FOR DYNAMIC MEMORY

    !// IF THE ARRAY IS ALREADY ALLOCATED, JUST EXIT
    if (allocated(pp)) return

    num_pairs = num_orbitals*(num_orbitals + 1)/2

    !// ALLOCATE THE ARRAY FOR THE ONE ELECTORN INTEGRALS
    allocate(pp(num_pairs), stat = allocatestatus)
    call allocatecheck(allocatestatus,"pp      ")

    !// NOW WE CAN READ IN FROM THE APPROPRIATE INTEGRAL FILE THE DATA WE NEED.
    rewind ioii
    read(ioii) pp

  end subroutine getpp

  !*************************************************************
  subroutine getpppp_cho(cho_data)

    !// THIS SUBROUTINE GETS THE (PP|PP) INTEGRALS FROM CHOLESKY VECTORS

    use cholesky_structs
    implicit none

    integer::i,idum
    integer::allocatestatus
!    real(kind=real8) :: diff
    real(real8), allocatable, dimension(:)::pppp_vec
    real(kind=real8), external::ddot
    type(cholesky_data)::cho_data
    integer::threadID

#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif

    
    !// ALLOCATE ARRAY FOR (PP|PP)
    allocate(pppp_vec(numcho))


    if (allocated(pppp)) return 

    allocate(pppp(num_orbitals),stat = allocatestatus)
    call allocatecheck(allocatestatus,"pppp    ")

    pppp_vec = 0.0D0
    pppp = 0.0D0

    do i = 1, num_orbitals
       idum = i*(i+1)/2
       idum = cho_data%mo_ind_inv(idum)
       call for_double_buf_readblock(mo_int_no, idum, pppp_vec(1:numcho), threadID)
       pppp(i) = ddot(numcho,pppp_vec,1,pppp_vec,1)
       !call check_integral( pppp(i),i,i,i,i)
    enddo

    deallocate(pppp_vec)
  end subroutine getpppp_cho
  !*************************************************************
  subroutine getpprr_cho(cho_data)

    !// THIS SUBROUTINE GETS THE (PP|RR) INTEGRALS FROM CHOLESKY VECTORS

    use cholesky_structs
    implicit none

    integer::i,j,icount,idum
    integer::allocatestatus
!    real(real8) :: diff
    real(real8),allocatable,dimension(:,:)::pprr_vec
    real(kind=real8), external::ddot
    type(cholesky_data)::cho_data
    integer::threadID

#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif

    allocate(pprr_vec(numcho,num_orbitals))

    if (allocated(pprr)) return

    allocate(pprr(num_orbitals*(num_orbitals - 1)/2),stat=allocatestatus)
    call allocatecheck(allocatestatus,"pprr    ")

    pprr = 0.0D0
    pprr_vec = 0.0D0

    do i = 1, num_orbitals
       idum = i*(i+1)/2
       idum = cho_data%mo_ind_inv(idum)

       if (idum .ne. 0) then
          call for_double_buf_readblock(mo_int_no, idum, pprr_vec(1:numcho,i), threadID)
       endif

    enddo

    ! Construct the (pp|rr) integrals

    icount = 0
    do i = 1,num_orbitals
       do j = 1,i-1
          icount = icount + 1
          pprr(icount) = ddot(numcho,pprr_vec(1,i),1,pprr_vec(1,j),1)
          !call check_integral( pprr(icount), i,i,j,j)
       enddo
    enddo

    deallocate(pprr_vec)

  end subroutine getpprr_cho
  !**************************************************************
  subroutine getprpr

    !// THIS SUBROUTINE GETS THE (PR|PR)  INTEGRALS.

    implicit none

    integer::allocatestatus
    integer::deallocatestatus
    integer::nc2              !// NUMBER OF ORBITALS CHOSE 2
    real(real8), dimension(:), allocatable::tempprpr   !// FOR TEMPORARY STORAGE OF INTEGRALS

    !// IF THE ARRAY IS ALREADY ALLOCATED, JUST EXIT
    if (allocated(prpr)) return


    nc2 = num_orbitals*(num_orbitals - 1)/2

    allocate(prpr(nc2),stat=allocatestatus)
    call allocatecheck(allocatestatus,"pprr    ")
    prpr = 0

    allocate(tempprpr(2*nc2), stat=allocatestatus)
    call allocatecheck(allocatestatus,"tempprpr")
    tempprpr = 0

    rewind ioiijj
    read(ioiijj) tempprpr

    prpr = tempprpr(nc2+1:2*nc2:1)

    deallocate(tempprpr,stat=deallocatestatus)
    call deallocatecheck(deallocatestatus,"tempprpr")

  end subroutine getprpr
  !****************************************************************
  subroutine getprpr_cho(cho_data)

    !// THIS SUBROUTINE GETS THE (PR|PR) INTEGRALS FROM CHOLESKY VECTORS

    use cholesky_structs
    implicit none

    integer::i,j,icount,idum
    integer::allocatestatus
!    real(real8) :: diff
    real(real8),allocatable,dimension(:)::prpr_vec
    real(kind=real8), external::ddot
    type(cholesky_data)::cho_data
    integer::threadID

#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif
    
    allocate(prpr_vec(numcho))
    prpr_vec = 0.0

    !// ALLOCATE ARRAY FOR (PR|PR)

    if (allocated(prpr)) return

    if (.not.allocated(prpr)) then
       allocate(prpr(num_orbitals*(num_orbitals-1)/2),stat=allocatestatus)
       call allocatecheck(allocatestatus,"prpr    ")
    endif

    prpr = 0.0D0

    icount = 0
    do i = 1, num_orbitals
       do j = 1,i-1  

          prpr_vec = 0.0D0
          icount = icount + 1
          idum = i*(i-1)/2+j
          idum = cho_data%mo_ind_inv(idum)

          if ( idum .ne. 0) then
             call for_double_buf_readblock(mo_int_no, idum, prpr_vec(1:numcho), threadID)
          endif

          prpr(icount) = ddot(numcho,prpr_vec,1,prpr_vec,1)
          !call check_integral(prpr(icount), i,j,i,j)

       enddo
    enddo

    deallocate(prpr_vec)
  end subroutine getprpr_cho
  !****************************************************************
  subroutine get_ijaj_cho(cho_data,i,j,iteration)

    use cholesky_structs

    implicit none

    integer::a,i,j
    integer::ierror
    integer::ij_count,ia_count,ja_count
    integer::iteration
    integer::allocatestatus

    real(real8)::int1,int2
    real(real8),allocatable,dimension(:)::cho1,cho2
    real(kind=real8), external::ddot
    type(cholesky_data)::cho_data
    integer::threadID

#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif
    
    if (.not.allocated(integral_buffer)) then
       allocate(integral_buffer(3,integral_buffer_size),stat=allocatestatus)
       call allocatecheck(allocatestatus,"integral")
    endif

    integral_buffer = 0.0D0


    if (iteration .eq. 1) then
       allocate(cho1(numcho),cho2(numcho))

       cho1 = 0.0D0
       ij_count = i*(i-1)/2+j
       ij_count = cho_data%mo_ind_inv(ij_count)

       if (ij_count .ne. 0) then
          call for_double_buf_readblock(mo_int_no, ij_count, cho1, threadID)
       endif

       do a = 1,num_orbitals-num_internal

          if (ignorable_pair(i,a+num_internal) .or. ignorable_pair(j,a+num_internal) ) cycle 

          cho2 = 0.0D0
          ja_count = a+num_internal
          ja_count = ja_count*(ja_count-1)/2+j
          ja_count = cho_data%mo_ind_inv(ja_count)

          if (ja_count .ne. 0) then
             call for_double_buf_readblock(mo_int_no, ja_count, cho2, threadID)
             int1 = ddot(numcho,cho1,1,cho2,1) 
             integral_buffer(2,a) = int1 
          endif
          write(unit=360,iostat=ierror) int1
       enddo

       do a = 1,num_orbitals-num_internal

          if (ignorable_pair(i,a+num_internal) .or. ignorable_pair(j,a+num_internal) ) cycle

          cho2 = 0.0D0
          ia_count = a+num_internal
          ia_count = ia_count*(ia_count-1)/2+i
          ia_count = cho_data%mo_ind_inv(ia_count)

          if (ia_count .ne. 0) then
             call for_double_buf_readblock(mo_int_no, ia_count, cho2, threadID)
             int2 = ddot(numcho,cho1,1,cho2,1)
             integral_buffer(1,a) = int2
          endif
          write(unit=361,iostat=ierror) int2
       enddo
       deallocate(cho1,cho2)
    else

       do a = 1,num_orbitals-num_internal

          if (ignorable_pair(i,a+num_internal)  .or. ignorable_pair(j,a+num_internal)) cycle

          read(unit=360) integral_buffer(2,a)
          read(unit=361) integral_buffer(1,a)
       enddo

    endif

  end subroutine get_ijaj_cho
!  !*****************************************************************
  subroutine get_ipap_cho(cho_data,i,iteration)

    use cholesky_structs

    implicit none

    integer::icount,ip_count,ia_count,pp_count,ap_count
    integer::i,i1,a,p
    integer::allocatestatus
    integer::iteration
    integer::ierror

    real(real8),allocatable,dimension(:,:)::cho_1
    real(real8),allocatable,dimension(:)::cho_2
    real(real8)::int1,int2
    real(kind=real8), external::ddot
    type(cholesky_data)::cho_data
    integer::threadID
    
#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif
    
    ! START AUTOGENERATED INITIALIZATION 
    icount = 0
    p = 0
    a = 0
    int1 = 0.0
    int2 = 0.0
    ierror = 0
    ip_count = 0
    i1 = 0
    ap_count = 0
    pp_count = 0
    ! END AUTOGENERATED INITIALIZATION 

    if (.not.allocated(integral_buffer)) then
       allocate(integral_buffer(2,integral_buffer_size),stat=allocatestatus)
       call allocatecheck(allocatestatus,"integral")
    endif

    i1 = i

    if (iteration .eq. 1) then
    
       integral_buffer = 0 
       allocate(cho_2(numcho),cho_1(numcho,num_orbitals))

       cho_1 = 0.0D0
       do p = 1, num_internal
          pp_count = p*(p+1)/2
          pp_count = cho_data%mo_ind_inv(pp_count)  ! pp_cho should be present
          call for_double_buf_readblock(mo_int_no, pp_count, cho_1(1:numcho,p), threadID)
       enddo

       do a = num_internal+1,num_orbitals

          if (ignorable_pair(a,i1) ) cycle 

          ia_count = a*(a-1)/2+i1
          ia_count = cho_data%mo_ind_inv(ia_count)

          if (ia_count .eq. 0) cycle
          call for_double_buf_readblock(mo_int_no, ia_count, cho_1(1:numcho,a), threadID)
       enddo

       ! Form the (ia|pp)

       do a = num_internal+1,num_orbitals
          if (ignorable_pair(a,i1) ) cycle
          do p = 1, num_internal

             !       if (ignorable_pair(i1,p) .eq. 0 .and. ignorable_pair(i1,a) .eq. 0) cycle
             !       if (ignorable_pair(a,p) .eq. 0) cycle

             icount = (a-num_internal-1)*num_orbitals+p
             int1 = ddot(numcho,cho_1(1,a),1,cho_1(1,p),1)
             integral_buffer(1,icount) = int1
             write(unit=335,iostat=ierror) int1
          enddo
       enddo

       ! Now deal with (ip|ap)

       cho_1 = 0.0D0
       cho_2 = 0.0D0

       ! Get the ip cho
       do p = 1, num_internal

          !if (sep_orb(p,i1) .eq. 0) cycle

          ip_count = max(p,i1)  
          ip_count = ip_count*(ip_count-1)/2+min(p,i1)
          ip_count = cho_data%mo_ind_inv(ip_count)

          if (ip_count .eq. 0) cycle
          call for_double_buf_readblock(mo_int_no, ip_count, cho_1(1:numcho,p), threadID)
       enddo

       ! Get the ap cho

       do a = num_internal+1,num_orbitals

          if (ignorable_pair(a,i1) ) cycle 

          cho_2 = 0.0D0

          do p = 1,num_internal

             !if (sep_orb(p,i1) .eq. 0 ) cycle
             !       if (ignorable_pair(a,p) .eq. 0) cycle

             ap_count = a*(a-1)/2+p
             ap_count = cho_data%mo_ind_inv(ap_count)

             if (ap_count .eq. 0) cycle
             call for_double_buf_readblock(mo_int_no, ap_count, cho_2(1:numcho), threadID)

             !Form the (ip|ap)

             icount = (a-num_internal-1)*num_orbitals+p
             int2 = ddot(numcho,cho_1(1,p),1,cho_2(1),1)
             integral_buffer(2,icount) = int2
             write(unit=336,iostat=ierror) int2
          enddo

       enddo
       
       deallocate(cho_1,cho_2)

    else   ! subsequent iterations

       ! First deal with (ia|pp)
       ! read the (ia|pp)

       integral_buffer = 0.0D0

       do a = num_internal+1,num_orbitals
          if (ignorable_pair(a,i1)) cycle
          do p = 1, num_internal

             !       if (ignorable_pair(i1,p) .eq. 0 .and. ignorable_pair(i1,a) .eq. 0) cycle
             !       if (ignorable_pair(a,p) .eq. 0) cycle

             icount = (a-num_internal-1)*num_orbitals+p
             read (unit=335) integral_buffer(1,icount)
          enddo
       enddo

       ! Now deal with (ip|ap)
       ! read the (ip|ap)
       do a = num_internal+1,num_orbitals
          if (ignorable_pair(a,i1) ) cycle
          do p = 1,num_internal

             ap_count = a*(a-1)/2+p
             ap_count = cho_data%mo_ind_inv(ap_count)
             if (ap_count .eq. 0) cycle

             !if (sep_orb(p,i1) .eq. 0 ) cycle
             !       if (ignorable_pair(a,p) .eq. 0) cycle

             icount = (a-num_internal-1)*num_orbitals+p
             read(unit=336) integral_buffer(2,icount)
          enddo
       enddo

    endif

  end subroutine get_ipap_cho

subroutine make_ipjp(cho_data,i,j,counter)

    use cholesky_structs

    implicit none
    
    integer,intent(in)::i,j,counter
    type(cholesky_data),intent(in)::cho_data

    integer::ij_count,ip_count,jp_count,pp_count
    integer::p
    integer::allocatestatus

    real(real8),allocatable,dimension(:,:)::cho_1,cho_2
    real(real8),allocatable,dimension(:):: ij_cho,pp_cho
    real(real8)::int1,int2
    real(kind=real8), external::ddot
    integer::threadID

#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif
        
    if (.not.allocated(integral_buffer)) then
       allocate(integral_buffer(3,integral_buffer_size),stat=allocatestatus)
       call allocatecheck(allocatestatus,"integral")
    endif
    integral_buffer = 0
    
    ! allocate some stuff
    allocate(cho_1(numcho,num_internal),cho_2(numcho,num_internal))
    allocate(ij_cho(numcho),pp_cho(numcho))
    
    ! For the (ij|pp)

    ! Get the ij cho
    integral_buffer(1,1:num_orbitals) = 0.0D0

    ij_cho = 0.0D0
    ij_count= i*(i-1)/2+j
    ij_count = cho_data%mo_ind_inv(ij_count)

    if (ij_count .ne. 0) then
       call for_double_buf_readblock(mo_int_no, ij_count, ij_cho, threadID)
    endif

    ! Get the pp cho
    pp_cho = 0.0D0
    do p = 1, num_internal
       pp_count = p*(p+1)/2
       pp_count = cho_data%mo_ind_inv(pp_count)  ! pp_cho should be present
       call for_double_buf_readblock(mo_int_no, pp_count, pp_cho, threadID)
       ! Form the ijpp
       int1 = ddot(numcho,ij_cho,1,pp_cho,1)
       integral_buffer(1,p) = int1
       call for_double_buf_writeelement(cd_ijpp_no,counter*num_internal+p,int1,threadID)
    enddo


    ! For the (ip|jp)

    integral_buffer(2,1:num_orbitals) = 0.0D0
    ! Get the ip_cho
    cho_1 = 0.0D0
    do p = 1, num_internal
       !     if (ignorable_pair(i,p) .eq. 0 .and. ignorable_pair(j,p) .eq. 0) cycle
       ip_count = max(p,i)
       ip_count = ip_count*(ip_count-1)/2+min(p,i)
       ip_count = cho_data%mo_ind_inv(ip_count)

       if (ip_count .eq. 0) cycle
       call for_double_buf_readblock(mo_int_no, ip_count, cho_1(:,p), threadID)
    enddo

    cho_2 = 0.0D0
    do p = 1, num_internal
       !     if (ignorable_pair(i,p) .eq. 0 .and. ignorable_pair(j,p) .eq. 0) cycle
      jp_count = max(p,j)
      jp_count = jp_count*(jp_count-1)/2+min(p,j)
      jp_count = cho_data%mo_ind_inv(jp_count)

      if (jp_count .eq. 0) cycle
      call for_double_buf_readblock(mo_int_no, jp_count, cho_2(:,p), threadID)
    enddo

    do p = 1,num_internal
      !     if (ignorable_pair(i,p) .eq. 0 .and. ignorable_pair(j,p) .eq. 0) cycle
      int2 = ddot(numcho,cho_1(:,p),1,cho_2(:,p),1)
      integral_buffer(2,p) = int2
      call for_double_buf_writeelement(cd_ipjp_no,counter*num_internal+p,int2,threadID)
    enddo
       
    deallocate(cho_1,cho_2,ij_cho,pp_cho)
  end subroutine make_ipjp
  
  subroutine read_ipjp(buffer,counter)
  implicit none
  
  real(real8),dimension(:,:),allocatable,intent(inout)::buffer
  integer,intent(in)::counter
  integer::p
  integer::threadID

#ifdef TIGER_USE_OMP
    threadID = OMP_get_thread_num()+1
#else
    threadID = 1
#endif
  
  ! For the (ij|pp)
  buffer(1,1:num_orbitals) = 0.0D0
  do p = 1, num_internal
     ! Read the ijpp
     !     if (ignorable_pair(i,p) .eq. 0 .and. ignorable_pair(j,p) .eq. 0) cycle
     call for_double_buf_readelement(cd_ijpp_no,counter*num_internal+p,buffer(1,p),threadID)
     !read(unit=330) buffer(1,p)
  enddo

  ! For the (ip|jp)
  buffer(2,1:num_orbitals) = 0.0D0
  do p = 1,num_internal
     ! Read the (ip|jp)
     !     if (ignorable_pair(i,p) .eq. 0 .and. ignorable_pair(j,p) .eq. 0) cycle
     call for_double_buf_readelement(cd_ipjp_no,counter*num_internal+p,buffer(2,p),threadID)
     !read(unit=331) buffer(2,p)
  enddo
  
  end subroutine read_ipjp
 

end module get_integrals_mod
