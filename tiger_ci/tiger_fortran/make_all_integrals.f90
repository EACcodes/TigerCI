! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module make_all_integrals
  use global_var_mod
  use cholesky_structs
  use io_unit_numbers
  use fortran_timing
  use abcd_mod
  use iabc_mod
  use iaib_mod
  use iajb_mod
  use iajk_mod
  use iiab_mod
  use ijab_mod
  use ijka_mod
  use ikaj_mod
  use IOBuffer
  implicit none

contains

  subroutine makeIntegrals(cho_data)
    implicit none

    type(cholesky_data)::cho_data

#ifdef TIGER_FINE_TIMES
    type(clock) :: timer
#endif

    ! open some forbuf files which with the current setup might otherwise cause problems
    call for_double_buf_openfile(for_buf_int_poolID,cd_ijpp_no, &
              scratch_directory // 'cd_ijpp.dat', len(scratch_directory) + 11)
    call for_double_buf_openfile(for_buf_int_poolID,cd_ipjp_no, &
              scratch_directory // 'cd_ipjp.dat',len(scratch_directory) + 11)

    if (.not. integralDirect .or. sphere_based_integral_truncations) then

       write(*,*) "DEBUG: Doing NON-DIRECT (ab|cd) and (ia|bc) integrals."

#ifdef DEBUG_TIGER
       write(ioOutput,*) "DEBUG: Entering abcd."
       flush(ioOutput)
#endif
#ifdef TIGER_FINE_TIMES
       call start_clock(timer)
#endif
       if(.not. fullyIntegralDirect) then
          call make_abcd(cho_data) ! this is internally OMP parallelized and expensive
       endif
#ifdef DEBUG_TIGER
       write(ioOutput,*) "DEBUG: Done with abcd."
#endif
#ifdef TIGER_FINE_TIMES
       write(ioOutput,*) "Walltime for abcd integrals:", get_clock_wall_time(timer)
       call start_clock(timer)
#endif
       call make_iabc(cho_data) ! this is also internally OMP parallelized, yet less expensive
#ifdef DEBUG_TIGER
       write(ioOutput,*) "DEBUG: Done with iabc."
#endif
#ifdef TIGER_FINE_TIMES
       write(ioOutput,*) "Walltime for iabc integrals:", get_clock_wall_time(timer)
#endif

    else
       write(*,*) "DEBUG: Doing DIRECT (ab|cd) and (ia|bc) integrals."
    endif


#ifdef TIGER_FINE_TIMES
    call start_clock(timer)
#endif

    ! open a couple of files to not run into concurrency issues...
    if(.not. integralDirect .or. sphere_based_integral_truncations) then
       call for_double_buf_openfile(for_buf_int_poolID,cd_iajb_unsorted_no,&
          scratch_directory // 'cd_iajb_unsorted.dat',len(scratch_directory) + 20)
       call for_double_buf_openfile(for_buf_int_poolID,cd_iajb_no,&
          scratch_directory // 'cd_iajb.dat',len(scratch_directory) + 11)
       call for_double_buf_openfile(for_buf_int_poolID,cd_ijab_no,&
          scratch_directory // 'cd_ijab.dat',len(scratch_directory) + 11)
    endif
    call for_double_buf_openfile(for_buf_int_poolID,cd_ijka_no,&
          scratch_directory // 'cd_ijka.dat',len(scratch_directory) + 11)
    call for_double_buf_openfile(for_buf_int_poolID,cd_ikaj_no,& 
          scratch_directory // 'cd_ikaj.dat',len(scratch_directory) + 11)
    call for_double_buf_openfile(for_buf_int_poolID,cd_iajk_no,&
          scratch_directory // 'cd_iajk.dat',len(scratch_directory) + 11)
 !   endif
    call for_double_buf_openfile(for_buf_int_poolID,cd_iaib_no,&
          scratch_directory // 'cd_iaib.dat',len(scratch_directory) + 11)
    call for_double_buf_openfile(for_buf_int_poolID,cd_iiab_no,&
          scratch_directory // 'cd_iiab.dat',len(scratch_directory) + 11)

    if(.not. fullyIntegralDirect .and. .not. integralDirect) then
       call for_double_buf_openfile(for_buf_int_poolID,cd_ijkl_no,&
          scratch_directory // 'cd_ijkl.dat',len(scratch_directory) + 11)
    endif

    ! these should all be relatively cheap (and are therefore executed in parallel)
    !$omp parallel &
    !$omp default(none) &
    !$omp shared(cho_data,integralDirect,sphere_based_integral_truncations,fullyIntegralDirect, &
    !$omp cdVecsInMemory)
    !$omp sections
    !$omp section
    if(.not. integralDirect .or. sphere_based_integral_truncations) then
       if(.not. fullyIntegralDirect) then
          call make_iajb(cho_data)
#ifdef DEBUG_TIGER
          write(ioOutput,*) "DEBUG: Done with iajb."
#endif
       endif
    endif
    !$omp section
    if(.not. integralDirect .or. sphere_based_integral_truncations) then
       if(.not. fullyIntegralDirect) then
          call make_ijab(cho_data)
#ifdef DEBUG_TIGER
          write(ioOutput,*) "DEBUG: Done with ijab."
#endif
       endif
    endif
    !$omp section
!    if(.not. integralDirect .or. sphere_based_integral_truncations) then
       call make_iajk(cho_data)
#ifdef DEBUG_TIGER
       write(ioOutput,*) "DEBUG: Done with iajk."
#endif
!    endif
    !$omp section
!    if(.not. integralDirect .or. sphere_based_integral_truncations) then
       call make_ijka(cho_data)
#ifdef DEBUG_TIGER
       write(ioOutput,*) "DEBUG: Done with ijka."
#endif
!    endif
    !$omp section
!    if(.not. integralDirect .or. sphere_based_integral_truncations) then
       call make_ikaj(cho_data)
#ifdef DEBUG_TIGER
       write(ioOutput,*) "DEBUG: Done with ikaj."
#endif
!    endif
    !
    !
    !
    !
    !$omp section
    call make_iiab(cho_data)
#ifdef DEBUG_TIGER
    write(ioOutput,*) "DEBUG: Done with iiab."
#endif
    !$omp section
    call make_iaib(cho_data)
#ifdef DEBUG_TIGER
    write(ioOutput,*) "DEBUG: Done with iaib."
#endif
    !$omp end sections
    !$omp end parallel

    ! remove the unsorted integrals from the buffer
    if(.not.integralDirect .or. sphere_based_integral_truncations) then
       call for_double_buf_removefile(cd_iajb_unsorted_no)
    endif

#ifdef TIGER_FINE_TIMES
    write(ioOutput,*) "Walltime for other integrals:", get_clock_wall_time(timer)
#endif
  end subroutine makeIntegrals

  function makeOneIntegral(a,b,c,d,cho_data,threadID)

    use global_var_mod
    use cholesky_structs
    use io_unit_numbers
    use molecule_var_mod
    use IOBuffer

    implicit none
    integer,intent(in)::a,b,c,d
    type(cholesky_data),intent(in)::cho_data
    integer,intent(in)::threadID

    real(real8)::makeOneIntegral

    integer::cho_point1,cho_point2
    real(real8),parameter::zero=real(0.0,real8)
    real(real8),external::ddot
    real(real8),dimension(:),pointer::ab_cho,cd_cho

    cho_point1 = max(a,b)
    cho_point1 = cho_point1*(cho_point1-1)/2+min(a,b)
    cho_point1 = cho_data%mo_ind_inv(cho_point1)

    cho_point2 = max(c,d)
    cho_point2 = cho_point2*(cho_point2-1)/2+min(c,d)
    cho_point2 = cho_data%mo_ind_inv(cho_point2)

    if(fullyIntegralDirect) then
       ! we have them in memory
       ! check an adapted CS criterion
       if(cho_data%cho_norms(cho_point1)*cho_data%cho_norms(cho_point2) >= integral_threshold) then
          makeOneIntegral = ddot(numcho,cho_data%cho_vectors(:,cho_point1),1,cho_data%cho_vectors(:,cho_point2),1)
          !      makeOneIntegral = dot_product(cho_data%cho_vectors(:,cho_point1),cho_data%cho_vectors(:,cho_point2))
       else
          makeOneIntegral = zero
       endif
    else
       if(cdVecsInMemory) then

          makeOneIntegral = ddot(numcho,cho_data%cho_vectors(:,cho_point1),1,cho_data%cho_vectors(:,cho_point2),1)

       else
          ! ask IO Buffer
          call for_double_buf_get_constpointer(mo_int_no, cho_point1, ab_cho, threadID)
          call for_double_buf_get_constpointer(mo_int_no, cho_point2, cd_cho, threadID)

          makeOneIntegral = ddot(numcho,ab_cho,1,cd_cho,1)

          call for_double_buf_return_pointer(mo_int_no, cho_point1, threadID)
          call for_double_buf_return_pointer(mo_int_no, cho_point2, threadID)
       endif

    endif

  end function makeOneIntegral

  function makeOneIntegralPrescreened(a,b,c,d,cho_data,threadID)

    use global_var_mod
    use cholesky_structs
    use io_unit_numbers
    use molecule_var_mod
    use IOBuffer

    implicit none
    integer,intent(in)::a,b,c,d
    type(cholesky_data),intent(in)::cho_data
    integer,intent(in)::threadID

    real(real8)::makeOneIntegralPrescreened

    integer::cho_point1,cho_point2
    real(real8),parameter::zero=real(0.0,real8)
    real(real8),external::ddot
    real(real8),dimension(:),pointer::ab_cho,cd_cho

    cho_point1 = max(a,b)
    cho_point1 = cho_point1*(cho_point1-1)/2+min(a,b)
    cho_point1 = cho_data%mo_ind_inv(cho_point1)

    cho_point2 = max(c,d)
    cho_point2 = cho_point2*(cho_point2-1)/2+min(c,d)
    cho_point2 = cho_data%mo_ind_inv(cho_point2)
    
    if(cdVecsInMemory) then
       ! we have them in memory
       ! check an adapted CS criterion
       if(cho_data%cho_norms(cho_point1)*cho_data%cho_norms(cho_point2) >= integral_threshold) then
          makeOneIntegralPrescreened = ddot(numcho,cho_data%cho_vectors(:,cho_point1),1,cho_data%cho_vectors(:,cho_point2),1)
       else
          makeOneIntegralPrescreened = zero
       endif
    else
       ! ask IO Buffer
       ! check an adapted CS criterion
       if(cho_data%cho_norms(cho_point1)*cho_data%cho_norms(cho_point2) >= integral_threshold) then
          call for_double_buf_get_constpointer(mo_int_no, cho_point1, ab_cho, threadID)
          call for_double_buf_get_constpointer(mo_int_no, cho_point2, cd_cho, threadID)

          makeOneIntegralPrescreened = ddot(numcho,ab_cho,1,cd_cho,1)

         call for_double_buf_return_pointer(mo_int_no, cho_point1, threadID)
         call for_double_buf_return_pointer(mo_int_no, cho_point2, threadID) 
       else
          makeOneIntegralPrescreened = zero
       end if
    end if
  end function makeOneIntegralPrescreened


end module make_all_integrals

