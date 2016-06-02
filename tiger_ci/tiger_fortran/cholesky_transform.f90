! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> transforms the cholesky vectors from the AO to MO basis
!> \author Jeremy Tchwee
!> \author David Krisiloff
!>

module cholesky_transform
  use fortran_timing
  use global_var_mod
  use molecule_var_mod
  use utilities_mod
  use sort_utils
  use molecular_orbital_mod
  use wp_tov_mod
  use cholesky_decomposition
  use cholesky_structs
  use get_integrals
  use math_utils
  use c_io_finterface
#ifdef TIGER_USE_OMP
  use omp_lib
#endif

  implicit none

  private

  public :: transform_MidMem,cholesky_bomb,post_cd_int_file_prep,analyze_transform_mem,load_old_mo_cho_vec

  type trafoScratch
     real(real8), allocatable,dimension(:,:):: MOChoVecs
     real(real8), allocatable, dimension(:) :: AOChoVec
  end type trafoScratch

contains

  !> \brief brief Transforms cholesky vector to the MO basis, \f$  L_{\mu\nu}^J \rightarrow T_{ij}^J \f$
  !>
  !> \param cho_data    The data structure indexing the final MO cholesky vectors
  !> \param AOIndex     The data structure mapping the combined (ij) index to the individual i,j
  !> \param maxOrbIter  Maximum number of MO orbital targets per iteration of the transform
  !> \param activePairs     The number of basis products (ij) that are included in the decomposition after prescreening
  !> \param nCho        The number of cholesky vectors
  !> Each iteration of the transform we read in the set of AO cholesky vectors and transform 
  !> \f$  L_{\mu\nu}^J \rightarrow T_{i:}^J \f$ for a set of i. This should strike a decent balance
  !> between the amount of I/O performed and the amount of memory the intermediate pieces require.
  subroutine transform_MidMem(cho_data, AOIndex, maxOrbIter, activePairs, nCho)
    implicit none

    ! inputs
    type(simplePairIndex), intent(in) :: AOIndex
    integer, intent(in) :: maxOrbIter, activePairs, nCho
    type(cholesky_data), intent(inout) :: cho_data

    ! Variables
    integer :: ierr ! stores return codes
    integer :: iOrb ! loop counters for looping over the MOs
    integer :: last, i, batch
    integer :: numAOPairs, nOrb
#ifdef TIGER_USE_OMP
    integer :: thread
#else
    integer,parameter:: thread=1
#endif
    integer, dimension(:), allocatable :: orbital_targets
    real(kind = real8), allocatable, dimension(:,:,:) :: halfTransformed

    logical :: more_orbitals
    type(trafoScratch),dimension(:),allocatable::scr

#ifdef TIGER_FINE_TIMES
    type(clock) :: timer
#endif

#ifdef TIGER_USE_OMP
    integer::numthreads
    numthreads = numberOfThreads
#else
    integer,parameter::numthreads = 1
#endif

#ifdef TIGER_FINE_TIMES
    call start_clock(timer)
#endif


    ! Memory allocations before the transform begins
    !*************************************************
    numAOPairs = AOIndex%nPair

    allocate(scr(numthreads),stat=ierr)
    call allocatecheck(ierr,"trafoScratch")

    do i = 1,numthreads
       allocate(scr(i)%AOChoVec(numAOPairs), scr(i)%MOChoVecs(numcho, num_orbitals), stat = ierr)
       ! not in principle necessary, but valgrind is complaining about uninitialized values from somewhere
       scr(i)%AOChoVec = 0 
       scr(i)%MOChoVecs = 0 
       call allocatecheck(ierr, "blockAOChoVec")
    enddo

    ! allocate memory for the half transformed vectors
    allocate(orbital_targets(maxOrbIter),HalfTransformed(numcho, number_bas, maxOrbIter), stat=ierr)
    call allocatecheck(ierr, "HalfTransformed")
    ! not in principle necessary, but valgrind is complaining about uninitialized values from somewhere
    HalfTransformed = 0 
    orbital_targets = 0 

    ! Now for the fun part
    !***********************   
    ! We loop over the number of orbitals, each iteration building T^{:}_{iOrb,:}

    ! Each iteration we construct all \f$ T^{:}_{i,j} \f$ for all j and a subset of all i (the targeted molecular orbitals)
    last = 0 
    call get_next_orbital_batch(more_orbitals, last, orbital_targets, nOrb, maxOrbIter)
    batch = 1
    do while (more_orbitals)
    write(*,*) "Inside cholesky AO->MO transform iteration ", batch
#ifdef DEBUG_TIGER
       write(*,*) "   target orbitals = ", orbital_targets(1:nOrb)
#endif
       ! Read in the AO Cholesky vectors and perform a half transform to \f$ T^{:}_{i,:} \f$ for a set of i 
       call first_half_transform(AOIndex, HalfTransformed, scr, orbital_targets, nOrb, nCho, activePairs)
       
       ! For each \f$ T^{:}_{i,:} \f$ transform to each possible \f$ T^{:}_{i,j} \f$       
       !$omp parallel default(none) shared(HalfTransformed, scr, cho_data, nOrb, orbital_targets) private(i,iOrb, thread)   
#ifdef TIGER_USE_OMP
       thread = OMP_get_thread_num() + 1
#endif
       !$omp do schedule(static) 
       do i = 1, nOrb
          iOrb = orbital_targets(i)
          call second_half_transform(HalfTransformed(:,:,i), scr(thread)%MOChoVecs, iOrb)
          call write_MO_cho_vec(scr(thread)%MOChoVecs, iOrb, thread, cho_data)  
       end do
       !$omp end do
       !$omp end parallel

       call get_next_orbital_batch(more_orbitals, last, orbital_targets, nOrb, maxOrbIter)
       batch = batch + 1
    end do

    ! since we aren't truncating CD vectors at the moment, this index is rather easy
    do i = 1, num_orbitals * (num_orbitals+1)/2
       cho_data%mo_ind_inv(i) = i 
    end do

    ! deallocate memory
    do i = 1,numthreads
       deallocate(scr(i)%AOChoVec, scr(i)%MOChoVecs, stat = ierr)
       call deallocatecheck(ierr, "blockAOChoVec")
    enddo
    deallocate(scr, halfTransformed, stat = ierr)
    call allocatecheck(ierr, "deallocations in transform_MOLCAS_MidMem")

#ifdef TIGER_FINE_TIMES
    write(*, *) "Walltime for molcas MidMem transform:", get_clock_wall_time(timer)
#endif
  end subroutine transform_MidMem

  !> brief Opens a number of integral files need later on in the SGGA procedure
  subroutine post_CD_int_file_prep()

    use io_unit_numbers

    implicit none
    integer length_unf
    integer ierror

    !>\todo jmd: EVERY SINGLE ONE OF THOSE EVENTUALLY MUST BE MIGRATED TO THE FORTRAN BUFFER!!!

    ! Files w/ integer records
    inquire(iolength = length_unf) number_bas
    open(unit = iajk_bl_no, file = scratch_directory // 'iajk_bl.dat', access = 'direct', recl = length_unf, form = 'unformatted', iostat = ierror)
    call open_check(ierror)
    open(unit = ikaj_bl_no, file = scratch_directory // 'ikaj_bl.dat', access = 'direct', recl = length_unf, form = 'unformatted', iostat = ierror)
    call open_check(ierror)
    open(unit = ijka_bl_no, file = scratch_directory // 'ijka_bl.dat', access = 'direct', recl = length_unf, form = 'unformatted', iostat = ierror)
    call open_check(ierror)

    ! Files with no specified record size
    open(unit = cd_iapp_no, file = scratch_directory // 'cd_iapp.dat', form = 'unformatted', iostat = ierror)
    call open_check(ierror)
    open(unit = cd_ipap_no, file = scratch_directory // 'cd_ipap.dat', form = 'unformatted', iostat = ierror)
    call open_check(ierror)
    open(unit = rediaia_no, file = scratch_directory // 'rediaia.dat', form = 'unformatted', iostat = ierror)
    call open_check(ierror)
    open(unit = rediaib_no, file = scratch_directory // 'rediaib.dat', form = 'unformatted', iostat = ierror)
    call open_check(ierror)
    open(unit = cd_ijja_no, file = scratch_directory // 'cd_ijja.dat', form = 'unformatted', iostat = ierror)
    call open_check(ierror)
    open(unit = cd_ijia_no, file = scratch_directory // 'cd_ijia.dat', form = 'unformatted', iostat = ierror)
    call open_check(ierror)

    ! they may say test ... but they are used anyway ... sigh 
    !> \todo I probably should rethink these file names
    open(unit = cd_test_no, file = scratch_directory // 'cd_test.dat', form = 'unformatted', iostat = ierror)
    call open_check(ierror)

  end subroutine post_CD_int_file_prep

  !> \brief Checks the return code on a open statement
  subroutine open_check(ierr)
    implicit none
    integer, intent(in) :: ierr
    if (ierr /= 0) then
       write(*, *) "Error opening a file in cholesky_vectors.f90"
    end if
  end subroutine open_check

  !> \brief This subroutine uses MO CD vectors computed from a previous calculation. It allows us to restart
  !> a CD-LMRSDCI calculation more efficiently
  !> This routine requires a file CD_restart to available as well as the finished CD files. The CD_restart file has the format:
  !> (integer) record size of the mo_int.dat file
  !> (integer) number of cholesky vecs (numcho)
  !> (integer array) mo_ind_inv 
  subroutine load_old_mo_cho_vec(cho_data)
    use io_unit_numbers

    implicit none
    type(cholesky_data) :: cho_data
    integer :: ierror
    integer :: cho_vec_rec_size ! size of a record in the cholesky vector file
    logical :: fileExists

    write(*, *) " "
    write(*, *) "*******************************************"
    write(*, *) " Restarting from previous cholesky vectors "

    ! Read in the global data which is normally set by ChoVecMol
    allocate(cho_data%mo_ind_inv(num_orbitals * (num_orbitals + 1)/2), stat = ierror)
    call allocatecheck(ierror, "mo_ind_inv")
    inquire(file = "CD_restart", exist = fileExists)
    if (.not.fileExists) then
       call cholesky_bomb("Could not find file CD_restart. Unable to restart from old cholesky vectors")
    end if
    open(unit = cd_restart_no, file = scratch_directory // "CD_restart", form = 'unformatted')
    read(unit = cd_restart_no) cho_vec_rec_size
    read(unit = cd_restart_no) numcho
    read(unit = cd_restart_no) cho_data%mo_ind_inv
    close(unit = cd_restart_no)
    write(*, *) "   - Found cholesky specific restart information "

    ! Open up the cholesky vector file
    inquire(file = "mo_int.dat", exist = fileExists)
    if (.not.fileExists) then
       call cholesky_bomb("Could not find file mo_int.dat. Unable to restart from old cholesky vectors")
    end if

    ! Open up the integral files we'll need later
    call post_CD_int_file_prep

    write(*, *) "   - Done "
    write(*, *) "*******************************************"
    write(*, *) " "
  end subroutine load_old_mo_cho_vec

  !> \brief Provides a simple produce error message and then crash routine
  !> \param message The message to write out before crashing
  subroutine cholesky_bomb(message)
    use io_unit_numbers

    implicit none
    character(len = *) :: message
    write(*, *) "FATAL ERROR IN CHOLESKY TRANSFORM"
    write(*, *) message
    call flush(6)
    stop
  end subroutine cholesky_bomb

  !> \brief Reads in a set of our Cholesky vectors
  !> \param vec     Storage space for the cholesky vector
  !> \param index   The starting point for reading vectors. We read [start, start+block-1]
  !> \param activePairs     The number of basis products (ij) that are included in the decomposition after prescreening
  subroutine read_our_vec(vec, index, activePairs)
    implicit none

    ! inputs
    real(kind = real8), intent(out) :: vec(:)
    integer, intent(in) :: index, activePairs

    ! important ... zero vec first
    vec = 0
    if (CPP_DECOMPOSED_INTS) then 
        call for_double_buf_readblock(cvec_no, index, vec(:), 1)
    else
        call read_vector_with_c(scratch_directory // 'cvec.dat', vec(:), activePairs, index)
    end if
    vec(1:index-1) = 0.0 ! These elements (everything above the lower triangular) aren't part of the CD

  end subroutine read_our_vec


  !> \brief Transforms the first AO index to an MO index
  !> \param AOIndex     The data structure mapping the combined (ij) index to the individual i,j
  !> \param HalfTransformed An array to hold the half transformed (one AO index, one MO index) vectors
  !> \param scr An object holding the required scratch space: array to hold the CD vectors read in and array to hold the CD vectors read in (transposed from the MOLCAS ordering)
  !> \param orbitals    The different MO orbitals to transform to 
  !> \param nOrb        The number of different MO orbitals to transform to 
  !> \param nCho        Number of cholesky vectors
  !> \param activePairs     The number of basis products (ij) that are included in the decomposition after prescreening
  subroutine first_half_transform(AOIndex, HalfTransformed, scr, orbitals, nOrb, nCho, activePairs)
    implicit none

    ! inputs
    type(simplePairIndex), intent(in) :: AOIndex
    type(trafoScratch),dimension(:),intent(inout)::scr
    real(kind = real8), dimension(:,:,:) :: HalfTransformed

    integer, dimension(:) :: orbitals
    integer, intent(in) :: nOrb, nCho, activePairs
    
    ! local variables
    integer :: iPair
    integer :: mu, nu
    integer :: iOrb, i, iCho
    integer :: thread

    real(kind = real8) :: c_mu, c_nu

    ! Read in block-wise the AO cholesky vectors and half transform them
    HalfTransformed = 0.0

    !$omp parallel &
    !$omp default(none) &
    !$omp private(mu,nu,c_mu,c_nu,iPair,thread,i,iOrb) &
    !$omp shared(nCho,AOIndex,HalfTransformed,molecular_orbitals, &
    !$omp ignorable_pair,numcho,number_bas,orbitals,nOrb,scr,activePairs)
    !$omp do schedule(static)
    do iCho = 1, nCho
#ifdef TIGER_USE_OMP
       thread = OMP_get_thread_num() + 1
#else
       thread = 1
#endif
       !$omp critical
       call read_our_vec(scr(thread)%AOChoVec, iCho, activePairs)
       !$omp end critical

       do i = 1, nOrb
          iOrb = orbitals(i)
          ! Transform the first index
          do iPair = 1, AOIndex%nPair
             mu = AOIndex%pair2ij(iPair)%i
             nu = AOIndex%pair2ij(iPair)%j
             c_mu = molecular_orbitals(mu, iOrb)
             c_nu = molecular_orbitals(nu, iOrb)

             ! Note that we could transform either mu or nu --> iOrb
             if (mu /= nu) then
                 HalfTransformed(iCho, mu, i) = HalfTransformed(iCho, mu, i) + c_nu * scr(thread)%AOChoVec(iPair)
                 HalfTransformed(iCho, nu, i) = HalfTransformed(iCho, nu, i) + c_mu * scr(thread)%AOChoVec(iPair)
             else
                 HalfTransformed(iCho, mu, i) = HalfTransformed(iCho, mu, i) + c_mu * scr(thread)%AOChoVec(iPair)
             end if
          end do ! iPair
       end do ! iOrb
    end do ! iCho
    !$omp end do nowait
    !$omp end parallel
  end subroutine first_half_transform


  !> \brief Takes the half transformed (1 AO, 1MO index) vectors and completes the transformation to the MO basis
  !> \param HalfTransformed  The half transformed vectors
  !> \param MOChoVecs         The final, fully transformed vectors
  !> \param iOrb             The MO index of the half transformed vectors
  subroutine second_half_transform(HalfTransformed, MOChoVecs, iOrb)
    implicit none

    ! inputs
    integer,intent(in) :: iOrb
    real(kind = real8), dimension(:,:),intent(in) :: HalfTransformed
    real(kind = real8), dimension(:,:),intent(inout)::MOChoVecs

    call dgemm('N', 'N', numcho, iOrb, number_bas, 1.0, HalfTransformed, numcho, molecular_orbitals, number_bas, 0.0, MOChoVecs, numcho)
  end subroutine second_half_transform

  !> \brief Writes the final MO cholesky vectors to disk 
  !> \param MOChoVecs  The final, fully transformed vectors
  !> \param iOrb      The first MO index in the MOChoVecs
  !> \param thread    The thread number
  !> \param cho_data  The Cholesky data object
  subroutine write_MO_cho_vec(MOChoVecs, iOrb, thread, cho_data)
    implicit none

    ! inputs
    integer, intent(in) :: iOrb, thread
    real(kind = real8), dimension(:,:), intent(in) :: MOChoVecs
    type(cholesky_data), intent(inout)::cho_data

    ! local variables
    integer :: jOrb, ij

    !the whole formula is ij = iOrb * (iOrb - 1) /2 + jOrb
    ij = iOrb * (iOrb -1)/2 
    ! Write them out to disk
    do jOrb = 1, iOrb
       ij = ij + 1

       cho_data%cho_norms(ij) = sqrt(abs(dot_product(MOChoVecs(:, jOrb),MOChoVecs(:, jOrb))))
       if(cdVecsInMemory) then
          ! keep stuff in memory
          cho_data%cho_vectors(:,ij) = MOChoVecs(:, jOrb)
       endif
       call for_double_buf_writeblock(mo_int_no, ij, MOChoVecs(:, jOrb), thread)
    end do
  end subroutine write_MO_cho_vec

  !> \brief Figures out how much memory we have to work with for the CD transform
  !> \param numAOPairs Number of (ij) AO pairs that exist 
  !> \param numthreads   The number of threads during OMP execution
  !> \param numcho       Number of cholesky vectors
  !> \param nOrbTargets  Maximum number of MO indices targeting during one single half transform
  subroutine analyze_transform_mem(mem, numAOPairs, numthreads, numcho, nOrbTargets)
    implicit none

    ! inputs
    integer, intent(in) ::  numAOPairs, numthreads, mem, numcho
    integer, intent(out) :: nOrbTargets

    ! local variables
    integer :: requirements


    ! Here we figure out exactly how much memory we have to work with
    !
    ! Constant pieces
    ! ----------------
    !  (number_bas + 1) * number_bas / 2               AOIndex
    !  num_orbitals * (num_orbitals + 1)/2             mo_ind_inv
    !  numAOPairs * numthreads                          AOChoVecs
    !  numcho * num_orbitals * numthreads              MOChoVecs
    !
    ! What we can play with
    ! ---------------------
    ! numcho*number_bas*nOrbTargets                    Half Transformed
    ! 

    
    
    requirements = (number_bas + 1) * number_bas / 2 + num_orbitals * (num_orbitals + 1)/2 + &
                   numAOPairs * numthreads + numcho * num_orbitals * numthreads
                   
!#ifdef DEBUG_TIGER
    write(*,*) "Transform memory considerations"
    write(*,*) 'AOIndex', (number_bas + 1) * number_bas / 2 
    write(*,*) 'mo_ind_inv', num_orbitals * (num_orbitals + 1)/2
    write(*,*) 'AOChoVecs', numAOPairs * numthreads
    write(*,*) 'MOChoVecs', numcho * num_orbitals * numthreads
    write(*,*) 'Total required memory in GB = ', requirements * 8 / 1024. /1024. /1024. 
!#endif
    
    nOrbTargets = floor(float((mem - requirements)) / ( numcho * number_bas))
    if (nOrbTargets > num_orbitals) then
       ! We can't target more orbitals to transform than already exist
       nOrbTargets = num_orbitals
    end if

    if (nOrbTargets <= 0 ) then 
100    format(/,"Insufficient memory for the Cholesky Transform"/,"Memory (in # doubles)",/"-------------------")
103    format("1 Cholesky vector          ", I20)
104    format("Additional data structures ", I20)
105    format("Current CD memory          ", I20)
106    format(" [ current memory - additonal data structures ] / size of a cholesky vector <= 0 !!")
       write(*,100)
       write(*,103) numcho*num_orbitals
       write(*,104) requirements
       write(*,105) mem
       write(*,106)
       stop
    end if

    write(*,*) "Number of MO orbital targets per transform iteration = " , nOrbTargets

  end subroutine analyze_transform_mem

  !> \brief Figures out which orbitals we should target during the next transform iteration
  !> \param more_orbitals If true on exit we have more oribtals to transform
  !> \param last_orbital  The index of the last MO index transformed
  !> \param orbitals      On exit contains the indices to transform to 
  !> \param nOrb          On exit contains the number of MO indices to transform to 
  !> \param maxOrb        The maximum number of MO indicies to transform to (constrained by memory)
  subroutine get_next_orbital_batch( more_orbitals, last_orbital, orbitals, nOrb, maxOrb) 
    implicit none

    ! inputs
    integer, intent(in) :: maxOrb
    integer, intent(out) :: orbitals(:), nOrb
    integer, intent(inout) :: last_orbital
    logical, intent(out) :: more_orbitals

    ! local variables
    integer :: count

    if ( last_orbital == num_orbitals ) then 
       more_orbitals = .false. 
    else
       more_orbitals = .true. 
       count = 1 
       nOrb  = 0 
       do while (.true.) 
          if ( count > maxOrb) exit
          if ( count + last_orbital > num_orbitals) exit
          orbitals(count) = last_orbital + count
          nOrb = nOrb + 1
          count = count + 1
       end do
       last_orbital = orbitals(nOrb)
    end if

  end subroutine get_next_orbital_batch

end module cholesky_transform



