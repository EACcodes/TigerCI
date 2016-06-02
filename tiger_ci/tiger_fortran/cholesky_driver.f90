! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> \brief Driver for handling the cholesky vectors. Handles the decomposition of the AO 2-electron integrals 
!> and their transform to the MO basis set

module cholesky_driver_mod
    use global_var_mod
    use molecule_var_mod
    use utilities_mod
    use molecular_orbital_mod
    use cholesky_structs
    use cholesky_decomposition
    use cholesky_transform
    use fortran_timing
    use get_basis_data
    use IOBuffer
    use fetch_globals
    use io_unit_numbers
    implicit none
    
    contains
    
    
    !*************************************************************************
    !> \brief Controls which cholesky routines are run
    !> \param restart Restart from previously constructed cholesky vectors
    !> Note that this requires that the MO's already be constructed. It also (counter-intuitively) requires
    !> that the spheres and cylinders already be set up. This isn't an obvious requirement, but due to the fact
    !> that we use them to truncate vectors they really need to be there.
    subroutine cholesky_driver(cho_data, restart, mem)
        implicit none

        ! inputs
        logical, intent(in) :: restart
        integer, intent(in) :: mem
        type(cholesky_data) :: cho_data

        ! local variables 
        integer :: dataSize
        integer :: ierr
        integer :: maxOrbIter
        integer :: nAOPairs
        real(kind = real8), allocatable, dimension(:) :: cholesky_vector
        type(clock) :: transformCDTimer

        type(simplePairIndex) :: AOIndex

#ifndef TIGER_USE_OMP
        integer, parameter :: numberOfThreads = 1
#endif
        
    ! if we are restarting do that
        if (restart) then
            call load_old_MO_cho_vec(cho_data)
            return
        end if

        ! Run the two-electron integral cholesky decomposition
        !*************************************************************************************

        ! Note that if we are doing the CD ourselves, we do it in our_setup
        ! if cpp_decomposed_ints==true only prepares all quantities for
        ! transformation as decomposition done in C++
        call OUR_CD(AOIndex, numcho, mem, nAOPairs)


        ! Final memory allocations, file handling, etc... before the transform begins
        !****************************************************************************
        call start_clock(transformCDTimer)
        allocate(cho_data%mo_ind_inv(num_orbitals * (num_orbitals + 1)/2), &
           cho_data%cho_norms(num_orbitals * (num_orbitals + 1)/2),stat = ierr)
        call allocatecheck(ierr, "mo_ind_inv")
        cho_data%mo_ind_inv = 0
        if(cdVecsInMemory) then
           allocate(cho_data%cho_vectors(numcho,num_orbitals * (num_orbitals + 1)/2),stat = ierr)
           call allocatecheck(ierr, "mo_ind_inv")
        endif

        ! open the file to store the finished cholesky vectors
        allocate(cholesky_vector(1:numcho))
        inquire(iolength = dataSize) cholesky_vector(1:numcho) ! figure out how long this data structure should be
        !open(unit=mo_int_no, file= scratch_directory // 'mo_int.dat',form='unformatted',access='direct',recl=dataSize, iostat=ierr)
        deallocate(cholesky_vector)

        ! For restarting later ... if for any reason you are SURE you won't restart you could skip this step
        open(unit = cd_restart_no, file = scratch_directory // "CD_restart", form = 'unformatted')
        write(unit = cd_restart_no) dataSize

        ! Write out some summary information on the amount of data we are about to process
        write(*, *) " "
        140 format(/, 1x, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", /&
        , 1x, "!//                              ", /&
        , 1x, "!// ENTERED AO-> MO CHOLESKY     ", /&
        , 1X, "!//     VECTOR TRANSFORM         ", /&
        , 1x, "!//                              ", /&
        , 1x, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", /)
        write(*, 140)
        write(*, *) "Number of cholesky vectors.......", numcho
        write(*, *) "Number of (ij) pairs (AO)........", AOIndex%nPair

        write(*, *) "Creating buffer for transformed CD vectors of maximum size ", for_buf_maxmemCD
        call for_double_buf_construct(for_buf_maxmemCD/(numcho * real8), numcho, 0, numberOfThreads, for_buf_cd_poolID)
        if (restart) then
            call for_double_buf_reopenfile(for_buf_cd_poolID, mo_int_no, &
               scratch_directory // 'mo_int.dat', len(scratch_directory) + 10)
        else
            call for_double_buf_openfile(for_buf_cd_poolID, mo_int_no, &
               scratch_directory // 'mo_int.dat', len(scratch_directory) + 10)
        endif
        write(*, *)


        ! Here is where we call the routine that will do the actual math ....
        !*********************************************************************
        call analyze_transform_mem(mem, AOIndex%nPair, numberOfThreads, numcho, maxOrbIter)
        call transform_MidMem(cho_data, AOIndex, maxOrbIter, nAOPairs, numcho)


        ! We finished :) lets just clean up and do some quick set up for the rest of the code
        !*************************************************************************************
        write(*, *) "Finished AO->MO Cholesky vector transform"
        write(*, *) " "
        call flush(6)

        ! Print out data for restarting from the finished cholesky vectors
        write(unit = cd_restart_no) numcho
        write(unit = cd_restart_no) cho_data%mo_ind_inv
        close(unit = cd_restart_no)

        ! Open a bunch of integral files for later
        call post_CD_int_file_prep()

        ! Print out timing info
        call print_clock(transformCDTimer, "AO->MO Cholesky Transform")

        ! memory deallocations
        deallocate(AOIndex%pair2ij)

    end subroutine cholesky_driver

    !> \brief Sets up our cholesky decomposition of the two-electron integrals
    !> \param AOIndex   atomic orbital pair index
    !> \param nCho      Number of cholesky vectors
    !> \param mem      The total memory available for cholesky related things (in # real(kind=8))
    !> \param nAOPairs   The number of (ij) pairs after prescreening (forms the basis of the (ij|kl) matrix which was decomposed)
    subroutine OUR_CD(AOIndex, nCho, mem,nAOPairs)
        implicit none

        ! inputs
        integer, intent(in) :: mem
        integer, intent(out) :: nCho
        type(simplePairIndex), intent(inout) :: AOIndex
        integer,allocatable :: pair(:)
        
        ! local variables
        integer :: ierr
        integer :: i,unit_no
        integer :: nShell
        integer :: base
        integer, intent(out) :: nAOPairs
        integer, allocatable, dimension(:) :: nFuncInShell, bas2shell


        ! if integrals were not decomposed in C++ do it here
        if (.NOT. CPP_DECOMPOSED_INTS) then
            ! First step, we need to gather some basic basis set information
            allocate( nFuncInShell(number_bas), bas2shell(number_bas), stat = ierr)
            nFuncInShell = 0
            bas2shell = 0
            call allocatecheck(ierr, "arrays for basis_func_data")
            call nr_shells(nShell)
            call basis_func_data(nFuncInShell, number_bas, nShell)
            do i = 1, nShell
                base = sum(nFuncInShell(1:i - 1))
                bas2shell(base + 1: base + nFuncInShell(i)) = i
            end do
    
            ! now run our cholesky decomposition
            call new_cholesky(nShell, nFuncInShell, bas2shell, number_bas, nCho, cd_thresh, mem, nAOPairs)
        ! if integrals decomposed in C++ do here the AO-MO transformation 
        else 
            write(*,*) "Decomposition done in C++. Only prepare AO-MO transformation"
            call get_global_int("Nr of active Pairs", nAOPairs)
            call get_global_int("Density Fitting naux", numcho)
        end if
        
        ! read in the mu and nu index
        allocate(AOIndex%pair2ij(nAOPairs), stat=ierr)
        call allocatecheck(ierr, "AOIndex%pair2ij")
        AOIndex%nPair = nAOPairs
        
        allocate(pair(2))
        ! if already decomposed in C++ indices are in IO Buffer
        if (cpp_decomposed_ints) then
            call get_global_int("Pair Index Unit Number", unit_no)
            if (unit_no /= cvec_index_no) then
                write(*,*) "Fatal error in density fitting setup ... why is cvec_index_no different in C++ and F90?"
                stop
            end if

            do i = 1, nAOPairs
                call for_int_buf_readblock(cvec_index_no, i, pair, 1)
                AOIndex%pair2ij(i)%i = pair(1)
                AOIndex%pair2ij(i)%j = pair(2)
            end do
        ! if decomposition in Fortran, the IO Buffer was not used
        else
            do i = 1, nAOPairs
!               call setup_vector_io_with_c(scratch_directory // 'cvec_index.dat')
                call read_ints_with_c(scratch_directory // 'cvec_index.dat', pair, 2, i)
                AOIndex%pair2ij(i)%i = pair(1)
                AOIndex%pair2ij(i)%j = pair(2)
            enddo
        endif
        deallocate(pair)

        ! clean up memory if decomposition done here
        if (.NOT. cpp_decomposed_ints) then
            deallocate(nFuncInShell)
        endif
    end subroutine OUR_CD
end module
