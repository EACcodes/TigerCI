! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> Some utility functions for looking at CI Hamiltonians. Take care when using these
!> functions. You can quickly get ridiculous wall times if you try and build the complete
!> Hamiltonian for anything but the smallest of systems


module hamiltonian_utils
    use blocked_locks_mod
    use molecule_var_mod
    use global_var_mod
    use io_unit_numbers
    use random
    use sort_utils
    use c_sort_finterface
    use SGGA, only: htimesc
    use sgga_struct_mod, only: sgga_struct
    use graph_var_mod, only: total_csfs
    use c_io_finterface
    implicit none
contains

    !> \brief Builds the entire Hamiltonian matrix column by column (stored on disk)
    !> \param d  An sgga struct (memory state required for an SGGA calculation)
    subroutine build_a_hamiltonian(d)
        implicit none
        type(sgga_struct), intent(inout) :: d
        integer :: i
        type(blockedLockVectorType) :: civec, sigmavec
        character(len=15) :: fname = "hamiltonian.dat"

        ! allocate the needed memory
        call allocateBlockedLockVector(civec, total_csfs, blockSizeCI)
        call allocateBlockedLockVector(sigmavec, total_csfs, blockSizeSigma)
        civec%v = 0
        sigmavec%v = 0
        
        ! open the file to contain the final data
        !open(unit=ioH, file=scratch_directory // "hamiltonian.dat", form="unformatted")        
        call setup_vector_io_with_c(fname)
        ! build the entire Hamiltonian, one column at a time
        do i = 1, total_csfs
            civec%v = 0 
            civec%v(i) = 1.0
            sigmavec%v = 0 
            call htimesc(civec, sigmavec, d)
            !write(ioH) sigmavec%v
            call c_write_vector(fname, i, sigmavec%v, total_csfs)
        end do
        call clean_up_io_with_c()
        ! congrats you survived
        write(*,*) "Congrats !!! You have built an entire CI hamiltonian without breaking your computer !!!"
        close(unit=ioH)
    end subroutine build_a_hamiltonian
    
    
    !> \brief Randomly samples the diagonal dominance of the CI Hamiltonian. 
    !> \param n  the number of columns to sample
    !> \param d  An sgga struct (memory state required for an SGGA calculation)
    !>
    !> This routine selects a random sample of columns of H and prints out 
    !>  - H_ii
    !>  - \sum H_ij i/=j
    !>  - the ratio of the first two (H_ii/\sum H_ij i/=j)
    !>
    !> Since the CI Hamiltonian is diagonally dominant the ratio is always > 1. 
    !> But in at least one particular case, I'm interested by how much.
    subroutine diag_dom_sample(n,d)
        implicit none
        type(sgga_struct), intent(inout) :: d
        integer, intent(in) :: n 
       
        integer :: i, col, ierr
        integer, dimension(:), allocatable :: columns
        real(kind=real8) :: ratio
        real(kind=real8), dimension(:), allocatable :: diag, sum_not_diag
        type(blockedLockVectorType) :: civec, sigmavec
        
        ! allocate the needed memory
        call allocateBlockedLockVector(civec, total_csfs, blockSizeCI)
        call allocateBlockedLockVector(sigmavec, total_csfs, blockSizeSigma)
        civec%v = 0
        sigmavec%v = 0
        
        allocate(columns(total_csfs), stat=ierr)
        call allocatecheck(ierr, "columns in diag_dom_sample")
        do i = 1, total_csfs
            columns(i) = i 
        end do
        allocate(diag(n), sum_not_diag(n))
        
        ! randomly shuffle the list of columns
        call shuffle_int(columns)
        
        ! calculate the first n columns and process
        101 format ("Column ", I10 , " H_ii = " , E15.8, " sum(H(:,i)) - H(i,i) = ", E15.8, " ratio = " , E15.8) 
        do i = 1, n 
            col = columns(i)
            civec%v = 0 
            civec%v(col) = 1.0
            sigmavec%v = 0 
            call htimesc(civec, sigmavec, d)
            diag(i) = abs(sigmavec%v(col))
            ! before summing we sort the array. might be overkill, but if 
            ! we have a large number of really small numbers and one or two large ones
            ! then this might improve the accuracy
            sigmavec%v(col) = 0.0
            !call hybrid_quicksort_real(sigmavec%v, 1, total_csfs)
            call sort_real_array(sigmavec%v)
            sigmavec%v = abs(sigmavec%v)
            sum_not_diag(i) = sum(sigmavec%v) 
        end do
        
        ! summary of the results
        do i = 1, n
            col = columns(i)
            ratio = diag(i) / sum_not_diag(i)
            write(*,101) col, diag(i), sum_not_diag(i), ratio    
        end do  
        
        ! you finished congrats
        102 format("Congrats !!! You have sampled the CI hamiltonian without breaking your computer !!!")
        write(*,102)
        deallocate(columns)
    end subroutine diag_dom_sample


end module hamiltonian_utils
