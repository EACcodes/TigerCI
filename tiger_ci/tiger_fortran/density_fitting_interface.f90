! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module density_fitting_interface
    use global_var_mod
    use cholesky_transform
    use cholesky_structs
    use molecule_var_mod
    use fetch_globals
    use io_unit_numbers
    use two_electron_integrals
    use c_io_finterface
    use IOBuffer
    implicit none
    
    contains
    
    subroutine setup_density_fitting(cho_data)
        implicit none
        type(cholesky_data), intent(inout) :: cho_data
        integer :: i,npair,ierr
        real(kind=real8), allocatable :: tmp(:)
        
        ! get data about the number of DF integrals and the necessary io buffer info of where they were stored
        call get_global_int("Density Fitting Unit Number", i)
        if (i /= mo_int_no) then
            write(*,*) "Fatal error in density fitting setup ... why is mo_int_no different in C++ and F90?"
            stop
        end if
        call get_global_int("Density Fitting Pool ID", for_buf_cd_poolID)
        call get_global_int("Density Fitting naux", numcho)       ! for the rest of the code numcho is just number of auxiliary functions
        write(*,*) "Number of auxiliary functions = " , numcho
        
        ! open necessary files
        call post_cd_int_file_prep
        npair = num_orbitals * (num_orbitals+1) /2
        
        ! setup cho_data: the norms and if we are doing DIRECT mode also store the DF quantities 
        allocate(cho_data%mo_ind_inv(npair), cho_data%cho_norms(npair), tmp(numcho), stat=ierr)
        call allocatecheck(ierr,  "setup_density_fitting")
        do i = 1, npair
            cho_data%mo_ind_inv(i) = i 
        end do
        if (cdVecsInMemory) then 
           allocate(cho_data%cho_vectors(numcho, npair), stat=ierr)
           call allocatecheck(ierr, "Memory for storing the density fitting integrals in mem")
        end if

        do i = 1, npair
            call for_double_buf_readblock(mo_int_no, i, tmp, 1)
            cho_data%cho_norms(i) = sqrt(abs(dot_product(tmp,tmp)))
            if (cdVecsInMemory) then
               cho_data%cho_vectors(:,i) = tmp(:)
            end if
        end do
                
    end subroutine
    
end module
    
    
    
    
