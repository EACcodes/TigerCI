! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> A really simple set of functions for handling sparse vector
!> 
!> \author David Krisiloff
!> \date   2/4/2013


module sparse_mod
    use global_var_mod
    use utilities_mod
    
    !>  A sparse vector, stores the non-zero values in values(:) with indicies in index(:)
    !>  For example if V(4) = 3 is the only non-zero entry in V then index(1)=4, values(1)=3
    type sparse_vector
        integer :: length, total_length
        integer, allocatable :: index(:)
        real(kind=real8), allocatable :: values(:)
    end type
    
    contains
    
    !> \brief Deallocates the memory of a sparse vector
    !> \param v  The sparse vector
    subroutine destroy_sparse_vector(v)
        implicit none
        type(sparse_vector) :: v
        integer :: ierr
        deallocate(v%values, v%index, stat=ierr)
        call deallocatecheck(ierr, "deallocation of a sparse vector")
        v%length = 0
    end subroutine destroy_sparse_vector
    
    !> \brief Adds a new number to a sparse vector
    !> \param val  The new value to add to the vector
    !> \param pos  The position of the new value
    !> \param v    The sparse vector I'm adding to
    subroutine add_to_sparse_vector(val, pos, v)
        implicit none
        type(sparse_vector) :: v
        integer :: pos
        real(kind=real8) :: val
        
        if ( v%length == v%total_length ) then
            ! Out of room ... allocate more space
            call double_sparse_vector(v)
        end if
        v%index(v%length + 1) = pos
        v%values(v%length + 1) = val
        v%length = v%length + 1
    end subroutine add_to_sparse_vector
    
    !> \brief Doubles the size of the sparse vector
    !> \param v  The sparse vector
    subroutine double_sparse_vector (v)
        implicit none
        type(sparse_vector) :: v
        type(sparse_vector) :: tmp
        ! build a temporary vector to hold the information in this one while we reallocate
        call create_sparse_vector(tmp, v%total_length)
 
        ! copy everything from v to tmp
        tmp%length = v%length
        tmp%index(1:v%total_length) = v%index(1:v%total_length)
        tmp%values(1:v%total_length) = v%values(1:v%total_length)
        
        ! deallocate v and then reallocate it to twice the size
        call destroy_sparse_vector(v)
        call create_sparse_vector(v, tmp%total_length*2)
        
        ! copy everything from tmp back to v
        v%length = tmp%length
        v%index(1:tmp%total_length) = tmp%index(1:tmp%total_length)
        v%values(1:tmp%total_length) = tmp%values(1:tmp%total_length)
        
        ! clean up tmp
        call destroy_sparse_vector(tmp)
    end subroutine double_sparse_vector
    
    !> \brief Creates a new sparse vector
    !> \param v The unallocated vector
    !> \param n The initial size
    subroutine create_sparse_vector(v, n)
        implicit none
        integer :: n
        integer :: ierr
        type(sparse_vector) :: v
        
        v%length = 0 
        v%total_length = n
        allocate(v%index(n),v%values(n),stat=ierr)
        call allocatecheck(ierr, "creating a sparse vector")
    end subroutine create_sparse_vector
        
    !> \brief Takes a spare vector and decompresses it into a non-sparse vector
    !> \param sparse The sparse vector
    !> \param full_vector   The non-sparse vector
    subroutine inflate_sparse_vector(sparse, full_vector)
        implicit none
        type(sparse_vector) sparse
        real(kind=real8) :: full_vector(:)
        integer :: i
        full_vector = 0.0
        do i = 1 , sparse%length
            full_vector( sparse%index(i) ) = sparse%values(i)
        end do
    end subroutine inflate_sparse_vector
    
    !> \brief  A simple print out of a sparse vector
    !> \param v The sparse vector
    subroutine print_sparse_vector(v)
        implicit none
        type(sparse_vector) :: v
        integer :: i 
        do i = 1, v%length
            write(*,*) "element " , v%index(i) , " has value " , v%values(i)
        end do
    end subroutine print_sparse_vector
    
    !> \brief Copies sparse vector a to sparse vector b. Does not check that both are the correct size!!!!
    !> \param a Sparse vector a 
    !> \param b Sparse vector b
    subroutine unsafe_copy_sparse_vector(a, b)
        implicit none
        type(sparse_vector) :: a,b
        integer :: n
        n = a%length
        b%length = n
        b%index(1:n) = a%index(1:n)
        b%values(1:n) = a%values(1:n)
    end subroutine unsafe_copy_sparse_vector
    
    !> \brief Copies spare vector a to sparse vector b. This will resize b is its too small
    !> \param a Sparse vector a 
    !> \param b Sparse vector b
    subroutine copy_sparse_vector(a, b)
        implicit none
        type(sparse_vector), intent(in) :: a
        type(sparse_vector) :: b
        integer :: n 
        
        n = a%length
        if ( n > b%total_length ) then
            ! b is too small, need to delete it and reallocate it
            call destroy_sparse_vector(b)
            call create_sparse_vector(b, n)
        end if
        b%length = n
        b%index(1:n) = a%index(1:n)
        b%values(1:n) = a%values(1:n)
    end subroutine copy_sparse_vector
        
        
        
        
        
    
end module sparse_mod
