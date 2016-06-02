! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! The structures used in the cholesky decomposition code
! Placed here so the definition can be used in multiple 
! cholesky modules (by are nicely seperated from the other 
! parts of tiger)

module cholesky_structs
    use global_var_mod  
    use utilities_mod
    implicit none
    
    ! A single basis function (ij) pair
    type pair
     integer :: i,j
    end type pair

    ! The index I use to keep track of how each (ij) -> i,j and how each i,j -> (ij)
    type pairIndex
        integer :: nPair                                  ! number of pairs
     type(pair), dimension(:), allocatable :: pair2ij     ! maps the pair index to i,j
     integer, dimension(:,:),  allocatable :: ij2pair     ! maps i,j to the pair index (ij)
    end type pairIndex
    
    ! A one way indexing for pair number -> (ij). Useful because the transform doesn't need the inverse mapping
    type simplePairIndex
        integer :: nPair
        type(pair), dimension(:), allocatable :: pair2ij     ! maps the pair index to i,j
    end type simplePairIndex
        
    ! An index for the MO cholesky vectors. This data structure is frequently used by parts of TigerCI outside
    ! the cholesky codes
    type cholesky_data
     integer,dimension(:),allocatable:: mo_ind_inv
     real(real8),dimension(:),allocatable::cho_norms
     real(real8),dimension(:,:),allocatable::cho_vectors
    end type cholesky_data
    
  !> \brief Contains all the half transformed cholesky vectors \f$ T^{:}_{iOrb,:} \f$ for a particular \f$ iOrb \f$
  type half_transformed_block
      integer :: MOIndex                                                !> The transformed MO index \f$ iOrb \f$
      real(kind=real8), dimension(:,:), allocatable :: T(:,:)           !> The \f$ T^{:}_{iOrb,:} \f$
  end type
  
!  !> \brief Keeps track of the block reading of cholesky vectors
!  type choleskyBlocking
!      integer :: nBlocks                                                !> number of blocks
!      integer, dimension(:), allocatable :: blockSize                                !> size of each block (in number of vectors)
!      integer, dimension(:), allocatable :: offsets                                  !> the offsets (in number of vectors) of each block
!  end type 
  
  contains
  
  subroutine copy_cholesky_data(src, dest)
      implicit none
      type(cholesky_data), intent(in) :: src
      type(cholesky_data), intent(inout) :: dest
      integer :: n, m, ierr
      call deallocate_cholesky_data(dest)
      n = size(src%mo_ind_inv)
      allocate(dest%mo_ind_inv(n), stat=ierr)
      call allocatecheck(ierr, "dest%mo_ind_inv")
      n = size(src%cho_norms)
      allocate(dest%cho_norms(n), stat=ierr)
      call allocatecheck(ierr, "dest%cho_norms")
      n = size(src%cho_vectors,1)
      m = size(src%cho_vectors,2)
      allocate(dest%cho_vectors(n,m), stat=ierr)
      call allocatecheck(ierr, "dest%cho_vectors")     
      dest%mo_ind_inv = src%mo_ind_inv
      dest%cho_norms = src%cho_norms
      dest%cho_vectors = src%cho_vectors
  end subroutine
  
  subroutine deallocate_cholesky_data(d)
      type(cholesky_data) :: d
      call try_deallocate_int_1D(d%mo_ind_inv, "d%mo_ind_inv")
      call try_deallocate_real_1D(d%cho_norms, "d%cho_norms")
      call try_deallocate_real_2D(d%cho_vectors,"d%cho_vectors")
  end subroutine

    
end module cholesky_structs
