! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module davidson_restart
    use global_var_mod
    use davidson_io_mod
    use math_utils
    implicit none

    type restarted_root
        integer :: root
        integer :: nVec
        integer, allocatable :: newSubspaceIndex(:)
        real(kind=real8), allocatable :: vectors(:,:)
    end type

    contains

    ! **********************************************************************************
    !> \brief Function determines if we need to restart the Davidson subspace 
    !>        because we are running out of disk space
    !> \param nvectors The current number of vectors in the Davidson subspace
    !> \param nroots   The number of roots being solved for
    !> \param vector_size The size of a single vector (in # of real numbers)
    function need_to_restart_subspace(nvectors, nroots, vector_size)
        implicit none
        logical :: need_to_restart_subspace
        integer, intent(in) :: nvectors, vector_size, nroots
        real(kind=real8) :: mem_required, used_space

        ! never restart if we have less than 2 vectors per root 
        ! (this avoids restarting if we don't have 2 iterations to restart with)
        if (nvectors <= 4 * nroots) then
            need_to_restart_subspace = .false.
            return
        end if

        ! total hard drive space required per vector in GB (SI)
        mem_required = (vector_size * 8. /1000. /1000. /1000.)

        ! given the assumed disk space (in GB) figure out if we need
        ! to restart. Remember we actually store 2 vectors per iteration (c,Hc)
        used_space = mem_required * nvectors * 2
        if (used_space >= ASSUMED_DAV_DISK_SPACE) then
            need_to_restart_subspace = .true.
        else
            need_to_restart_subspace = .false.
        end if

#ifdef DEBUG_TIGER
        write(*, *) "----------------------------"
        write(*, *) "Davidson Restart DEBUG      "
        write(*, *) "----------------------------"
        write(*, *) "Size of each vector (#real)= ", vector_size
        write(*, *) "Size of each vector (GB)   = ", mem_required
        write(*, *) "Size of the subspace       = ", nvectors
        write(*, *) "Used disk space (GB)       = ", used_space
        write(*, *) "Assumed disk size (GB)     = ", ASSUMED_DAV_DISK_SPACE
        write(*, *) "Am I going to restart ?      ", need_to_restart_subspace
        write(*, *) " "
#endif
    end function
    
     ! **********************************************************************************
     !> \brief Restarts the davidson subspace using the last two best guesses for the CI vector
     !> \param nroots Number of roots being solved for
     !> \param converged  A list of which roots are already converged
     !> \param subspace_size Current size of the davidson subspace
     !> \param current_wavefuncs  The current best guess for the CI vector (written in the davidson subspace), one for each root
     !> \param previous_wavefuncs  The previous best guess for the CI vector (written in the davidson subspace), one for each root
     !> \param scratch1 Scratch space (size=total_csfs)
     !> \param scratch2 Scratch space (size=total_csfs)
     !> \param H   The Hamiltonian matrix on the current davidson subspace
     !> \param G   The G matrix on the current davdison subspace (optional)
     subroutine restart_davidson_subspace(nroots, converged, subspace_size, current_wavefuncs, &
                                          previous_wavefuncs, scratch1, scratch2, H, G)
         implicit none
        
         ! inputs
         integer, intent(in) :: nroots
         real(kind=real8), intent(inout) :: current_wavefuncs(:,:), previous_wavefuncs(:,:)
         real(kind=real8), intent(inout) :: scratch1(:), scratch2(:)
         integer, intent(inout) :: subspace_size
         real(kind=real8), intent(inout) :: H(:,:)
         logical, intent(in) :: converged(:)
         real(kind=real8), intent(inout), optional :: G(:,:)
        
         ! local variables
         integer :: iroot, new_subspace_size, ierr, p
         type(restarted_root), allocatable :: rootdata(:)
        
         ! (1) set up all the info describing how we want to restart 
         !     and if necessary, orthogonalize the two wavefunctions
         allocate(rootdata(nroots), stat=ierr)
         call allocatecheck(ierr, "rootdata")
         new_subspace_size = 0 
         do iroot = 1, nroots
             call setup_root(new_subspace_size, rootdata(iroot), converged(iroot), &
                             current_wavefuncs(:,iroot), previous_wavefuncs(:,iroot), subspace_size)
         end do

         ! (2) Collapse the c and H*c vectors to the new subspace {current_wavefuncs, previous_wavefuncs}
         call collapse_vectors(rootdata, cfile, ci_final, scratch1, scratch2, subspace_size, nroots)
         call collapse_vectors(rootdata, Hcfile, ci_final, scratch1, scratch2, subspace_size, nroots)
        
         ! (3) Build H in the new subspace
         call build_new_operator(H, rootdata, nroots, subspace_size, new_subspace_size)
         if (present(G)) then 
             call build_new_operator(G, rootdata, nroots, subspace_size, new_subspace_size)
         end if
        
         ! (4) Reset the current wavefuncs to reflect the new subspace
         do iroot = 1, nroots
             current_wavefuncs(:,iroot) = 0.0
             ! the current best guess for the ith root wavefunc is always the first vector in rootdata
             p = rootdata(iroot)%newSubspaceindex(1)
             current_wavefuncs(p,iroot) = 1.0
         end do
         
#ifdef DEBUG_TIGER
         write(*,*) "new wavefuncs"
         do iroot = 1, nroots
             write(*,*) current_wavefuncs(:,iroot)
         end do
#endif
         
         ! (5) Set the new subspace size
         subspace_size = new_subspace_size  
     end subroutine  restart_davidson_subspace
     
     
     ! **********************************************************************************
     !> \brief Sets up the data describing how we want to restart a single root
     !> \param id   The current maximum id for the new subspace (i.e. (id+1) isn't occupied yet)
     !> \param d  Data for this root 
     !> \param converged True if this root is already converged
     !> \param curr_wave  Current best guess for this wavefunction
     !> \param prev_wave  Previous best guess for this wavefunction
     !> \param n  The size of the wavefunctions
     subroutine setup_root(id, d, converged, curr_wave, prev_wave, n)
         implicit none
         integer, intent(inout) :: id 
         integer, intent(in) :: n
         type(restarted_root) :: d
         logical, intent(in) :: converged
         real(kind=real8), intent(in) :: curr_wave(:), prev_wave(:)
         integer :: ierr
         
         if ( .not. converged ) then 
             ! We restart using the last two "best guesses" for the wavefunction
             d%nVec = 2 
             allocate(d%newSubspaceIndex(2), d%vectors(n,2), stat=ierr)
             call allocatecheck(ierr, "d%newSubspaceIndex(2)")
             d%newSubspaceIndex(1) = id + 1
             d%newSubspaceIndex(2) = id + 2
             id = id + 2
             d%vectors(1:n,1) = curr_wave(1:n)
             d%vectors(1:n,2) = prev_wave(1:n)
             call orthonormalize(d%vectors(:,1), d%vectors(:,2))
         else
             ! If this root already converged just use the converged result
             ! (this avoids using 2 basically linear dependent vectors which may
             ! create issues)
             d%nVec = 1
             allocate(d%newSubspaceIndex(1), d%vectors(n,1), stat=ierr)
             call allocatecheck(ierr, "d%newSubspaceIndex(1)")
             d%newSubspaceIndex(1) = id + 1
             id = id + 1
             d%vectors(1:n,1) = curr_wave(1:n)
         end if
         
#ifdef DEBUG_TIGER
         write(*,*) "Debuggng data setup for davidson subspace restart"
         write(*,*) "Number of vectors for this root = " , d%nVec
         write(*,*) "The index of this vector(s) in the new subspace is/are = " , d%newSubspaceIndex
         write(*,*) "New Vectors .... "
         do ierr =1, d%nVec
             write(*,*) d%vectors(:,ierr)
         end do   
#endif
         
     end subroutine
 
     ! **********************************************************************************
     !> \brief Orthonormalize two vectors
     subroutine orthonormalize(v1, v2)
         implicit none
         real(kind=real8), intent(inout) :: v1(:), v2(:)
         real(kind=real8) :: norm
         ! just in case normalize v1
         norm = sqrt(dot_product(v1,v1))
         v1 = v1 / norm
         ! now make v2 orthogonal to v1
         v2 = v2 - dot_product(v1,v2) * v1
         ! normalize the new v2
         norm = sqrt(dot_product(v2,v2))
         v2 = v2 / norm 
     end subroutine  
     
     ! **********************************************************************************
     !> \brief Collapses the current Davidson subspace vectors into the new (1 or 2) vector(s) per root
     !> \param d  The root data 
     !> \param dav_file The file containing the davidson vectors, on input contains the old subspace, 
     !>                 on exit contains the new subspace
     !> \param scratch_file A file used for scratch space
     !> \param v1  A scratch vector (size of the wavefunction)
     !> \param v2 A scratch vector (size of the wavefunction)
     !> \param subspace_size Size of the old davidson subspace
     !> \param nroot Number of roots being solved for
     subroutine collapse_vectors(d, dav_file, scratch_file, v1, v2, subspace_size, nroot)
          implicit none
          type(restarted_root), intent(in) :: d(:)
          character(len=*), intent(in) :: dav_file, scratch_file
          real(kind=real8), intent(inout) :: v1(:), v2(:)
          integer, intent(in) :: subspace_size, nroot    
     
          integer :: iroot, ivec, p
          
         ! (1) build the vectors in the CI basis instead of the davidson subspace and place them into the scratch file
         ! Note that we need the scratch file because we can't overwrite the dav_file until we finish building the new
         ! vectors for every root !!!
         do iroot = 1, nroot
             do ivec = 1, d(iroot)%nVec
                 p = d(iroot)%newSubspaceIndex(ivec)
                 call build_write_contracted_vector(d(iroot)%vectors(:,ivec), v1, v2, dav_file, scratch_file, p, subspace_size)
             end do
         end do
         
        ! (3) place the new vectors in the correct file
         do iroot = 1, nroot
             do ivec = 1, d(iroot)%nVec
                 p = d(iroot)%newSubspaceIndex(ivec)
                 call swap_vectors_on_disk(v1, scratch_file, p, dav_file, p)
             end do
         end do
     end subroutine 
    
     ! **********************************************************************************
     !> \brief swaps a vector from one disk location to another
     !> \param A scratch vector (size of the wavefunction)
     !> \param src_file The source file which is read from 
     !> \param src _pos  The position in the source file to read from 
     !> \param dest_file The file where the vector is copied to 
     !> \param dest_pos  The position in dest_file the vector is copied to
     subroutine swap_vectors_on_disk(scratch, src_file, src_pos, dest_file, dest_pos)
         implicit none
         real(kind=real8), intent(out) :: scratch(:)
         character(len=*), intent(in) :: src_file, dest_file
         integer, intent(in) :: src_pos, dest_pos
        
         call davidson_io_read(scratch, src_file, src_pos)
         call davidson_io_write(scratch, dest_file, dest_pos)
     end subroutine

     ! **********************************************************************************
     !> \brief Builds a symmetric operator in the new davidson subspace, given the 
     !>        representation in the old subspace
     !> \param op  The operator in the old subspace, on exit contains the new representation
     !> \param d   The root data containing all data on the restarted vectors
     !> \param nroot The number of roots
     !> \param old_subspace_size Size of the old subspace
     !> \param new_subspace_size Size of the new subspace
     subroutine build_new_operator(op, d, nroot, old_subspace_size, new_subspace_size)
         implicit none
         real(kind=real8), intent(inout) :: op(:,:)
         type(restarted_root), intent(in) :: d(:)
         integer, intent(in) :: nroot, old_subspace_size, new_subspace_size
         
         integer :: iroot, ivec, p, ierr
         real(kind=real8), allocatable ::new_op(:,:), U(:,:)
        
         ! We are going to build the new operator with the classic unitary transform 
         ! O_new = U^T * O * U
         allocate(U(old_subspace_size,new_subspace_size), new_op(new_subspace_size, new_subspace_size), stat=ierr)
         call allocatecheck(ierr, "temporaries in build_new_operator")  
         do iroot = 1, nroot
             do ivec = 1, d(iroot)%nVec
                 p = d(iroot)%newSubspaceIndex(ivec)
                 U(:,p) = d(iroot)%vectors(:,ivec)
             end do
         end do
         
         ! this should be small enough that the temporaries it generates won't harm us
         new_op = matmul(transpose(U), matmul(op(1:old_subspace_size,1:old_subspace_size), U)) 
#ifdef DEBUG_TIGER
         write(*,*) "Debug the operator transform ...."
         call error_matrix_write(U, "transform U")
         call error_matrix_write(op(1:old_subspace_size,1:old_subspace_size), "old operator")
         call error_matrix_write(new_op(1:new_subspace_size,1:new_subspace_size), "new operator")
#endif
         op(1:new_subspace_size,1:new_subspace_size) = new_op(1:new_subspace_size,1:new_subspace_size)      
     end subroutine
     
     
     ! **********************************************************************************
     !> \brief We encounter multiple times vectors written as \f$v=\sum_ic_i\phi_i\f$ where we have the coefficients and
     !>        the \f$\phi\f$ are on disk. This routines builds the vector and writes it to a given file and location
     !> \param coeff The vector (written as a linear combination of some set of vectors)
     !> \param vector Space for the final vector (size of the wavefunction)
     !> \param scratch Scratch space (size of the wavefunction)
     !> \param basis_file  A file containing the set of vectors (davidson subspace) we wrote coeff with
     !> \param dest_file   A file where the final vector is written to
     !> \param dest_file_pos  The positon in the dest_file where we write
     !> \param n The number of vectors in basis_file and therefore the size of coeff
     subroutine build_write_contracted_vector(coeff, vector, scratch, basis_file, dest_file, dest_file_pos, n)
         implicit none
         ! inputs
         real(kind=real8), intent(in) :: coeff(:)
         real(kind=real8), intent(out) :: vector(:), scratch(:)
         character(len=*), intent(in) :: basis_file, dest_file
         integer, intent(in) :: dest_file_pos, n
        
         ! local variable
         integer :: i 
        
         vector = 0.0
         do i = 1, n 
             call davidson_io_read(scratch, basis_file, i)
             vector = vector + coeff(i) * scratch 
         end do
         call davidson_io_write(vector, dest_file, dest_file_pos)
     end subroutine

end module
