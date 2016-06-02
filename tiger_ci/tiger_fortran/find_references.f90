! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module find_references_mod
  use global_var_mod
  use utilities_mod
  use new_tree_search_mod
  use new_tree_search_structs
  use stack_mod
  use sort_utils
  use c_sort_finterface
  implicit none

  contains

    !> \brief Finds all the CSFs with reference occupations
    !> \param positions An array on exit containing the positions of all the reference CSFSs. Must be allocatable (we'll deallocate and reallocate it
    !>                  during the routine to the correct size). Final array is returned in ascending order
    !> \param nRefCSFs  Number of CSFs corresponding to reference occupations (depending on spin functions this isn't always equal to num_ref)
    subroutine find_all_reference_CSFs(positions, nRefCSFs)
      implicit none

      ! inputs
      integer, intent(out) :: nRefCSFs
      integer, intent(out), allocatable :: positions(:)

      ! local variables 
      integer :: weight, address, spin_dim, spin, pos
      integer :: i, ierr
      type(orbital_path) :: path 
      type(IntegerStackType) :: ref_address_stack
      type(graph_search_state) :: G

      ! Search all the paths in the graph with all N electrons in the internal space
      call allocate_orbital_path(path)
      call init_tree_search(G, 0, 0, num_internal, num_elec)
      do while( get_internal_path(path,G) )
         ! check to see if these orbital occupations match a reference
         if (is_reference(path%occupations(1:num_internal))) then 
            weight = sum(path%arc_weights)+1
            address = internal_index_vector0(weight)
            spin_dim = fsn(v_singles(weight))
            ! loop over each spin function associated with this orbital configuration
            do spin = 1, spin_dim
               pos = address + (spin-1)
               ! we found a reference CSF and its address, add it to my address_stack
               call push_stack(ref_address_stack, pos)
            end do
         end if
      end do

      ! all my addresses are on the stack, put them back into positions
      nRefCSFs = get_stack_size(ref_address_stack)
      if (allocated(positions)) then 
         deallocate(positions, stat=ierr)
         call deallocatecheck(ierr, "positions in find_all_reference_CSFs")
      end if
      allocate(positions(nRefCSFs), stat=ierr)
      call allocatecheck(ierr, "positions in find_all_reference_CSFs")
      do i = 1, nRefCSFs
         address = pop_stack(ref_address_stack)
         positions(i) = address
      end do

      ! return the final list sorted 
      !call insertion_sort_int(positions)
      call sort_int_array(positions)
    end subroutine find_all_reference_CSFs
         
end module find_references_mod
