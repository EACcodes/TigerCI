! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module number_virt_per_hole
  use global_var_mod
  use tree_search_mod
  use graph_var_mod
  use allowed_virtuals_mod
  use new_tree_search_mod
  use new_tree_search_structs
  
  implicit none
  contains

    subroutine count_virts_per_hole()
      type(orbital_path) :: path 
      type(graph_search_state) :: G
      
      integer :: n_holes 
      integer :: n_virts
      integer :: elec, virts, weight
      integer :: path_address_ab, path_address_aa

      ! write out each step to disk
      !open (unit = 752, file = "particles-holes.txt", form = "formatted")

      ! allocate the path
      call allocate_orbital_path(path)
      call zero_orbital_path(path)

      ! Run a search through the internal configurations. For this analysis
      ! we are going to focus just on states w/ 2 virtuals
      elec =num_elec - 2
      n_holes = 0 
      n_virts = 0 
      call init_tree_search(G, 0, 0, num_internal, elec)
      do while(get_internal_path(path,G))

         weight = sum(path%arc_weights(0:num_internal)) + 1

         ! don't include this configuration if it doesn't exist
         path_address_ab = internal_index_vector2(weight)
         path_address_aa = internal_index_vector3(weight)
         if ( (path_address_ab<0) .or. (path_address_aa<0)) cycle

         ! record info about the number of virtuals we can excite to for this hole pair
         virts = num_allowed_virtuals(weight, "D")
         !write(752,*) weight, virts
         n_holes = n_holes + 1
         n_virts = n_virts + virts
      end do

      write(*,*) " " 
      write(*,*) "Finished # virtual orbitals / ij hole analysis"
      write(*,*) "Total number of holes = " , n_holes
      write(*,*) "Average number of virtual orbitals per hole (for just the D vertex!)= " , n_virts/real(n_holes,kind=real8)
      write(*,*) ""
      !close(unit=752)
         
    end subroutine count_virts_per_hole


end module number_virt_per_hole
