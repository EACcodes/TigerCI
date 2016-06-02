! Copyright (c) 1998-2004, University of California, Los Angeles, Emily A. Carter
!                       2004-2016, Princeton University, Emily A. Carter
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of
!     conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!     of conditions and the following disclaimer in the documentation and/or other
!     materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!     used to endorse or promote products derived from this software without specific
!     prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
! CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
! OF THE POSSIBILITY OF SUCH DAMAGE.
!*****************************************************************
!//  LOCIST_MOD: THIS MODULE CONTAINS THE ROUTINES FOR DOING THE LOCAL
!//  CI CALCULATION.  THIS MODULE ALSO TAKES CARE OF THE VALENCE
!//  CI IF ONE SHOULD ASK FOR IT.
!//
!//  WRITTEN BY DEREK WALTER, 2000
!//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS.
!//  ERASED ALL MODULE VARIABLES: JOHANNES M. DIETERICH (2012)
!*****************************************************************

module locist_mod

  use global_var_mod
  use molecule_var_mod
  use utilities_mod
  use graph_var_mod
  use allowed_virtuals_mod
  use wp_tov_mod
  use new_tree_search_mod
  use locist_var_mod,only:locist_scratch

#ifdef TIGER_USE_OMP
  use omp_lib
#endif

  implicit none
  
  public
  
  contains

  !*****************************************************************
  !> \brief This controls the local truncation (WP/TOV). Here is where we actually throw things out
  subroutine local_driver(loc_scr)
    implicit none

    type(locist_scratch)::loc_scr

    ! Produce some output
    write(*,*) "Parameters for truncation:              "
    write(*,*) "------------------------------------------------"
    write(*,*) "Reference CI ............................",  (reference_ci_flag == 1) 
    write(*,*) "Valence CI   ............................",  (valence_ci_flag == 1)
    write(*,*)
    write(*,*) "Occupation Threshold for Internal Orbs...",internal_threshold 
    write(*,*) "Default Radius       for WP..............",wp_default_radius 
    write(*,*) "Radius  Multiplier   for WP..............",wp_multiplier 
    write(*,*) "Cylinder Radius      for WP..............",wp_cylinder_radius 
    write(*,*)
    write(*,*) "Default Radius       for TOV (Occupied) .",tov_occupied_default_radius 
    write(*,*) "Radius  Multiplier   for TOV (Occupied) .",tov_occupied_multiplier 
    write(*,*)
    write(*,*) "Occupation Threshold for Virtual Orbs....",virtual_threshold 
    write(*,*) "Default Radius       for TOV (Virtual) ..",tov_virtual_default_radius 
    write(*,*) "Radius  Multiplier   for TOV (Virtual) ..",tov_virtual_multiplier 
    write(*,*) "Cylinder Radius      for TOV.............",tov_cylinder_radius 
    write(*,*)

    ! Truncation in the internal space (WP approximation or reference/valence CI)
    write(*,*) "Internal truncation:              "
    write(*,*) "------------------------------------------------"
    call internal_truncation

    ! Do the Truncation of Virtuals (TOV) approximation 
    write(*,*) "Virtual truncation:              "
    write(*,*) "------------------------------------------------"
    call virtual_truncation(loc_scr)
  
    
  end subroutine local_driver


  !*****************************************************************
  !> \brief This routine removes configurations based on just the internal space
  !>
  !> This is where we handle the WP truncation and a few other special cases (reference CI/valence CI)
  subroutine internal_truncation
    implicit none

    ! local variables
    integer :: elec
    procedure(reference_ci_throw_out), pointer :: truncation_routine
    type(orbital_path) :: path 
    type(graph_search_state) :: internal_graph
    integer :: ierr
    integer :: num_out

    ! allocate the path variable
    allocate(path%arc_weights(0:num_orbitals), path%occupations(0:num_orbitals), path%singles(0:num_orbitals), stat=ierr)
    call allocatecheck(ierr, "a path in internal_truncation")

    ! set the truncation routine we want to run. this depends on the type of calculation we are running
    if ( reference_ci_flag == 1 ) then 
       truncation_routine => reference_ci_throw_out
    elseif ( valence_ci_flag == 1 ) then 
       truncation_routine => valence_ci_throw_out
    else 
       truncation_routine => WP_throw_out
    end if

    ! search for all the configurations in the internal graph with num_elec in the internal space
    path%occupations = 0
    path%arc_weights = 0
    path%singles = 0
    path%num_singles = 0
    elec = num_elec
    call init_tree_search(internal_graph, 0, 0, num_internal, elec)
    do while ( get_internal_path(path, internal_graph)) 
       call truncation_routine(path, internal_index_vector0)
    end do
    num_out = count( internal_index_vector0 < 0)
    write(ioOutput,140) num_out,real(num_out,real8)/real(size(internal_index_vector0),real8)

    ! search for all the configurations in the internal graph with num_elec-1 in the internal space
    path%occupations = 0
    path%arc_weights = 0
    path%singles = 0
    path%num_singles = 0
    elec = num_elec-1
    call init_tree_search(internal_graph, 0, 0, num_internal, elec)
    do while ( get_internal_path(path, internal_graph)) 
       call truncation_routine(path, internal_index_vector1)
    end do
    num_out = count( internal_index_vector1 < 0 )
    write(ioOutput,130) num_out,real(num_out,real8)/real(size(internal_index_vector1),real8)

    ! search for all the configurations in the internal graph with num_elec-2 in the internal space
    ! recall we break this into 2 parts: configurations with 2 identical virtuals are 2 different virtuals
    path%occupations = 0
    path%arc_weights = 0
    path%singles = 0
    path%num_singles = 0
    elec = num_elec-2
    call init_tree_search(internal_graph, 0, 0, num_internal, elec)
    do while ( get_internal_path(path, internal_graph)) 
       call truncation_routine(path, internal_index_vector2)
    end do
    num_out = count( internal_index_vector2 < 0 )
    write(ioOutput,120) num_out,real(num_out,real8)/real(size(internal_index_vector2),real8)
    
    path%occupations = 0
    path%arc_weights = 0
    path%singles = 0
    path%num_singles = 0
    elec = num_elec-2
    call init_tree_search(internal_graph, 0, 0, num_internal, elec)
    do while ( get_internal_path(path, internal_graph)) 
       call truncation_routine(path, internal_index_vector3)
    end do
    num_out = count( internal_index_vector2 < 0 )
    write(ioOutput,110) num_out,real(num_out,real8)/real(size(internal_index_vector3),real8) ! works since  the index2 and index3 are the same



110 format(1x,"Number of N-2 double paths thrown out by WP....",i10,/,&
         1x,"% thrown out by WP.............................",f10.3,/)

120 format(1x,"Number of N-2 single paths thrown out by WP....",i10,/,&
         1x,"% thrown out by WP.............................",f10.3,/)      

130 format(1x,"Number of N-1 paths thrown out by WP...........",i10,/,&
         1x,"% thrown out by WP.............................",f10.3,/)      

140 format(1x,"Number of V paths thrown out by WP.............",i10,/,&
         1x,"% thrown out by WP.............................",f10.3,/)

  end subroutine internal_truncation

  !*****************************************************************
  !> \brief Checks if we want to remove a path based on the Weak Pairs approximation
  !> \param path  The path (i.e. configuration) we are considering 
  !> \param index A path index (if we remove the path the index entry for path is set to < 0 )
  subroutine WP_throw_out(path, index)
    implicit none

    ! input
    type(orbital_path), intent(in) :: path 
    integer, dimension(:) :: index

    ! local variables
    integer :: i
    integer :: num_pair, hole_pair_list(num_ref, 2)
    integer :: weight
    logical :: double_exc_only, truncate

    ! check to see if this is a reference (we don't truncate those)
    if (is_reference(path%occupations(1:num_internal)) ) return

    ! Find the (i,j) hole in the path (possibly more than 1 in an MR calculation)
    weight = sum(path%arc_weights(0:num_internal)) + 1
    if ( index(weight) == kill_number ) return                 ! don't bother going on if we already removed it 
    call find_hole_pairs(path, double_exc_only, num_pair, hole_pair_list)
    
    ! If this path could be construction via a single excitation from a reference, stop
    ! (we don't truncate single excitations)
    if ( .not. double_exc_only ) return 
    
    ! If the hole pair (i,j) does not interact under WP remove this path
    ! (In an MR case: if ALL of the hole pairs don't interact then remove the path. Otherwise keep it)
    truncate = .true.
    do i = 1, num_pair
        if ( are_interacting_WP(hole_pair_list(i,1),  hole_pair_list(i,2)) ) then
            truncate = .false.
        end if
    end do
    if (truncate) then
        index(weight) = kill_number
    end if
  end subroutine WP_throw_out

  !*****************************************************************
  !> \brief Checks if we want to remove a path because it isn't a valence configuration
  !> \param path  The path (i.e. configuration) we are considering 
  !> \param index A path index (if we remove the path the index entry for path is set to < 0 )
  subroutine valence_ci_throw_out(path, index)
    implicit none

    ! input
    type(orbital_path), intent(in) :: path 
    integer, dimension(:) :: index

    ! local variables
    integer :: elec, weight

    elec = sum(path%occupations(0:num_internal))
    weight = sum(path%arc_weights(0:num_internal)) + 1

    if ( elec /= num_elec ) then 
       ! this path contains an excitation into the virtual space ... not a valence configuration
       index(weight) = kill_number
    end if

  end subroutine valence_ci_throw_out

  !*****************************************************************
  !> \brief Checks if we want to remove a path because it isn't a reference
  !> \param path  The path (i.e. configuration) we are considering 
  !> \param index A path index (if we remove the path the index entry for path is set to < 0 )
  subroutine reference_ci_throw_out(path, index)
    implicit none

    ! input
    type(orbital_path), intent(in) :: path 
    integer, dimension(:) :: index

    ! local variables
    logical :: reference_check
    integer :: weight

    weight = sum(path%arc_weights(0:num_internal)) + 1
    reference_check = is_reference(path%occupations(1:num_internal))

    if ( .not. reference_check ) then 
       ! this is not a reference ... remove it
       index(weight) = kill_number
    end if

  end subroutine reference_ci_throw_out

  !*****************************************************************
  !> \brief Checks if either given internal configuration has been removed from our calculation
  !> \param weight1  The internal weight of the first internal configuration
  !> \param weight1  The internal weight of the second internal configuration
  !> \param elec1 The number of electrons in the first internal configuration
  !> \param elec2 The number of electrons in the first internal configuration
  function skip_these_internals(weight1,weight2,elec1,elec2)

    !// TESTS TO SEE IF TWO SPHERES OVERLAP

    logical::skip_these_internals
    integer, intent(in)::weight1, weight2, elec1, elec2

    skip_these_internals = .false.

    if ( skip_this_internal(weight1,elec1) .or. skip_this_internal(weight2,elec2)) then
       skip_these_internals = .true.
    end if

  end function skip_these_internals

  !*****************************************************************
  !> \brief Checks if an internal configuration has been removed from our calculation
  !> \param weight  The internal weight of the internal configuration
  !> \param elec The number of electrons in the internal configuration
  function skip_this_internal(weight,elec)

    !// TESTS TO SEE IF TWO SPHERES OVERLAP
    logical::skip_this_internal
    integer, intent(in)::weight,  elec

    skip_this_internal = .false.

    if (elec == num_elec) then

       if (internal_index_vector0(weight + 1) < 0) then
          skip_this_internal = .true.
          return
       endif

    elseif (elec == num_elec - 1) then

       if (internal_index_vector1(weight + 1) < 0) then
          skip_this_internal = .true.

          return
       endif

    elseif (elec == num_elec - 2) then

       if (internal_index_vector2(weight + 1) < 0) then
          skip_this_internal = .true.
          return
       endif

    endif

  end function skip_this_internal


  !*****************************************************************
  !> \brief Sets up the truncation of virtuals (TOV) approximation
  !> \param loc_scr   Scratch space
  !>
  !> We accomplish the TOV approximation by determining which virtuals each
  !> internal configuration can excite to. We then record these and later on
  !> during the matrix vector product we can check the allowed virtuals lists.
  subroutine virtual_truncation(loc_scr)
    implicit none

    ! input
    type(locist_scratch)::loc_scr

    ! local variables
    integer :: ierr      
    integer :: elec, weight, kept
    integer :: total, removed
    integer, dimension(:), allocatable :: virtuals

    type(orbital_path)::path
    type(graph_search_state) :: graph

    ! Allocate a path variable
    allocate(path%arc_weights(0:num_orbitals), path%occupations(0:num_orbitals), path%singles(0:num_orbitals), stat=ierr)
    call allocatecheck(ierr,"%lam_sin")

    ! Allocate space for recording virtual orbitals
    allocate(virtuals(num_external), stat=ierr)
    call allocatecheck(ierr, "array for virtuals")
    
    ! Allocate space for storing the allowed excitations
    call init_allowed_virtuals()

    ! Check all the internal paths with N-1 electrons
    elec = num_elec - 1
    total = 0 
    path%occupations = 0
    path%singles = 0
    path%arc_weights = 0
    path%num_singles = 0
    call init_tree_search(graph, 0, 0, num_internal, elec)
    do while (get_internal_path(path, graph))
        ! Find the allowed excitations and record them
        weight = sum(path%arc_weights(0:num_internal)) + 1
        if (internal_index_vector1(weight) < 0) cycle 
        call find_allowed_virtuals(path, loc_scr, virtuals, removed)
        total = total + removed
    end do
    write(*,*) "Number of (N-1) excitations removed              ", total
    
    ! Check all the internal paths with N-2 electrons
    ! First double excitations to one orbital and then double excitations to two different virtuals
    elec = num_elec - 2
    total = 0 
    path%occupations = 0
    path%singles = 0
    path%arc_weights = 0
    path%num_singles = 0
    call init_tree_search(graph, 0, 0, num_internal, elec)
    do while (get_internal_path(path, graph))
        ! Find the allowed excitations and record them
        weight = sum(path%arc_weights(0:num_internal)) + 1
        if (internal_index_vector2(weight) < 0) cycle
        call find_allowed_virtuals(path, loc_scr, virtuals, removed)
        total = total + removed                                         ! double excitations to 1 orbital
    end do
    write(*,*) "Number of (N-2) excitations removed to 1 virtual " , total
    
    elec = num_elec - 2
    total = 0 
    path%occupations = 0
    path%singles = 0
    path%arc_weights = 0
    path%num_singles = 0
    call init_tree_search(graph, 0, 0, num_internal, elec)
    do while (get_internal_path(path, graph))
        ! Find the allowed excitations and record them
        weight = sum(path%arc_weights(0:num_internal)) + 1
        if (internal_index_vector3(weight) < 0) cycle
        call find_allowed_virtuals(path, loc_scr, virtuals, removed)
        kept = num_external - removed
        total = total + ( num_external*(num_external+1)/2 - kept*(kept+1)/2 )      ! double excitations to 2 orbitals 
    end do
    write(*,*) "Number of (N-2) excitations removed to 2 virtuals = " , total
    write(*,*)

    ! Clean up
    deallocate(path%arc_weights,path%occupations, path%singles, virtuals, stat=ierr)
    call deallocatecheck(ierr, "lambda path in virtual_truncation ")

  end subroutine virtual_truncation


  !> \brief For a given path finds the allowed virtual excitations under the TOV approximations
  !> \param path  The internal configuration
  !> \param loc_scr Scratch space
  !> \param virtuals Scratch space for storing virtual orbital indicies
  !> \param removed  On exit the number of virtual orbitals removed (NOT CONFIGURATIONS!)
  !> 
  !> This routine finds the holes (i,j) where we excite from (or multiple holes from different references in the MR case). 
  !> We then determine which virtuals are close enough to interact under the TOV approximation, (i,j) -> a is allowed. We 
  !> then store that list for future use. 
  subroutine find_allowed_virtuals(path, loc_scr, virtuals, removed)
    implicit none

    ! inputs
    type(orbital_path):: path
    type(locist_scratch)::loc_scr
    integer, dimension(:) :: virtuals
    integer, intent(out) :: removed

    ! local variables
    integer :: weight, elec
    integer :: num_pair, hole_pair_list(num_ref,2)
    integer :: iPair, i, j, a, b
    integer :: num_virt
    logical :: only_double_exc, add
    
    ! step (1) some preliminaries ...
    loc_scr%number_paths = loc_scr%number_paths + 1
    weight = sum(path%arc_weights(0:num_internal)) + 1
    elec   = sum(path%occupations(1:num_internal))
    removed = 0 
    virtuals = 0

    ! step (2) find the (ij) holes we excite from 
    call find_hole_pairs(path, only_double_exc, num_pair, hole_pair_list)

    ! step (3) if this isn't a double excitation then its a single excitation
    ! (since we never call this routine with a reference). Set the list of 
    ! allowed virtuals to be all virtuals (we don't truncate single excitations
    if ( .not. only_double_exc ) then
        do a = num_internal + 1, num_orbitals
            virtuals(a-num_internal) = a
        end do
        num_virt = num_external
    else
        ! step (4) Its a double excitation. Find the allowed virtual orbital for this/these hole pair(s)
        num_virt = 0 
        do iPair = 1, num_pair
           i = hole_pair_list(iPair,1)
           j = hole_pair_list(iPair,2)
           if ( .not. are_interacting_WP(i,j) ) cycle
           ! for the hole pair (i,j) find the virtual orbitals which interact under TOV
           do a = num_internal+1, num_orbitals
              if ( are_interacting_TOV(i,j,a) ) then 
                  ! check to see if this virtual was already added
                  add = .true.
                  do b = 1, num_virt
                      if ( virtuals(b) == a ) then 
                          add = .false.
                      end if
                  end do
                  if ( add ) then
                      num_virt = num_virt + 1
                      virtuals(num_virt) = a 
                  end if
              end if
           end do
        end do
    end if

    ! step (5) add this weight and set of allowed excitations 
    call add_new_virtual_list(weight, elec, num_virt, virtuals)
    
    ! step (6) record how many virtuals (if any) were removed
    removed = num_external - num_virt
  end subroutine find_allowed_virtuals
  
  
  
  !*****************************************************************
  !> \brief Given an internal path this routine determines all the possible hole pairs
  !> \param path                The internal configuration
  !> \param only_double         On exit is true ONLY IF this configuration can only be formed by a double excitation
  !> \param num_pair            On exit contains the number of hole pairs found
  !> \param hole_pair_list      On exit contains all the possible hole pairs 
  !>
  !> For a single reference calculation determining where the holes are in a path is relatively easy. For a 
  !> multireference problem this isn't as easy. Its possible to find scenarios where a single configuration
  !> could have arisen from different references - meaning that we could associate multiple hole pairs to 
  !> one configuration. This routine returns a list of all possible hole pairs for the given configuration.
  subroutine find_hole_pairs(path, only_double, num_pair, hole_pair_list)
    implicit none

    ! inputs
    type(orbital_path), intent(in) :: path 
    logical, intent(out) :: only_double 
    integer, intent(out) :: num_pair
    integer, dimension(:,:) :: hole_pair_list

    ! local variables
    integer :: i, j
    integer :: num_holes, num_particles
    integer :: array(1:num_internal)

    num_pair = 0 
    hole_pair_list = 0 
    only_double = .true.

    do i = 1, num_ref
       ! this is a reference exit immediately ... I will never care about holes here
       if ( is_reference(path%occupations(1:num_internal))) then 
          only_double = .false.
          num_pair = 0 
          return 
       end if

       array = references(1:num_internal,i) - path%occupations(1:num_internal) 
       num_holes = sum(array,array>0)                   ! All the positive entries are holes
       num_particles = abs(sum(array, array<0))         ! All the negative entries are particles

       if ( num_holes > 2 ) cycle            ! No triple or higher excitations
       if ( num_particles > 2 ) cycle        ! No triple or higher excitations, plus one electron should be in the virtual space

       ! since we are using the hole pair info for local truncation and because we never 
       ! truncate single excitations, if I can represent this configuration as a single 
       ! excitation from a reference -> stop looking for holes and just exit
       if ( num_holes == 1 ) then 
          only_double = .false.
          num_pair = 0 
          return 
       end if

       ! At this point I have already exited if I had 0 holes (reference) or 1 hole (single excitation)
       ! this is a double excitation. Lets look for where the holes are 
       num_pair = num_pair + 1
       do j = 1, num_internal
          ! 2 holes in one orbital
          if ( array(j) == 2 ) then
             hole_pair_list(num_pair,1) = j
             hole_pair_list(num_pair,2) = j
          end if

          ! 1 hole in the orbital
          if ( array(j) == 1 ) then
             ! if I haven't found a hole yet
             if ( hole_pair_list(num_pair,1) == 0 ) then
                hole_pair_list(num_pair,1) = j 
             else ! otherwise fill in the second entry
                hole_pair_list(num_pair,2) =j
             end if
          end if
       end do ! loop over internals
    end do ! loop over references

    ! Sanity check: I never exited because this config was a reference or single excitation
    ! therefore I should have found at least 1 hole pair (there are no triple or higher excitations!)
    if ( num_pair == 0 ) then 
       write(*,*) "Fatal error in find_hole_pairs"
       write(*,*) "For some reason I think this configuration isn't"
       write(*,*) " a reference, single or double excitation ????"
       write(*,*) path%occupations(1:num_internal)
       stop
    end if

  end subroutine find_hole_pairs

end module locist_mod

