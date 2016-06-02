! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module new_tree_search_structs
    use global_var_mod
    use utilities_mod
    use molecule_var_mod
    implicit none

    !// HERE IS OUR DERIVED DATA TYPE.  
    type orbital_path

        !// ARRAYS
        integer, dimension(:), allocatable :: occupations !// STORES OCCUPATION PATTERN
        integer, dimension(:), allocatable :: arc_weights !// WEIGHT OF EACH ARC IN PATTERN
        integer, dimension(:), allocatable :: singles !// ARRAY TO REORDER SINGLES 
        integer, dimension(:,:), allocatable :: constraints !// HOLDS LOOP INFORMATION

        !// SCALARS
        integer :: num_singles !// TOTAL NUMBER OF SINGLES IN PATH
        integer :: weight !// TOTAL WEIGHT OF PATH

        integer :: level1 !// USED TO HOLD INFO ABOUT LOOP LEVELS
        integer :: level2
        integer :: level3
        integer :: level4

        integer :: i_segment !// ALSO STORES INFORMATION ABOUT THE LEVELS
        integer :: j_segment
        integer :: k_segment
        integer :: l_segment

        integer :: loop_type !// KEEPS TRACK OF THE LOOP TYPE WE ARE MAKING
        integer :: rt_loop_weight !// USED TO HOLD INFORMATION ABOUT THE "OTHER" HALF OF THE LOOP
        integer :: start_elec !// NUMBER OF ELECTRONS AT HEAD OF LOOP

        real(real8) :: element1 !// THESE REALS ARE USED TO HOLD PARTIAL MATRIX ELEMENTS
        real(real8) :: element2 !// OR SOME INTEGRALS
        real(real8) :: element3

        real(real8), dimension(:), allocatable :: integrals !// HOLDS INTEGRAL COMBOS.

        logical, dimension(:), allocatable :: encountered !// THIS ARRAY KEEPS TRACK OF WHETHER
        !// NOT WE HAVE ENCOUNTERED ONE OF THE 
        !// CONSTRAINED LEVELS.    

    end type orbital_path


    type graph_search_state
        integer :: level !// THE LEVEL WE ARE ON
        integer :: current_vertex !// THE VERTEX WE ARE ON
        integer :: ending_vertex !// THE VERTEX WE END ON
        integer :: step_type !// THE TYPE OF STEP WE ARE TAKING
        integer :: path_elecs !// THE NUMBER OF ELECTRONS IN THE PATH
        integer :: top_level !// TOP LEVEL OF THE PATH
        integer :: bottom_level !// BOTTOM LEVEL OF PATH
        integer :: top_occ !// NUMBER OF ELECTRONS AT TOP LEVEL
        integer :: bottom_occ !// NUMBER OF ELECTRONS AT BOTTOM LEVEL      
    end type graph_search_state

contains

  ! ALL OF THESE FUNCTIONS HANDLE JUST THE PARTS OF THE ORBITAL PATH
  ! WHICH ARE ACTUALLY THE 'PATH' 
  ! 
  ! I'M IGNORING THE TACKED ON PARTS  (CONSTRAINTS/INTEGRALS)
  !
  !> \todo Please fix the orbital_path data type (see above)


  logical function is_orbital_path_allocated(path)
    implicit none
    type(orbital_path) :: path 
    if (allocated(path%occupations) .and. allocated(path%singles) .and. allocated(path%arc_weights)) then
       is_orbital_path_allocated= .true.
    else
       is_orbital_path_allocated = .false.
    end if
  end function is_orbital_path_allocated

  subroutine allocate_orbital_path(path)
    implicit none
    ! inputs 
    type(orbital_path) :: path

    ! local variables
    integer :: ierr

    if ( .not. allocated(path%arc_weights) ) then 
       allocate(path%arc_weights(0:num_orbitals), stat=ierr)
       call allocatecheck(ierr, "path%arc_weights")
    end if
    if ( .not. allocated(path%occupations) ) then 
       allocate(path%occupations(0:num_orbitals), stat=ierr)
       call allocatecheck(ierr, "path%occupations")
    end if
    if ( .not. allocated(path%singles) ) then 
       allocate(path%singles(0:num_orbitals), stat=ierr)
       call allocatecheck(ierr, "path%singles")
    end if
    call zero_orbital_path(path)
  end subroutine allocate_orbital_path

  subroutine deallocate_orbital_path(path)
    ! input
    type(orbital_path) :: path
    ! local variable
    integer :: ierr
    deallocate(path%arc_weights, path%occupations, path%singles, stat=ierr)
    call deallocatecheck(ierr, "a path")
  end subroutine deallocate_orbital_path


  subroutine zero_orbital_path(path)
    implicit none
    ! input
    type(orbital_path) :: path

    ! scalars
    path%num_singles = 0
    path%weight = 0

    path%level1 = 0
    path%level2 = 0
    path%level3 = 0
    path%level4 = 0

    path%i_segment = 0 
    path%j_segment = 0
    path%k_segment = 0
    path%l_segment = 0

    path%loop_type = 0
    path%rt_loop_weight = 0 
    path%start_elec = 0 

    path%element1 = 0
    path%element2 = 0
    path%element3 = 0

    ! vectors
    path%occupations = 0 
    path%arc_weights = 0 
    path%singles = 0 
  end subroutine zero_orbital_path


  logical function are_paths_equivalent(path1, path2)
    implicit none
    ! inputs
    type(orbital_path) :: path1, path2

    if (path1%num_singles/=path2%num_singles) then
       write(*,*) "PATHS ARE NOT EQUIVALENT: num_singles"
       are_paths_equivalent = .false.
    else if (.not. arraysAreEqual(path1%singles, path2%singles) ) then
       write(*,*)"PATHS ARE NOT EQUIVALENT: singles"
       are_paths_equivalent = .false.
    else if (.not. arraysAreEqual(path1%arc_weights, path2%arc_weights)) then
       write(*,*)"PATHS ARE NOT EQUIVALENT: arc_weights"
       write(*,*) path1%arc_weights
       write(*,*) path2%arc_weights
       are_paths_equivalent = .false.
    else if (.not. arraysAreEqual(path1%occupations,path2%occupations))  then
       write(*,*)"PATHS ARE NOT EQUIVALENT: occupations"
       are_paths_equivalent = .false.
    else
       are_paths_equivalent = .true.
    end if

  end function are_paths_equivalent

  subroutine write_path(the_path, astring)

    !// THIS SUBROUTINE WRITES OUT ALL THE RELEVANT INFORMATION ON
    !// THE PATH.  IT IS MAINLY FOR DEBUGGING


    implicit none

    type(orbital_path)::the_path
    character(len=*)::astring
    integer::i

    write(*,*) "*******************************************"

    write(*,*) "Path info: ", astring
    write(*,*)

    write(*,*) "Levels:"
    write(*,*) "level 1:", the_path%level1
    write(*,*) "level 2:", the_path%level2
    write(*,*) "level 3:", the_path%level3
    write(*,*) "level 4:", the_path%level4
    write(*,*)

    write(*,*) "Segments: "
    write(*,*) "i_segment: ",the_path%i_segment
    write(*,*) "j_segment: ",the_path%j_segment
    write(*,*) "k_segment: ",the_path%k_segment
    write(*,*) "l_segment: ",the_path%l_segment
    write(*,*)

    write(*,*) "Reals: "
    write(*,*) "element1: ",the_path%element1
    write(*,*) "element2: ",the_path%element2
    write(*,*) "element3: ",the_path%element3

    write(*,*) "Occupations:"
    write(*,500) (the_path%occupations(i), i = 1, num_orbitals)

    write(*,*) "Arc weights:"
    write(*,500) (the_path%arc_weights(i), i = 1, num_orbitals)

    write(*,*) "Singles: "
    write(*,500) (the_path%singles(i), i = 1, num_orbitals)

    write(*,*) "Total number of singles: ", the_path%num_singles
    write(*,*) "Total weight: ", the_path%weight
    write(*,*) "the_path%rt_loop_weight: ", the_path%rt_loop_weight
    write(*,*) "the_path%loop_type: ", the_path%loop_type
    write(*,*)
    write(*,*) "Constraints for loop: ", the_path%loop_type
    if ( allocated(the_path%constraints) ) then 
       write(*,500) (the_path%constraints(1,i), i = 1, size(the_path%constraints,2))
       write(*,500) (the_path%constraints(2,i), i = 1, size(the_path%constraints,2))
       write(*,500) (the_path%constraints(3,i), i = 1, size(the_path%constraints,2))
    end if
    call flush(6)

500 format(1x, 13i4, /)

  end subroutine write_path
end module new_tree_search_structs
