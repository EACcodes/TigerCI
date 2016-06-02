! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief This module contains the spheres and cylinders we use for the WP/TOV approximations and all associated sphere/cylinder
!> overlaps with sphere/cylinder functions.
!> \author David Krisiloff
!> 
!> Most of this code is a heavily refactored version of the locist_mod

!> \todo The cylinder overlaps cylinder code is INCORRECT. It only checks if the hemispherical ends overlap.

module sphere_cylinder
  use global_var_mod
  use molecule_var_mod
  use utilities_mod
  use molecular_orbital_mod
  use math_utils
  use sort_utils
  use fortran_timing
  use c_sort_finterface
  use get_basis_data
  use get_cylinder_data
  implicit none
  
  private
  
  public :: make_spheres_cylinders, sphere_overlaps_chain, chains_overlap, spheres_overlap
  public :: wp_spheres, tov_spheres, wp_cylinders, tov_cylinders, is_cylinder
  
  
  !! Data structures for the spheres and cylinders
  !***********************************************
  !> A sphere around an orbital
  type sphere
     real(real8), dimension(3)::position     ! x,y,z coordinates of the center of the sphere
     real(real8)::radius                     ! radius of the sphere
  end type sphere

  !> A cylinder with hemispherical caps. Stored as two spheres.
  type cylinder
     type(sphere)::sphere_1,sphere_2
  end type cylinder

  !> A chain of cylinders (not necessarily contiguous, although they generally will be in practice
  type chain
      integer(kind=8)::num_cylinders
      type(cylinder), dimension(:), allocatable :: cylinders
  end type chain
  
  ! Note that the following data structures are allocated in a deliberately unusual way. 
  ! (see allocate_spheres_cylinders for details)
  type(sphere), dimension(:), allocatable :: wp_spheres     !< spheres for the WP approximation
  type(sphere), dimension(:), allocatable :: tov_spheres    !< spheres for the TOV approximation
  type(chain), dimension(:), allocatable :: wp_cylinders    !< cylinders for the WP approximation
  type(chain), dimension(:), allocatable :: tov_cylinders   !< cylinders for the TOV approximation
  logical, dimension(:), allocatable :: is_cylinder         !< is this orbital a cylinder?

  !> Data struct for handling the mulliken pop and subsequent orbital analysis
  type orbital_analysis 
     real(kind=real8), dimension(:,:), allocatable :: density_matrix       !< orbital density
     real(kind=real8), dimension(:), allocatable  :: mulliken_pop          !< final mulliken population 
     real(kind=real8), dimension(:), allocatable  :: trace                 !< scratch space for the trace of the density matrix
     integer, dimension(:), allocatable :: mulliken_index                  !< index for the mulliken population
     integer :: last_important                                             !< index of the last important center
     real(kind=real8) :: max_dist                                                   !< maximum distance between the important centers
     real(kind=real8) :: center(3)                                                  !< calculated center of the orbital
  end type orbital_analysis


contains

  !> \brief Calculates all the necessary spheres and cylinders
  !>
  !> Note that this, by definition, requires that the molecular orbitals are already calculated
  subroutine make_spheres_cylinders()
    implicit none

    ! Variables
    integer :: ierr
    integer :: iOrb
    real(kind=real8) :: threshold
    real(kind=real8):: sphere_radius

    real(kind=real8), dimension(:,:), allocatable :: geometry                ! geometry of the system
    real(kind=real8), dimension(:,:), allocatable :: overlap                 ! AO basis function overlap 
    integer, dimension(:), allocatable :: basis2atom            ! maps the AO basis function to its atomic center

    type(orbital_analysis) :: orb_data                          ! holds data for the mulliken population
    type(clock) :: timer                                        ! timing data for this method



    ! Start of the routine
    !*********************
    call start_clock(timer)
    write(ioOutput,130)
    130 format(/,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/&
                ,1x,"!//                              ",/&
                ,1x,"!// ENTERED SPHERE/CYLINDER      ",/&
                ,1X,"!//       SETUP                  ",/&
                ,1x,"!//                              ",/&
                ,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/)
                
    
    ! Allocations for external info
    allocate(geometry(3,num_atoms), stat=ierr)
    call allocatecheck(ierr, "geometry")
    geometry = 0 

    allocate(overlap(number_bas,number_bas),stat=ierr)
    call allocatecheck(ierr, "overlap ")
    overlap = 0 

    allocate(basis2atom(number_bas), stat=ierr)
    call allocatecheck(ierr, "basis2atom ")
    basis2atom = 0 

    call get_coordinates(num_atoms, geometry)    ! get the geometry of the system
    call get_overlap(number_bas, overlap)        ! get the AO overlap matrix
    call allocate_spheres_cylinders()            ! allocate the module-global sphere and cylinder data structures
#ifdef TIGER_GAMESS
    call gms_basis_map(number_bas,basis2atom)    ! get the AO basis to atom center map, as per GAMESS notion
#else
    call get_basis_map(number_bas, basis2atom)    ! get the AO basis to atom center map
#endif

    ! Allocations for scratch space
    allocate( orb_data%density_matrix(number_bas, number_bas), orb_data%trace(number_bas), & 
         orb_data%mulliken_pop(num_atoms), orb_data%mulliken_index(num_atoms), &
         stat=ierr)
    call allocatecheck(ierr, "orb_data ")

    ! Here is where the real work is done. Note that we split everything into 3 loops one for occupied orbitals,
    ! one for active orbitals and one for virtual orbitals. I think this layout makes everything much clearer as
    ! opposed to having one large loop over all orbitals and lots of if statements

    ! Internal Orbitals
    !  - WP sphere
    !  - TOV sphere
    do iOrb = 1, num_internal 
        if (is_cylinder(iOrb)) then
            ! If the user manually set the endpoints use those instead
            call manual_endpoints(iOrb, geometry)
        else
            ! Analyze
            threshold = internal_threshold                               ! threshold for the importance test 
            call zero_orbital_analysis(orb_data)
            if (sphereprint) call detailed_printout_header(iOrb)
            call mulliken_pop(iOrb, orb_data, basis2atom, overlap, sphereprint)     ! do the mulliken population
            call mulliken_sort(orb_data, sphereprint)                               ! sort the mulliken population
            call find_important_centers(orb_data, threshold, sphereprint)           ! find the important centers
            call find_center(orb_data, geometry, sphereprint)                       ! find the center of the orbital
            call find_max_dist(orb_data, geometry, sphereprint)                     ! find the maximum distance between the important centers

            ! Set up the spheres
            sphere_radius = max( orb_data%max_dist * wp_multiplier , wp_default_radius )
            wp_spheres(iOrb)%radius = sphere_radius
            wp_spheres(iOrb)%position(:) = orb_data%center(:)

            sphere_radius = max( orb_data%max_dist * tov_occupied_multiplier , tov_occupied_default_radius )
            tov_spheres(iOrb)%radius = sphere_radius
            tov_spheres(iOrb)%position(:) = orb_data%center(:)
        end if
    end do
    
    ! Virtual Orbitals
    ! - TOV sphere
    do iOrb = num_internal+1 , num_orbitals

        if (is_cylinder(iOrb)) then
            ! If the user manually set the endpoints use those instead
            call manual_endpoints(iOrb, geometry)
        else
            ! Analyze
            threshold = virtual_threshold                          ! threshold for the importance test 
            if (sphereprint) call detailed_printout_header(iOrb)
            call zero_orbital_analysis(orb_data)
            call mulliken_pop(iOrb, orb_data, basis2atom, overlap, sphereprint) ! do the mulliken population
            call mulliken_sort(orb_data, sphereprint)                           ! sort the mulliken population
            call find_important_centers(orb_data, threshold, sphereprint)       ! find the important centers
            call find_center(orb_data, geometry, sphereprint)                   ! find the center of the orbital
            call find_max_dist(orb_data, geometry, sphereprint)                 ! find the maximum distance between the important centers

            ! Set up the spheres         
            sphere_radius = max( orb_data%max_dist * tov_virtual_multiplier , tov_virtual_default_radius )
            tov_spheres(iOrb)%radius = sphere_radius
            tov_spheres(iOrb)%position(:) = orb_data%center(:)
        end if
    end do

    ! print out a table of the spheres/cylinders
    call pretty_table(geometry)


    ! End of the routine ... clean up
    !********************************
    deallocate(geometry, overlap, basis2atom, stat=ierr) 
    call deallocatecheck(ierr, "geometry ...")

    deallocate(orb_data%density_matrix,orb_data%trace, orb_data%mulliken_pop, orb_data%mulliken_index,stat=ierr)
    call deallocatecheck(ierr, "orb_data% ")
    
    call print_clock(timer, "sphere&cylinder setup")
  end subroutine make_spheres_cylinders
  
  
! set up the cylinder endpoints for a particular cylinder
subroutine set_cylinder(cyl, geometry, end_1, end_2, radius)
    implicit none
    
    type(cylinder), intent(out) :: cyl
    real(kind=real8), dimension(:,:), intent(in) :: geometry
    integer(kind=8), intent(in) :: end_1, end_2
    real(kind=real8), intent(in) :: radius
  
    cyl%sphere_1%position = geometry(:,end_1)
    cyl%sphere_1%radius   = radius
    cyl%sphere_2%position = geometry(:,end_2)
    cyl%sphere_2%radius   = radius

end subroutine set_cylinder

! get manual cylinder endpoints
subroutine manual_endpoints(iOrb, geometry)
    implicit none
   
    real(kind=real8), dimension(:,:), intent(in) :: geometry
    integer(kind=8), intent(in) :: iOrb
    
    integer(kind=8) :: i, count, end_1, end_2
    integer(kind=8), dimension(:), allocatable :: cyls
    
    allocate(cyls(2*tov_cylinders(iOrb)%num_cylinders))
    call get_cylinders(iOrb, cyls)
    
    count = 1
    do i = 1, wp_cylinders(iOrb)%num_cylinders
        end_1 = cyls(count)
        count = count + 1
        end_2 = cyls(count)
        count = count + 1
        call set_cylinder(wp_cylinders(iOrb)%cylinders(i), geometry, end_1, end_2, wp_cylinder_radius)
        call set_cylinder(tov_cylinders(iOrb)%cylinders(i), geometry, end_1, end_2, tov_cylinder_radius)
    end do    
    
end subroutine manual_endpoints

  !> \brief Calculates a mulliken population for an individual orbital 
  !> \param iOrb The index of the orbital
  !> \param orb_data On exit contains the mulliken population and associated index
  !> \param basis2atom AO basis function to atom center map
  !> \param overlap  AO overlap matrix
  subroutine mulliken_pop( iOrb, orb_data, basis2atom, overlap, sphereprint)
    implicit none

    ! Inputs
    integer, intent(in) :: iOrb
    type(orbital_analysis) :: orb_data
    real(kind=real8), intent(in) :: overlap(:,:)
    integer, intent(in) :: basis2atom(:)
    logical :: sphereprint

    ! Variables 
    integer ::iBas, jBas, iAtom

    ! Zero everything just in case
    orb_data%mulliken_pop = 0.0
    orb_data%mulliken_index = 0
    orb_data%trace = 0.0 
    orb_data%density_matrix = 0.0
    
    ! If this routine was ever a performance concern the following loops could be fused
    ! (the compiler might already be doing it when optimizing)

    ! Calculate the density
    do iBas = 1, number_bas
       do jBas = 1, number_bas
          orb_data%density_matrix(iBas,jBas) = molecular_orbitals(iBas, iOrb) * molecular_orbitals(jBas, iOrb)
       end do
    end do

    ! Trace the density * overlap
    do iBas = 1, number_bas
       do jBas = 1, number_bas
          orb_data%trace(iBas) = orb_data%trace(iBas) + orb_data%density_matrix(iBas,jBas) * overlap(iBas,jBas)
       end do
    end do

    ! Generate the mulliken population by adding up the contributions on a per atom basis
    do iBas = 1, number_bas
       iAtom = basis2atom(iBas)
       orb_data%mulliken_pop(iAtom) = orb_data%mulliken_pop(iAtom) + orb_data%trace(iBas)
    end do

    ! Finally set up the index, at the moment its very simple
    do iAtom = 1 , num_atoms
       orb_data%mulliken_index(iAtom) = iAtom
    end do
    
    if (sphereprint) call detailed_printout_mulliken(orb_data)
  end subroutine mulliken_pop

  !> \brief Sorts the mulliken populations into ascending order
  !> \param orb_data On exit contains the mulliken population and associated index       
  !>
  !> This became its own routine when we decided to move towards sorting on abs()
  !> and were considering other possible algorithms
  subroutine mulliken_sort(orb_data, sphereprint)
    implicit none

    ! Inputs
    type(orbital_analysis) :: orb_data
    logical :: sphereprint

    ! Variables
    logical, dimension(:), allocatable :: negative          ! records which pops where negative before the abs()
    integer :: iAtom, ierr, jAtom
    integer, dimension(:), allocatable :: tmp, tmp2

    ! Record the negative populations
    allocate(negative(num_atoms),tmp(num_atoms),tmp2(num_atoms),stat=ierr)
    call allocatecheck(ierr,"negatives ")
    negative = .false.
    do iAtom = 1, num_atoms
       if ( orb_data%mulliken_pop(iAtom) < 0 ) negative(iAtom) = .True.
    end do

    ! Abs() the populations
    orb_data%mulliken_pop = abs( orb_data%mulliken_pop )

    ! Sort using an insertion sort (this is at most an array with ~ 100 atoms )
    do iAtom = 1, num_atoms
       tmp(iAtom) = iAtom
    end do
    call sort_real_array_with_index( orb_data%mulliken_pop, tmp)  ! tmp records the permutation created by the sort so we can update the index
    tmp2 = orb_data%mulliken_index
    do iAtom = 1, num_atoms
       orb_data%mulliken_index(iAtom) = tmp2(tmp(iAtom))
    end do
    
    ! Turn the negative populations back into negative populations
    do iAtom =1, num_atoms
       jAtom = orb_data%mulliken_index(iAtom)
       if ( negative(jAtom) ) then 
          orb_data%mulliken_pop(iAtom) = -1.0 * orb_data%mulliken_pop(iAtom)
       end if
    end do
    
    if(sphereprint) call detailed_printout_sorted_mulliken(orb_data)
    deallocate(tmp,tmp2,negative)
  end subroutine mulliken_sort

  !> \brief Finds all the most important centers for the orbital 
  !> \param orb_data data for the orbital analysis
  !> \param threshold threshold for determining which centers are significant
  subroutine find_important_centers(orb_data, threshold, sphereprint)
    implicit none

    ! Inputs
    type(orbital_analysis) :: orb_data
    real(kind=real8), intent(in) :: threshold
    logical, intent(in) :: sphereprint

    ! Variables
    real(kind=real8) :: norm_sum
    integer :: iAtom

    ! We decide which atoms are important by adding up the mulliken pops starting from the largest
    ! once we hit the threshold we stop. (Note that if we continued on the sum total 
    ! would be 1, the orbitals are normalized)
    norm_sum = 0.0
    orb_data%last_important = -1  
    do iAtom = num_atoms, 1, -1
       norm_sum = norm_sum + orb_data%mulliken_pop(iAtom)
       if ( norm_sum > threshold ) then 
          orb_data%last_important = iAtom 
          exit
       end if
    end do

    if (orb_data%last_important == -1 ) then 
       ! If we ever get here we have a serious problem ( the orbital wasn't normalized to 1 or
       ! the memory got corrupted somewhere). 
       write(*,*) "FATAL ERROR in find_important_centers."
       write(*,*) "One of these doesn't make sense, or we have memory corruption"
       write(*,*) "Calculated norm = " , norm_sum , "threshold = " , threshold
       stop
    end if
    
    if (sphereprint) call detailed_printout_import_centers(orb_data, threshold, norm_sum)
    
  end subroutine find_important_centers


  !> \brief Finds the center of the orbital by a weighted averaged of the nuclear positions
  !> \param orb_data contains the orbital analysis
  !> \param geometry the nuclear coordinates of all atomic centers
  subroutine find_center(orb_data, geometry, sphereprint)
    implicit none

    ! Inputs
    type(orbital_analysis) :: orb_data
    real(kind=real8), dimension(:,:) :: geometry
    logical, intent(in) :: sphereprint

    ! Variables
    integer :: iAtom, index
    real(kind=real8) :: partial_sum, weight

    ! We perform a weighted average of the nuclear positions of all important atomic centers
    ! To do this we take the sum of all the important contributions(partial_sum) and normalize by that amount
    orb_data%center = 0.0 
    partial_sum = sum( orb_data%mulliken_pop(orb_data%last_important:num_atoms) )
    do iAtom = num_atoms, orb_data%last_important, -1 
       weight = orb_data%mulliken_pop(iAtom) / partial_sum 
       index = orb_data%mulliken_index(iAtom)
       orb_data%center = orb_data%center + weight * geometry(:,index)          ! Note this is an array operation 
    end do
    if (sphereprint) call detailed_printout_center(orb_data)
  end subroutine find_center


  !> \brief Finds the maximum distance between the important centers in the orbital. Used for the sphere radius 
  !> \param orb_data contains the orbital analysis
  !> \param geometry the nuclear coordinates of all atomic centers
  subroutine find_max_dist(orb_data, geometry, sphereprint)
    implicit none

    ! Inputs
    type(orbital_analysis) :: orb_data
    real(kind=real8), dimension(:,:) :: geometry
    logical, intent(in) :: sphereprint

    ! Variables 
    integer :: iAtom, jAtom                          ! index for counting over atoms 
    integer :: atomA, atomB                          ! the actual atom index from the mulliken_index
    real(kind=real8) :: d

    orb_data%max_dist = -1.0
    do iAtom = orb_data%last_important, num_atoms
       atomA = orb_data%mulliken_index(iAtom)
       do jAtom = iAtom, num_atoms
          atomB = orb_data%mulliken_index(jAtom)
          d = distance( geometry(:,atomA), geometry(:,atomB) )
          if ( d > orb_data%max_dist ) orb_data%max_dist = d 
       end do
    end do
    if (sphereprint) call detailed_printout_dist(orb_data)
  end subroutine find_max_dist

  !> \brief Allocates the module-global arrays for storing the spheres and cylinders. 
  !> 
  !> Please pay careful attention to the atypical array bounds when allocating
  subroutine allocate_spheres_cylinders()
    implicit none

    ! Variables 
    integer :: ierr
    integer :: iOrb, iCyl, num_cyl

    ! Allocate the spheres and cylinders. Note that I'm always going to be accessing these 
    ! arrays by orbital number, i.e. wp_cylinders(num_internal). However I don't want to allocate
    ! unnecessary memory (not so much because of bloat but rather I want it to be explicitly clear
    ! that I don't have cylinders for virtual orbitals, etc ...) So please take care with these
    ! unusual allocation statements.

    ! note that despite this fanciness ... tov_spheres still contains spheres for the active orbitals. 
    ! THESE ARE NEVER USED.
    allocate( wp_spheres(1:num_internal), tov_spheres(1:num_orbitals), & 
        wp_cylinders(1:num_internal), tov_cylinders(1:num_orbitals), &
        is_cylinder(1:num_orbitals), stat = ierr )
    call allocatecheck(ierr, "spheres&chains")
    is_cylinder(:) = .false.
    
    do iOrb = 1, num_internal
        call zero_sphere( tov_spheres(iOrb) )
        call zero_sphere( wp_spheres(iOrb) )
        
        call num_cylinders(iOrb, num_cyl)
        if (num_cyl <= 0) cycle
        
        is_cylinder(iOrb) = .true.
        
        wp_cylinders(iOrb)%num_cylinders = num_cyl
        tov_cylinders(iOrb)%num_cylinders = num_cyl
        
        allocate(tov_cylinders(iOrb)%cylinders(num_cyl),&
                 wp_cylinders(iOrb)%cylinders(num_cyl),stat=ierr)
        call allocatecheck(ierr, "cylinders")
        
        do iCyl = 1, num_cyl
            call zero_cylinder( tov_cylinders(iOrb)%cylinders(iCyl) )
            call zero_cylinder( wp_cylinders(iOrb)%cylinders(iCyl) )
        end do
    end do
    
    do iOrb = num_internal+1, num_orbitals
        call zero_sphere( tov_spheres(iOrb) )
        
        call num_cylinders(iOrb, num_cyl)
        if (num_cyl == 0) cycle
        
        is_cylinder(iOrb) = .true.
        
        tov_cylinders(iOrb)%num_cylinders = num_cyl
        
        allocate(tov_cylinders(iOrb)%cylinders(num_cyl),stat=ierr)
        call allocatecheck(ierr, "cylinders")
        
        do iCyl = 1, num_cyl
            call zero_cylinder( tov_cylinders(iOrb)%cylinders(iCyl) )
        end do
    end do
  end subroutine allocate_spheres_cylinders


  !> \brief Helper function to initialize an orbital to zero
  !> \param A a sphere
  subroutine zero_sphere( A ) 
    implicit none
    type(sphere) :: A
    real(kind=real8) :: zero = 0.0
    A%position = zero
    A%radius = zero
  end subroutine zero_sphere

  !> \brief Helper function to initialize an orbital to zero
  !> \param A a cylinder
  subroutine zero_cylinder(A)
    implicit none
    type(cylinder) :: A
    call zero_sphere(A%sphere_1)
    call zero_sphere(A%sphere_2)
  end subroutine zero_cylinder

  !> \brief Deallocates the module-global arrays for storing the spheres and cylinders
  subroutine deallocate_spheres_cylinders
    implicit none
    integer :: ierr
    deallocate(wp_spheres, tov_spheres, wp_cylinders, tov_cylinders, stat=ierr)
    call deallocatecheck(ierr, "spheres&cylinders")
  end subroutine deallocate_spheres_cylinders


  !> \brief Determines if two spheres overlap
  !> \param sphere1 The first sphere
  !> \param sphere2 The second sphere
  logical function spheres_overlap(sphere1,sphere2)
    implicit none

    !// TESTS TO SEE IF TWO SPHERES OVERLAP
    real(real8)::sphere_distance
    type(sphere), intent(in)::sphere1,sphere2
    ! START AUTOGENERATED INITIALIZATION 
    sphere_distance = 0.0
    spheres_overlap = .false.
    ! END AUTOGENERATED INITIALIZATION 

    spheres_overlap = .false.
    sphere_distance = distance(sphere1%position, sphere2%position)
    if (sphere_distance <= (sphere1%radius + sphere2%radius)) spheres_overlap = .true.

  end function spheres_overlap


  !> \brief Determines if a single point is inside a sphere
  !> \param the_sphere The sphere
  !> \param the_point The point
  !> Please note that this is only used for producing a nice output about sphere size. It isn't used
  !> anywhere else ... which is a good thing.
  function point_in_sphere(the_sphere,the_point)
    implicit none

    !// TESTS TO SEE IF A POINT IS CONTAINED WITHIN A SPHERE
    logical::point_in_sphere

    real(real8)::sphere_distance
    type(sphere), intent(in)::the_sphere
    real(real8),dimension(3)::the_point
    ! START AUTOGENERATED INITIALIZATION 
    sphere_distance = 0.0
    point_in_sphere = .false.
    ! END AUTOGENERATED INITIALIZATION 

    point_in_sphere = .false.
    sphere_distance = distance(the_sphere%position, the_point)
    if (sphere_distance <= the_sphere%radius) point_in_sphere = .true.

  end function point_in_sphere

  !> \brief Determines if a single point is inside a cylinder
  !> \param the_cylinder The cylinder
  !> \param the_point The point
  !> Please note that this is only used for producing a nice output about cylinder size. It isn't used
  !> anywhere else ... which is a good thing. Also this logic is a bit convoluted (but props to whoever 
  !> went through the trig). There are simpler algorithms available (looking up collision detection), but
  !> for now this works.
  function point_in_cylinder(the_cylinder,the_point)
    implicit none

    !// TESTS TO SEE IF A POINT IS CONTAINED WITHIN A 
    !// CYLINDER WITH ROUNDED EDGES
    logical::point_in_cylinder

    real(real8)::sphere_distance,norm1,norm2,angle1,angle2
    real(real8)::pi,line_distance
    type(cylinder), intent(in)::the_cylinder
    real(real8),dimension(3)::the_point
    real(real8), dimension(3)::v1,v2
    ! START AUTOGENERATED INITIALIZATION 
    pi = 0.0
    angle2 = 0.0
    angle1 = 0.0
    point_in_cylinder = .false.
    sphere_distance = 0.0
    line_distance = 0.0
    norm2 = 0.0
    v2 = 0.0
    v1 = 0.0
    norm1 = 0.0
    ! END AUTOGENERATED INITIALIZATION 

    pi = 4.0d0*atan(1.0d0)

!!!!!!!!!!
    !//
    !// TEST TO SEE IF POINT IS INSIDE END
    !// SPHERES
    !//
!!!!!!!!!!
    sphere_distance = distance(the_cylinder%sphere_1%position, the_point)
    if (sphere_distance <= the_cylinder%sphere_1%radius) then
       point_in_cylinder = .true.
       return
    endif

    sphere_distance = distance(the_cylinder%sphere_2%position, the_point)
    if (sphere_distance <= the_cylinder%sphere_2%radius) then
       point_in_cylinder = .true.
       return
    endif

!!!!!!!!!!!!
    !//
    !// TEST TO SEE IF POINT IS IN
    !// INTERMEDIATE REGION BETWEEN 
    !// SPHERES
    !//
!!!!!!!!!!!!

    ! Based on graphics software methods
    point_in_cylinder = point_in_center(the_cylinder, the_point)

  end function point_in_cylinder

  ! Function to check whether a point is inbetween the two 
  ! spheres making up a cylinder
  function point_in_center(the_cylinder, the_point)
    implicit none

    logical::point_in_center

    real(kind=real8) length_squared, radius_squared, dot_prod, dist_squared
    type(cylinder), intent(in)::the_cylinder
    real(kind=real8),dimension(3), intent(in)::the_point
    real(kind=real8), dimension(3) :: cyl_1, cyl_2, xyz_dist, xyz_pt

    cyl_1 = the_cylinder%sphere_1%position
    cyl_2 = the_cylinder%sphere_2%position
    radius_squared = the_cylinder%sphere_1%radius * the_cylinder%sphere_1%radius
    xyz_dist = cyl_2 - cyl_1
    length_squared = xyz_dist(1) * xyz_dist(1) + xyz_dist(2) * xyz_dist(2) + xyz_dist(3) * xyz_dist(3)
    xyz_pt = the_point - cyl_1

    dot_prod = dot_product(xyz_pt, xyz_dist)

    if (dot_prod < 0.0 .or. dot_prod > length_squared) then
       point_in_center = .false.
       return
    else
       dist_squared = dot_product(xyz_pt, xyz_pt) - (dot_prod * dot_prod)/length_squared
       if (dist_squared > radius_squared) then
          point_in_center = .false.
          return
       else
          point_in_center = .true.
          return
       end if
    end if

  end function point_in_center

  !> \brief Determines if a sphere overlaps with a chain
  !> \param the_chain The chain
  !> \param the_sphere The sphere 
  function sphere_overlaps_chain(the_sphere, the_chain)
      implicit none
      
      logical::sphere_overlaps_chain
      integer::iCyl
      type(sphere)::the_sphere
      type(chain)::the_chain
      
      sphere_overlaps_chain = .true.
      
      do iCyl = 1, the_chain%num_cylinders
          if (sphere_overlaps_cylinder(the_sphere, the_chain%cylinders(iCyl))) return
      end do
      
      sphere_overlaps_chain = .false.
      
  end function sphere_overlaps_chain
  
  !> \brief Determines if a sphere overlaps with a cylinder
  !> \param the_cylinder The cylinder
  !> \param the_sphere The sphere
  function sphere_overlaps_cylinder(the_sphere,the_cylinder)

    !// TESTS TO SEE IF A SPHERE OVERLAPS WITH A CYLINDER

    implicit none

    logical::sphere_overlaps_cylinder
    real(real8)::line_distance
    real(real8), dimension(3)::v1,v2
    real(real8)::angle1,angle2,norm1,norm2,pi
    type(sphere), intent(in)::the_sphere
    type(cylinder), intent(in)::the_cylinder
    ! START AUTOGENERATED INITIALIZATION 
    norm1 = 0.0
    v1 = 0.0
    line_distance = 0.0
    v2 = 0.0
    norm2 = 0.0
    angle2 = 0.0
    angle1 = 0.0
    pi = 0.0
    sphere_overlaps_cylinder = .false.
    ! END AUTOGENERATED INITIALIZATION 

    pi = 4.0d0*atan(1.0d0)

!!!!!!!!!!!!!
    !//
    !// CHECK IF SPHERE OVERLAPS WITH
    !// END SPHERES OF CYLINDER
    !//
!!!!!!!!!!!!!

    sphere_overlaps_cylinder = .false.

    if (spheres_overlap(the_sphere,the_cylinder%sphere_1)) then
       sphere_overlaps_cylinder = .true.
       return
    endif

    if (spheres_overlap(the_sphere,the_cylinder%sphere_2)) then
       sphere_overlaps_cylinder = .true.
       return
    endif

!!!!!!!!!!!!!
    !//
    !// CHECK TO SEE IF THE CYLINDER OVERLAPS
    !// WITH THE MIDSECTION OF THE SPHERE
    !//
!!!!!!!!!!!!!
    sphere_overlaps_cylinder = .false.

    v1 = the_sphere%position - the_cylinder%sphere_1%position
    v2 = the_cylinder%sphere_2%position - the_cylinder%sphere_1%position
    norm1 = sqrt(dot_product(v1,v1))
    norm2 = sqrt(dot_product(v2,v2))
    angle1 = acos(dot_product(v1,v2)/(norm1*norm2))

    v1 = the_sphere%position - the_cylinder%sphere_2%position
    v2 = the_cylinder%sphere_1%position - the_cylinder%sphere_2%position
    norm1 = sqrt(dot_product(v1,v1))
    norm2 = sqrt(dot_product(v2,v2))
    angle2 = acos(dot_product(v1,v2)/(norm1*norm2))

    !// CHECK TO SEE IF POINT IS LOCATED SOMEWHERE BETWEEN A PAIR
    !// OF PLANES PASSING THROUGH THE CENTERS OF THE SPHERES
    !// AND PERPENDICULAR TO THE LINE CONNECTING THE
    !// SPHERES
    if ( ((angle1 < pi/2.0d0).or.(angle1 > 3.0d0*pi/4.0d0)).and.&
         ((angle2 < pi/2.0d0).or.(angle2 > 3.0d0*pi/4.0d0)) ) then

       !// OK, NOW CHECK TO SEE IF THE POINT IS WITHIN
       !// THE SPECIFIED RADIUS OF THE LINE CONNECTING 
       !// THE TWO SPHERES
       line_distance = norm1*sin(angle2)

       if (line_distance < the_sphere%radius + the_cylinder%sphere_1%radius) then
          sphere_overlaps_cylinder = .true.
          return
       endif

    endif
  end function sphere_overlaps_cylinder

  !> \brief Determines if two chains overlap 
  !> \param cylinder_1 The first chain
  !> \param cylinder_2 The second chain  
  function chains_overlap(chain_1, chain_2)
      implicit none
      
      type(chain), intent(in)::chain_1, chain_2
      integer::iCyl,jCyl
      logical::chains_overlap
      
      chains_overlap = .true.
      
      do iCyl = 1, chain_1%num_cylinders
          do jCyl = 1, chain_2%num_cylinders
            if (cylinders_overlap(chain_1%cylinders(iCyl), chain_2%cylinders(jCyl))) return
          end do
      end do
      
      chains_overlap = .false.
      
  end function chains_overlap
  
  !> \brief Determines if two cylinders overlap 
  !> \param cylinder_1 The first cylinder
  !> \param cylinder_2 The second cylinder
  function cylinders_overlap(cylinder_1,cylinder_2)

    !// TESTS TO SEE IF TWO CYLINDERS OVERLAP
    !// WARNING - THIS ISN'T A GENERAL TEST, IT
    !// SIMPLY CHECKS TO SEE IF THE CAPPING SPHERES
    !// ON ONE CYLINDER OVERLAPS WITH THE OTHER 
    !// CYLINDER

    implicit none

    type(cylinder), intent(in)::cylinder_1,cylinder_2
    logical::cylinders_overlap
    ! START AUTOGENERATED INITIALIZATION 
    cylinders_overlap = .false.
    ! END AUTOGENERATED INITIALIZATION 

    cylinders_overlap = .true.

    call flush(ioOutput)

    if (sphere_overlaps_cylinder(cylinder_2%sphere_1,cylinder_1)) return
    if (sphere_overlaps_cylinder(cylinder_2%sphere_2,cylinder_1)) return

    cylinders_overlap = .false.

  end function cylinders_overlap

  !> \brief Prints out a table explaining how the orbital spheres/cylinders were setup
  subroutine pretty_table( geometry)
      implicit none

      ! inputs
      real(kind=real8), dimension(:,:) :: geometry

      ! variables
      integer :: iOrb, iAtom, iCyl

      ! output formats
      701 format (" |", I9, "|",  E13.2, "|")
      702 format (" |", I9, "|",  F13.2, "|")

       write(*,*)
       write(*,*)
       write(*,*)
       write(*,*) "Sphere/Cylinder setup complete. "
       write(*,*) "All radii are in bohr"
       write(*,*) "(1) Weak Pair Spheres "
       write(*,*) "------------------------------------------------------------------"
       write(*,*) "| Orbital |   Radius    |    Atoms in the sphere "
       write(*,*) "------------------------------------------------------------------"
       do iOrb = 1, num_internal
           if (.not. is_cylinder(iOrb)) then
                ! If statement to deal with the nonlocal case of huge radii
                 if ( wp_spheres(iOrb)%radius > 100   ) then
                     write(*, 701 , advance='NO') iOrb, wp_spheres(iOrb)%radius
                 else
                     write(*, 702 , advance='NO') iOrb, wp_spheres(iOrb)%radius
                 end if
                 
                 ! Print out which atom centers fall within this sphere
                 do iAtom = 1 , num_atoms
                     if ( point_in_sphere( wp_spheres(iOrb), geometry(:,iAtom))) then
                         write(*,'("  ", I4)', advance='NO') iAtom
                     end if
                 end do
                 write(*,*) ""
           endif
       end do


       write(*,*) " "
       write(*,*) "(2) Weak Pair Cylinders "
       write(*,*) "------------------------------------------------------------------"
       write(*,*) "| Orbital |   Radius    |    Atoms in the cylinder "
       write(*,*) "------------------------------------------------------------------"
       do iOrb = 1 , num_internal
           if (is_cylinder(iOrb)) then
                do iCyl = 1, wp_cylinders(iOrb)%num_cylinders
                 ! If statement to deal with the nonlocal case of huge radii
                  if ( wp_cylinders(iOrb)%cylinders(iCyl)%sphere_1%radius > 100   ) then
                      write(*, 701 , advance='NO') iOrb, wp_cylinders(iOrb)%cylinders(iCyl)%sphere_1%radius
                  else
                      write(*, 702 , advance='NO') iOrb, wp_cylinders(iOrb)%cylinders(iCyl)%sphere_1%radius
                  end if

                  ! Print out which atom centers fall within this sphere
                  do iAtom = 1 , num_atoms
                      if ( point_in_cylinder( wp_cylinders(iOrb)%cylinders(iCyl), geometry(:,iAtom))) then
                          write(*,'("  ", I4)', advance='NO') iAtom
                      end if
                  end do
                  write(*,*) ""
                end do
           endif
       end do

       write(*,*) " "
       write(*,*) "(3) Truncation of Virtuals Spheres "
       write(*,*) "------------------------------------------------------------------"
       write(*,*) "| Orbital |   Radius    |    Atoms in the sphere "
       write(*,*) "------------------------------------------------------------------"
       do iOrb = 1, num_orbitals
           if (.not. is_cylinder(iOrb)) then
            ! If statement to deal with the nonlocal case of huge radii
             if ( tov_spheres(iOrb)%radius > 100   ) then
                 write(*, 701 , advance='NO') iOrb, tov_spheres(iOrb)%radius
             else
                 write(*, 702 , advance='NO') iOrb, tov_spheres(iOrb)%radius
             end if

             ! Print out which atom centers fall within this sphere
             do iAtom = 1 , num_atoms
                 if ( point_in_sphere(tov_spheres(iOrb), geometry(:,iAtom))) then
                     write(*,'("  ", I4)', advance='NO') iAtom
                 end if
             end do
             write(*,*) ""
           endif
       end do
       
       write(*,*) " "
       write(*,*) "(4) Truncation of Virtuals Cylinders "
       write(*,*) "------------------------------------------------------------------"
       write(*,*) "| Orbital |   Radius    |    Atoms in the cylinder "
       write(*,*) "------------------------------------------------------------------"
       do iOrb = 1, num_orbitals
           if (is_cylinder(iOrb)) then
                do iCyl = 1, wp_cylinders(iOrb)%num_cylinders           
                ! If statement to deal with the nonlocal case of huge radii
                 if ( tov_cylinders(iOrb)%cylinders(iCyl)%sphere_1%radius > 100   ) then
                     write(*, 701 , advance='NO') iOrb, tov_cylinders(iOrb)%cylinders(iCyl)%sphere_1%radius
                 else
                     write(*, 702 , advance='NO') iOrb, tov_cylinders(iOrb)%cylinders(iCyl)%sphere_1%radius
                 end if

                 ! Print out which atom centers fall within this sphere
                 do iAtom = 1 , num_atoms
                     if ( point_in_cylinder( tov_cylinders(iOrb)%cylinders(iCyl), geometry(:,iAtom))) then
                         write(*,'("  ", I4)', advance='NO') iAtom
                     end if
                 end do
                 write(*,*) ""
                end do
           endif
       end do

       write(*,*) " "
       write(*,*) " "
       write(*,*) " "
       write(*,*) " "
  end subroutine pretty_table
  
  subroutine zero_orbital_analysis(orb_data)
      implicit none
      type(orbital_analysis) :: orb_data
      orb_data%density_matrix = 0.0
      orb_data%mulliken_pop = 0.0
      orb_data%trace = 0.0
      orb_data%mulliken_index = 0
      orb_data%last_important = 0
      orb_data%max_dist = 0.0
      orb_data%center = 0.0
  end subroutine zero_orbital_analysis
  
  
  subroutine detailed_printout_header(orb_number)
      implicit none
      integer, intent(in) :: orb_number
      201 format ("Producing detailed information for orbital ", I9)
      199 format ("****************************************************")
      write(*,*) 
      write(*,199)
      write(*,201) orb_number
      write(*,199)
  end subroutine detailed_printout_header
      
  subroutine detailed_printout_mulliken(orb_data)
      implicit none
      integer :: i
      type(orbital_analysis), intent(in) :: orb_data
      202 format ("-----------------------------------------------------")
      203 format ("Atom ", I9, " has mulliken population = " , E15.8)
      write(*,202)
      do i = 1, size(orb_data%mulliken_pop)
          write(*,203) i, orb_data%mulliken_pop(i)
      end do
  end subroutine detailed_printout_mulliken
  
  subroutine detailed_printout_import_centers(orb_data, threshold, total)
      implicit none
      type(orbital_analysis), intent(in) :: orb_data
      real(kind=real8), intent(in) :: threshold, total
      integer :: i 
      204 format ("-----------------------------------------------------")
      205 format ("There are ", I9, " important atom(s) for this sphere")
      206 format ("Atom ", I9)
      210 format ("Threshold for importance is ", F8.2, "%, reached ", F8.2, "% with last atom")
      write(*,204)
      write(*,205) num_atoms  - orb_data%last_important + 1
      do i = num_atoms, orb_data%last_important, -1
          write(*,206) i
      end do
      write(*,210) threshold*100, total*100
  end subroutine detailed_printout_import_centers
  
  subroutine detailed_printout_sorted_mulliken(orb_data)
      implicit none
      type(orbital_analysis), intent(in) :: orb_data
      integer :: i,j
      207 format ("-----------------------------------------------------")
      208 format ("Atoms listed in descending importance (absolute mulliken charge)")
      209 format ("Atom ", I9)
      write(*,207)
      write(*,208)
      do i = size(orb_data%mulliken_index), 1, -1
          j = orb_data%mulliken_index(i)
          write(*,209) j 
      end do
  end subroutine detailed_printout_sorted_mulliken
  
  subroutine detailed_printout_center(orb_data)
      implicit none
      type(orbital_analysis), intent(in) :: orb_data
      211 format ("-----------------------------------------------------")
      212 format ("Sphere located at x= ", F7.2, " y = " , F7.2, " z = " ,F7.2)
      write(*,211)
      write(*,212) orb_data%center(1), orb_data%center(2), orb_data%center(3)
  end subroutine detailed_printout_center
  
  subroutine detailed_printout_dist(orb_data)
      implicit none
      type(orbital_analysis), intent(in) :: orb_data
      213 format ("-----------------------------------------------------")
      214 format ("Radius of the sphere (before default size / multipliers) = ", F7.2) 
      write(*,213)
      write(*,214) orb_data%max_dist
  end subroutine detailed_printout_dist
  
                    
end module sphere_cylinder

