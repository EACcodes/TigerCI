! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> \brief This module provides a centralized interface to the WP and TOV approximations
!> \author David Krisiloff
!> \date 2/2012
!>
!> Most of this module is simply routines collected from "semi-random" locations throughout Tiger

module wp_tov_mod
  use global_var_mod
  use molecule_var_mod
  use utilities_mod
  use sphere_cylinder

  public :: build_ignorable_pair_matrix, are_interacting_WP, are_interacting_TOV

  logical,dimension(:,:),allocatable, public::ignorable_pair


contains
  !**********************************************************************************************
  !> \brief Generates the WP/TOV matrix 
  !> 
  !> This is a weird data structure. If for two orbitals ignorable_pair(p,q) is true it means that 
  !> when considering all double excitations that p and q are never in the same double excitation
  !> If they are both interal orbitals then (p,q) is removed during the WP approximation 
  !> If p is internal and q virtual than there is no hole pair (i,p) which interacts with q
  !> If both are virtuals then there is no hole pair (i,j) which interacts with both p and q
  subroutine build_ignorable_pair_matrix()
    implicit none
    integer :: i,j,k,a,b,ierr
    integer :: num_active_start
    logical :: i_interacts, j_interacts


    allocate(ignorable_pair(num_orbitals,num_orbitals),stat=ierr)
    call allocatecheck(ierr, "ignorable_pair  ")


    if (sphere_based_integral_truncations) then
       ignorable_pair = .true.
       num_active_start = num_internal - num_active + 1
       do i = 1, num_internal
          ! Do the weak pairs test
          do j = 1,i

             if ((.not. is_cylinder(i)) .and. (.not. is_cylinder(j))) then  ! i,j both spheres

                if (spheres_overlap(wp_spheres(i), wp_spheres(j))) then
                   ignorable_pair(i,j) = .false.
                   ignorable_pair(j,i) = .false.

                   do a = num_internal+1,num_orbitals

                      ! i and j interact under the WP approximation. if the virtual interacts with i or j 
                      ! we say interacts with the pair (ij)
                      i_interacts = spheres_overlap(tov_spheres(i),tov_spheres(a))
                      j_interacts = spheres_overlap(tov_spheres(j),tov_spheres(a))
                      if (i_interacts .or. j_interacts) then
                         ignorable_pair(i,a) = .false.
                         ignorable_pair(a,i) = .false.
                         ignorable_pair(j,a) = .false.
                         ignorable_pair(a,j) = .false.

                         do b = num_internal+1,a
                            i_interacts = spheres_overlap(tov_spheres(i),tov_spheres(b))
                            j_interacts = spheres_overlap(tov_spheres(j),tov_spheres(b))
                            if (i_interacts .or. j_interacts) then
                               ignorable_pair(a,b) = .false.
                               ignorable_pair(b,a) = .false.
                            endif

                         enddo

                      endif
                   enddo

                endif

             else if (.not. is_cylinder(i)) then  ! i sphere, j cylinder

                if ( sphere_overlaps_chain( wp_spheres(i), wp_cylinders(j)) ) then
                   ignorable_pair(i,j) = .false.
                   ignorable_pair(j,i) = .false.

                   do a = num_internal+1,num_orbitals
                      i_interacts = sphere_overlaps_chain( tov_spheres(a), tov_cylinders(j) )
                      j_interacts = spheres_overlap(tov_spheres(i),tov_spheres(a))

                      if(i_interacts .or. j_interacts) then
                         ignorable_pair(i,a) = .false.
                         ignorable_pair(a,i) = .false.
                         ignorable_pair(j,a) = .false.
                         ignorable_pair(a,j) = .false.

                         do b = num_internal+1,a
                            i_interacts = sphere_overlaps_chain( tov_spheres(b), tov_cylinders(j) )
                            j_interacts = spheres_overlap(tov_spheres(i),tov_spheres(b))
                            if (i_interacts .or. j_interacts) then
                               ignorable_pair(a,b) = .false.
                               ignorable_pair(b,a) = .false.
                            endif

                         enddo

                      endif
                   enddo

                endif
             else if (.not. is_cylinder(j)) then  ! j sphere, i cylinder

                if ( sphere_overlaps_chain( wp_spheres(j), wp_cylinders(i)) ) then
                   ignorable_pair(i,j) = .false.
                   ignorable_pair(j,i) = .false.

                   do a = num_internal+1,num_orbitals
                      i_interacts = sphere_overlaps_chain( tov_spheres(a), tov_cylinders(i) )
                      j_interacts = spheres_overlap(tov_spheres(j),tov_spheres(a))

                      if(i_interacts .or. j_interacts) then
                         ignorable_pair(i,a) = .false.
                         ignorable_pair(a,i) = .false.
                         ignorable_pair(j,a) = .false.
                         ignorable_pair(a,j) = .false.

                         do b = num_internal+1,a
                            i_interacts = sphere_overlaps_chain( tov_spheres(b), tov_cylinders(i) )
                            j_interacts = spheres_overlap(tov_spheres(j),tov_spheres(b))
                            if (i_interacts .or. j_interacts) then
                               ignorable_pair(a,b) = .false.
                               ignorable_pair(b,a) = .false.
                            endif

                         enddo

                      endif
                   enddo

                endif
                
             else ! both i and j are cylinders

                if ( chains_overlap( wp_cylinders(i), wp_cylinders(j) )) then 
                   ignorable_pair(i,j) = .false.
                   ignorable_pair(j,i) = .false.

                   do a = num_internal+1,num_orbitals
                      i_interacts = sphere_overlaps_chain( tov_spheres(a), tov_cylinders(i) )
                      j_interacts = sphere_overlaps_chain( tov_spheres(a), tov_cylinders(j) )
                      if ( i_interacts .or. j_interacts ) then
                         ignorable_pair(i,a) = .false.
                         ignorable_pair(a,i) = .false.
                         ignorable_pair(j,a) = .false.
                         ignorable_pair(a,j) = .false.

                         do b = num_internal+1,a
                            i_interacts = sphere_overlaps_chain( tov_spheres(b), tov_cylinders(i) )
                            j_interacts = sphere_overlaps_chain( tov_spheres(b), tov_cylinders(j) )
                            if ( i_interacts .or. j_interacts) then
                               ignorable_pair(a,b) = .false.
                               ignorable_pair(b,a) = .false.
                            endif

                         enddo

                      endif
                   enddo

                endif ! if i/j not wp

             endif  ! if both i/j are active

          enddo

       enddo
    else
       ignorable_pair = .false.
    end if


    ! POSSIBLE BUG CHECK
    if ( nonlocal_flag == 1)  then
       do i=1, num_orbitals
          do j=1 , num_orbitals
             if ( ignorable_pair(i,j) ) then
                write(*,*) "THIS IS A NONLOCAL CALCULATION"
                write(*,*) "HOWEVER THE ignorable_pair MATRIX INDICATES TRUNCATION !!"
                write(*,*) "FATAL ERROR"
                stop
             end if
          end do
       end do
    end if

    ! ANOTHER BUG CHECK
    do i=1, num_orbitals
       if ( ignorable_pair(i,i) ) then
          write(*,*) "ITS A REALLY BAD SIGN WHEN AN ORBITAL DOESN'T CORRELATE WITH ITSELF"
          write(*,*) "(EVEN WHEN USING THE LOCAL APPROXIMATION)"
          write(*,*) "FIRST OCCURANCE IS WITH ORBITAL " , i
          write(*,*) "FATAL ERROR"
          stop
       end if
    end do


    ! SOME USEFUL OUTPUT ABOUT HOW WELL LOCALIZED THINGS ARE
    if (nonlocal_flag == 0) then 
       open(unit=777, file=scratch_directory // "wp_tov.dat")
701    format (" ",I6, "|" )
702    format (" ",1x,I6)
703    format (" ",I6,1x,I6, "|")
       write(777,*) 
       write(777,*) "Local Correlation Summary"
       write(777,*) 
       write(777,*) "(1) Weak Pairs"
       write(777,*) "---------------------------------------"
       write(777,*) " Orbital | Correlating Orbitals "
       write(777,*) "---------------------------------------"
       do i = 1 , num_internal
          write(777,701,advance='no') i 
          do j = 1 , num_internal
             if (are_interacting_WP(i,j)) then
                write(777,702, advance='no') j
             end if
          end do
          write(777,*) ''
       end do

       write(777,*) 
       write(777,*) "(2) Truncation of Virtuals "
       write(777,*) "---------------------------------------"
       write(777,*) " Orbital Pair | Correlating Orbitals "
       write(777,*) "---------------------------------------"
       do i = 1, num_internal
          do j = i, num_internal
             if (are_interacting_WP(i,j)) then
                write(777,703,advance='no') i,j
                do k = num_internal+1 , num_orbitals
                   if (are_interacting_TOV(i,j,k)) then 
                      write(777,702,advance='no') k 
                   end if
                end do
                write(777,*) ''
             end if
          end do
       end do
    end if
  end subroutine build_ignorable_pair_matrix


  logical function are_interacting_WP(i,j)
    implicit none

    ! inputs
    integer, intent(in) :: i,j

    ! local variables
    logical :: i_active, j_active

    if ( (i>num_internal) .or. (j>num_internal) ) then
       write(*,*) "are_interacting_WP was called with virtual orbitals !!!!!!"
       write(*,*) "the weak pairs approximation applies only to internal orbitals"
       stop
    end if

    ! determine if i or j are active orbitals
    i_active = (i > num_inactive)
    j_active = (j > num_inactive)

    ! we check if the spheres/cylinders overlap under the Weak Pairs approximation
    if ( is_cylinder(i) .and. is_cylinder(j) ) then 
       are_interacting_WP = chains_overlap(wp_cylinders(i), wp_cylinders(j))
    elseif (is_cylinder(i)) then
       are_interacting_WP = sphere_overlaps_chain(wp_spheres(j), wp_cylinders(i))
    elseif (is_cylinder(j)) then 
       are_interacting_WP = sphere_overlaps_chain(wp_spheres(i), wp_cylinders(j))
    else
       are_interacting_WP = spheres_overlap(wp_spheres(i), wp_spheres(j))
    end if

  end function are_interacting_WP

  logical function are_interacting_TOV(i,j,virt)
    implicit none

    ! inputs
    integer, intent(in) :: i,j, virt

    ! local variables
    logical :: i_active, j_active

    if ( (i>num_internal) .or. (j>num_internal) ) then
       write(*,*) "are_interacting_TOV was called with virtual orbitals as the internal pair !!!!!!"
       write(*,*) "the internal pair needs to be made of internal orbitals"
       stop
    end if

    if ( virt < num_internal ) then 
       write(*,*) "are_interacting_TOV was called with a virtual orbital that is in fact an internal orbital"
       write(*,*) "that makes no sense"
       stop
    end if

    ! determine if i or j are active orbitals
    i_active = (i > num_inactive)
    j_active = (j > num_inactive)

    ! we check if either the cylinder/sphere for i or j interact with the sphere for virt
    if ( is_cylinder(i) .and. is_cylinder(j) ) then 
       are_interacting_TOV = sphere_overlaps_chain(tov_spheres(virt), tov_cylinders(i)) .or.  &
            sphere_overlaps_chain(tov_spheres(virt), tov_cylinders(j))
    elseif ( is_cylinder(i) ) then 
       are_interacting_TOV = sphere_overlaps_chain(tov_spheres(virt), tov_cylinders(i)) .or.  &
            spheres_overlap(tov_spheres(virt), tov_spheres(j))
    elseif ( is_cylinder(j) ) then 
       are_interacting_TOV = spheres_overlap(tov_spheres(virt), tov_spheres(i)) .or.  &
            sphere_overlaps_chain(tov_spheres(virt), tov_cylinders(j))
    else
       are_interacting_TOV = spheres_overlap(tov_spheres(virt), tov_spheres(i)) .or.  &
            spheres_overlap(tov_spheres(virt), tov_spheres(j))
    end if
  end function are_interacting_TOV

end module wp_tov_mod
