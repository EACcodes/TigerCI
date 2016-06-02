! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module initial_davidson_struct
  use global_var_mod
  use sparse_mod
  use davidson_io_mod
  implicit none

  type initial_dav_guess
     logical :: from_disk
     integer :: nroots, unit_num
     type(sparse_vector), dimension(:), allocatable :: vecs
     real(kind=real8), dimension(:), allocatable :: energies
  end type initial_dav_guess

contains


  !> \brief Returns the reference energy
  real(kind=real8) function get_reference_energy( obj, root ) 
    implicit none
    type(initial_dav_guess) :: obj
    integer :: root
    if ( .not. allocated(obj%energies) ) then 
       write(*,*) "fatal error ... I never computed the reference energies "
       stop
    else if ( size(obj%energies) < root ) then 
       write(*,*) "fatal error"
       write(*,*) "you asked for the energy of root = " , root
       write(*,*) "but I only have ", size(obj%energies), "roots"
       stop
    end if
    if ( allocated(user_ref_energy)) then 
       get_reference_energy = user_ref_energy(root)
    else
       get_reference_energy = obj%energies(root)
    end if
  end function get_reference_energy

  !> \brief Returns true if we are using a user provided reference energy
  logical function using_user_reference_energy()
    using_user_reference_energy =  allocated(user_ref_energy)
  end function using_user_reference_energy


  !> \brief Copies the initial vector (for the given root) into the given array 
  !> \param guess The initial_dav_guess storing the initial guess
  !> \param dest  The destination for the initial guess
  !> \param root  The root I want the initial guess for
  subroutine pull_initial_guess(guess, dest, root)
    implicit none

    ! inputs
    type(initial_dav_guess) :: guess
    real(kind=real8), dimension(:) :: dest
    integer, intent(in) :: root

    ! Two options, either the guess is in memory (as a sparse vector) or on disk
    if ( guess%from_disk ) then 
       call davidson_io_read(dest, ci_final, root) 
   else
       if (size(dest) <= 0 ) call cause_floating_point_exception
       call inflate_sparse_vector(guess%vecs(root), dest)
    end if
  end subroutine pull_initial_guess

  !> \brief Frees the memory associated with an initial_dav_guess type
  subroutine deallocate_initial_dav_guess(guess)
    implicit none

    ! inputs 
    type(initial_dav_guess) :: guess

    ! local variables
    integer :: i , n

    if ( allocated(guess%vecs) ) then 
       n = size(guess%vecs)
       do i = 1, n 
          call destroy_sparse_vector(guess%vecs(i))
       end do
       deallocate(guess%vecs)
    end if

    if ( allocated(guess%energies) ) then 
       deallocate(guess%energies)
    end if

  end subroutine deallocate_initial_dav_guess
end module initial_davidson_struct
