! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> \brief Handles the multiplication of G*c where G is the matrix of g values used in ACPF/ACPF2/AQCC calculations
module gtimesc_mod
  use utilities_mod
  use global_var_mod
  use new_tree_search_mod
  use allowed_virtuals_mod
  implicit none


  ! The compressed G matrix
  ! We store the diagonal of G as a compressed vector (everything else is zero)
  ! 
  ! ex: 111333 is stored as 
  ! 1 3
  ! 3 3
  ! where the first line (value array) stores the values
  ! and the second line (repeat array) stores the number of times the value 
  ! (in the same position) is repeated
  type CompressedG
     integer :: head
     real(kind=real8), dimension(:), allocatable :: value
     integer, dimension(:), allocatable :: repeat
  end type CompressedG

  type(CompressedG) :: g

  private
  public :: build_compressed_g, gtimesc, ginvtimesc


contains

  !> \brief Builds the compressed G matrix 
  !> \param work Some scratch space 
  subroutine build_compressed_g(work)
    ! inputs
    real(kind=real8), dimension(:) :: work 

    ! local variables
    real(kind=real8), dimension(5) :: g_values
    integer :: i

    ! first step setup the 5 g values
    call get_gvalues(g_values)

    ! next we are going build the diagonal of the g matrix
    call build_g_diagonal(g_values, work)

    ! now we compress the diagonal elements
    call allocate_compressed_g(g)
    do i = 1, total_csfs
       call add_val_to_compressed_g(work(i), g)
    end do
  end subroutine build_compressed_g

  !> \brief Performs the multiplication G*c
  !> \param c The vector being multiplied. On exit contains the product G*c
  subroutine gtimesc(c)
    ! inputs
    real(kind=real8), dimension(:) :: c

    ! local variable
    integer :: i,j,p
    real(kind=real8) :: g_value

    p = 1
    do i = 1, g%head
       g_value = g%value(i)
       do j = 1, g%repeat(i)
          c(p) = c(p) * g_value
          p = p + 1
       end do
    end do
  end subroutine gtimesc

  !> \brief Performs the multiplication G^(-1)*c
  !> \param c The vector being multiplied. On exit contains the product G^(-1)*c
  subroutine ginvtimesc(c)
    ! inputs
    real(kind=real8), dimension(:) :: c

    ! local variable
    integer :: i,j,p
    real(kind=real8) :: g_value

    p = 1
    do i = 1, g%head
       g_value = g%value(i)
       do j = 1, g%repeat(i)
          c(p) = c(p) / g_value
          p = p + 1
       end do
    end do
  end subroutine ginvtimesc


  !> \brief This routine calculates the diagonal elements of the G matrix.
  !> \param g  A list of the 5 g values
  !> \param diag On exit contains the diagonal elements
  !>
  !>
  !> In an orthonormal basis the diagonal elements are just the g values. All
  !> we do here is loop through the wavefunction figure out which configurations
  !> get which g values.
  subroutine build_g_diagonal(g, diag)

    ! inputs
    real(kind=real8), dimension(5) :: g
    real(kind=real8), dimension(:) :: diag

    ! local variables 
    integer::i,j                                !// loop control variables
    integer:: the_csf_start
    integer:: the_csf_end
    integer:: num_virtc2, num_virt
    type(orbital_path)::lambda_path
    type(graph_search_state) :: graph


    !// first we will treat all parts of the ci vector that require multiplication
    !// by g1, g2, or g3.  these parts correspond to csfs that have no spatial 
    !// occupation in the external space.  tree search will direct this computation 
    !// to the routine, "gtimesc_no_external".


    !// for all csfs w/o n-1 and n-2 csfs 
    call allocate_orbital_path(lambda_path)
    call init_tree_search(graph, 0, 0, num_internal, num_elec)
    do while ( get_internal_path(lambda_path, graph))
       call g_diag_no_external(lambda_path, diag, g)
    end do


    !// now we treat the singly excited csfs with the corresponding (g4) value.
    !// start by looping over the internal pathes that have one less electron 
    !// than does the reference.

    do i=1,size(internal_index_vector1) 
       if (internal_index_vector1(i) <= 0) cycle  
       num_virt = num_allowed_virtuals(i,"s")
       do j=1,fsn(nm1_singles(i)+1)
          the_csf_start = internal_index_vector1(i) + (j - 1)*num_virt
          the_csf_end   = the_csf_start + num_virt - 1
          diag(the_csf_start:the_csf_end) = g(4)
       enddo
    enddo


    !// treat the doubly excited csfs in the same fashion.  now we use the
    !// g5 value.  there are two cases here: 1) two virtuals singly occupied
    !// and 2) one virtual doubly occupied.  

    !// first case - two virtuals singly occupied

    do i=1,size(internal_index_vector2) 
       if (internal_index_vector2(i) <= 0) cycle  
       num_virt = num_allowed_virtuals(i,"d")
       num_virtc2 = num_virt*(num_virt-1)/2
       do j=1,fsn(nm2_singles(i)+2)
          the_csf_start = internal_index_vector2(i) + (j - 1)*num_virtc2 
          the_csf_end  = the_csf_start + num_virtc2 - 1
          diag(the_csf_start:the_csf_end) = g(5)
       enddo
    enddo

    !// second case - one virtual doubly occupied

    do i=1,size(internal_index_vector3) 
       if (internal_index_vector3(i) <= 0) cycle  
       num_virt = num_allowed_virtuals(i,"d")
       do j=1,fsn(nm2_singles(i))
          the_csf_start = internal_index_vector3(i) + (j - 1)*num_virt 
          the_csf_end  = the_csf_start + num_virt - 1 
          diag(the_csf_start:the_csf_end) = g(5)
       enddo
    enddo


    !// clean up before quitting.
    call deallocate_orbital_path(lambda_path)
  end subroutine build_g_diagonal

  !> in this subroutine we have the internal configuration path and
  !> there are no excitations into the external space. therefore we  
  !> determine the type of g value associated with the configuration and 
  !> it will be one of these three options:
  !>    g(1) -  reference csf
  !>    g(2) -  excitation only within the active space
  !>    g(3) -  excitation from internal to active space
  !*****************************************************************
  subroutine g_diag_no_external(lambda_path, diag, g)
    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  variable declaration !!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(orbital_path)::lambda_path
    real(kind=real8), dimension(:) :: diag
    real(kind=real8), dimension(5) :: g


    integer::start_elec             !// number of electrons in the internal space
    integer::internal_weight        !// weight of internal part of path
    integer::singles                !// number of singles in path
    integer::csf_address            !// address of internal space csf 
    integer::address                !// address of internal space csf 
    integer::lambda_dim             !// dimensions of mu paths
    integer::i,j                    !// loop variables
    integer::spin_lambda            !// for looping over spin functions
    real(real8)::two                !// the number two

    logical::ref_flag
    logical::test_flag
    ! START AUTOGENERATED INITIALIZATION 
    address = 0
    j = 0
    i = 0
    two = 0.0
    singles = 0
    ref_flag = .false.
    lambda_dim = 0
    csf_address = 0
    test_flag = .false.
    start_elec = 0
    spin_lambda = 0
    internal_weight = 0
    ! END AUTOGENERATED INITIALIZATION 



!!!!!!!!!!!!!!!!!!!!!!! variable initialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    start_elec = sum(lambda_path%occupations(0:num_internal))
    internal_weight = sum(lambda_path%arc_weights(0:num_internal))
    singles = lambda_path%num_singles
    !// set up dimensions
    lambda_dim = fsn(singles)
    two = real(2.0, real8)

!!!!!!!!!!!!!!!!!!!!!!!!!!! start of routine !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !// double check to make sure we have no electrons in our virtual orbitals
    if (start_elec == num_elec - 0) then

       !// set up the address for our path
       address = internal_index_vector0(internal_weight+1)

       if (address > 0 .and. lambda_dim > 0) then



          !// check to see if we have a type (3) g value
          if (sum(lambda_path%occupations(0:num_inactive)) .ne. num_inactive * 2) then
             !// if you get here then you have a g(3) type configuration

             do spin_lambda = 1, lambda_dim
                csf_address = address + spin_lambda - 1 

                diag(csf_address) = g(3) 
             enddo

          else
             !// we now have to see if we are dealing with a g(1) or g(2) type 
             !// configuration
             ref_flag = .false.
             do i=1, num_ref
                test_flag = .true.
                do j = num_inactive+1, num_internal
                   if (lambda_path%occupations(j) .ne. references(j,i)) then
                      test_flag = .false.
                   endif
                enddo
                if (test_flag) then
                   ref_flag = .true.
                endif
             enddo


             if (ref_flag) then
                do spin_lambda = 1, lambda_dim
                   csf_address = address + (spin_lambda-1) 
                   diag(csf_address) = g(1) 
                enddo
             else
                do spin_lambda = 1, lambda_dim
                   csf_address = address + (spin_lambda-1) 

                   diag(csf_address) = g(2) 
                enddo
             endif


          endif


       endif

    endif

  end subroutine g_diag_no_external


  !> this subroutine sets the g values based on the choosen calculation type.
  !> the values are stored in the g variable
  !*****************************************************************
  subroutine get_gvalues(values)


!!!!!!!!!!!!!!!!!!!!!!!!!!!  variable declaration !!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    real(kind=real8), dimension(5) :: values
    real(real8)::zero                     !//  the number zero
    real(real8)::one                      !//  the number one 
    real(real8)::two                      !//  the number two 
    real(real8)::four                     !//  the number four
    real(real8)::n                        !//  number of correlated electrons 
    integer :: num_corr_elec

!!!!!!!!!!!!!!!!!!!!!!! variable initialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zero = real(0.0,real8)      
    one  = real(1.0,real8)      
    two  = real(2.0,real8)      
    four = real(4.0,real8)      


    !*********************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start of routine !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*********************************************************************************

    num_corr_elec =  num_elec -  (2 * num_frozen) 
    n = real(num_corr_elec,real8)

    if (acpf_flag   ==  1) then 
       values(1:2)   = one
       values(3:5) =  two/n 
       write(iooutput,*) "acpf method"
       call flush(iooutput)

    elseif (acpf_flag ==  2) then   
       values(1:2) =  one 
       values(3:4) =  (four/n)*(one - (one/(two*(n-one)))) 
       values(5)   =  two/n
       write(iooutput,*) "acpf-2 method"
       call flush(iooutput)

    elseif (acpf_flag ==  9) then   
       values(1:5) =  one
       write(iooutput,*) "acpf test routine"
       call flush(iooutput)

    elseif (acpf_flag == 3) then
       values(1:2) = one
       values(3:5) =  (four/n)*(one - (one/(two*(n-one))))
       write(iooutput,*) "mraqcc method "
       call flush(iooutput)            

    elseif (acpf_flag ==4 ) then
       values(1:2) = one
       values(3:4) = custom_g_val
       values(5) = two/n
       write(iooutput,*) "custom MRACPF type functional "
       call flush(iooutput) 

    endif

    !// write the g values to the output file (for debugging and identification)
    write(iooutput,140) values
140 format(/,1x,"the following acpf g values have been used:",/&
         ,1x,"                                           ",/&
         ,1x,"g1 -",f8.5,/&
         ,1x,"g2 -",f8.5,/&
         ,1x,"g3 -",f8.5,/&
         ,1x,"g4 -",f8.5,/&
         ,1x,"g5 -",f8.5,/&
         ,1x," ",/)
  end subroutine get_gvalues

  !> \brief Allocates the compressed G matrix
  !> \param G The compressed G matrix
  subroutine allocate_compressed_g(g)
    type(CompressedG) :: g
    integer :: ierr
    allocate(g%value(1), g%repeat(1), stat=ierr)
    call allocatecheck(ierr, "compressed g matrix")
    g%value = 0
    g%repeat = 0 
    g%head = 0 
  end subroutine allocate_compressed_g

  !> \brief Adds a new value to a compressed G matrix
  !> \param val  The new value
  !> \param g    The compressed G matrix
  subroutine add_val_to_compressed_g(val, g)
    ! inputs
    type(CompressedG) :: g
    real(kind=real8), intent(in) :: val

    ! if we are just starting
    if (g%head == 0 ) then 
       g%head = 1
       g%value(1)  = val
       g%repeat(1) = 1
    else
       ! if this value is equivalent to the current head of the list 
       ! just increment the repeat counter
       if (val == g%value(g%head)) then 
          g%repeat(g%head) = g%repeat(g%head) + 1
       else
          ! we are adding a new value
          g%head = g%head + 1
          if (g%head > size(g%value)) call expand_compressed_g(g)
          g%value(g%head) = val
          g%repeat(g%head) = 1 
       end if
    end if
  end subroutine add_val_to_compressed_g

  !> \brief Doubles the size of the compressed G matrix
  !> \param G The compressed G matrix
  subroutine expand_compressed_g(g) 
    ! input
    type(CompressedG) :: g

    ! local variables
    type(CompressedG) :: new
    integer :: ierr, old_size, new_size

    old_size = size(g%value)
    new_size = old_size * 2

    allocate(new%value(new_size), new%repeat(new_size), stat=ierr)
    call allocatecheck(ierr, "reallocation of the compressed g matrix")
    new%value = 0 
    new%repeat = 0 

    new%value(:old_size) = g%value(:)
    new%repeat(:old_size) = g%repeat(:)

    call move_alloc(new%value, g%value)
    call move_alloc(new%repeat, g%repeat)
  end subroutine expand_compressed_g


end module gtimesc_mod



