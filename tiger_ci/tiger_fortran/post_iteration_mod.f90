! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> POST_ITERATION_MOD
!> 
!> This module performs an analysis of the "best guess" for the wavefunction 
!> at the end of a davidson iteration. There are two public interfaces
!>
!> post_iteration_driver:         performs the analysis of the best guess for the wavefunction
!> write_post_iteration_analysis: writes the analysis to the output file
!>
!>
!> This might seem like an odd way to do things. Why do we have two routines for what 
!> would usually be the job of one? It moves things into two different tasks, performing the analysis
!> and writing it to the output file. As it turns out, I'd prefer to write the analysis at a different
!> time then when we have all the information we need to perform the analysis.

module post_iteration_mod

  use global_var_mod
  use tree_search_mod
  use molecule_var_mod
  use graph_var_mod
  use allowed_virtuals_mod
  use cholesky_structs
  use new_tree_search_mod
  use new_tree_search_structs
  use ci_utilities_mod
  use io_unit_numbers
  use sort_utils
  use c_sort_finterface
  use davidson_io_mod
  implicit none

  private
  public :: post_iteration_driver, write_post_iteration_analysis

  type CSF_data_table
     real(real8), dimension(:), allocatable::top_coeff
     integer, dimension(:,:), allocatable::top_csf

     integer, dimension(:,:), allocatable::ref_csf
     real(real8), dimension(:), allocatable::ref_coeff
     integer:: num_ref_csf  
  end type CSF_data_table

  integer, parameter :: NUM_ENTRIES = 10          ! total number of entries in the output 

contains

  !> \brief Runs the analysis for an iteration of our Davidson solver
  !> \param wavefunc  The best guess for the wavefuncton from the current iteration of the solver
  !>
  !> This does two things. First it svaes a copy of the wavefunction to disk. Second it 
  !> produces a little analysis of the most important configurations in the wavefunction.
  subroutine post_iteration_driver(wavefunc, root)
    ! input
    real(kind=real8), dimension(:) :: wavefunc
    integer, intent(in) :: root

    !local variable
    type(CSF_data_table) :: table

    ! allocate some space to work in
    call allocate_data_table(table, NUM_ENTRIES, 2)

    ! write out a copy of wavefunc
    call davidson_io_write(wavefunc, ci_final, root)

    ! work out this iteration's analysis
    call single_CSF_analysis(wavefunc, table)

    ! deallocate our work space
    call deallocate_data_table(table)
  end subroutine post_iteration_driver

  !> \brief Allocates a CSF_data_table
  !> \param table The table to be allocated
  subroutine allocate_data_table(table, top_size, ref_size)
    ! inputs 
    type(CSF_data_table) :: table
    integer, intent(in) :: top_size, ref_size

    ! local variables 
    integer :: ierr

    allocate(table%top_csf(top_size, num_internal+4), table%top_coeff(top_size), stat=ierr)
    call allocatecheck(ierr, "top_coeff, top_csf")

    allocate(table%ref_csf(ref_size,num_internal+4), table%ref_coeff(ref_size), stat=ierr)
    call allocatecheck(ierr, "ref_csf, ref_coeff")

    table%top_csf = 0
    table%top_coeff = 0 
    table%ref_csf = 0
    table%ref_coeff = 0
    table%num_ref_csf = 0 
  end subroutine allocate_data_table

  !> \brief Deallocated a CSF_data_table
  !> \table The table
  subroutine deallocate_data_table(table)
    ! input
    type(CSF_data_table) :: table

    ! local variable
    integer :: ierr

    deallocate(table%top_csf, table%top_coeff, table%ref_coeff, table%ref_csf, stat=ierr)
    call deallocatecheck(ierr, "a CSF_data_table")
  end subroutine deallocate_data_table

  !> \brief Produces the analysis of the most important configurations in the wavefunction
  !> \param wavefunc The wavefunction
  !> \param table    A data table (scratch space)
  subroutine single_CSF_analysis(wavefunc, table)
    ! inputs
    real(kind=real8), dimension(:) :: wavefunc
    type(CSF_data_table) :: table

    ! local variables
    integer::i
    type(graph_search_state) :: G
    type(orbital_path) :: the_path
    logical :: opened_flag

    ! allocate the path 
    call allocate_orbital_path(the_path)

    ! search through all the configurations in the wavefunctions
    ! to find the most important
    do i = num_elec-2, num_elec
       call zero_orbital_path(the_path)
       call init_tree_search(G, 0, 0, num_internal, i)
       do while ( get_internal_path(the_path, G)) 
          call analyze_CSF(wavefunc, the_path, table)
       end do
    enddo

    ! check to see if the file I'm writing these data structs to is open
    inquire( unit=io_output_buff, opened=opened_flag)
    if ( .not. opened_flag) open(unit=io_output_buff, file= scratch_directory // "output_buff", form='unformatted')
    rewind(io_output_buff)

    ! Instead of producing output here I'm temporarily storing these data on disk
    ! that way I can produce the output in a clean and orderly fashion later on
    write(io_output_buff) size(table%top_coeff)
    write(io_output_buff) size(table%ref_coeff)
    write(io_output_buff) table%num_ref_csf
    write(io_output_buff) table%top_coeff
    write(io_output_buff) table%top_csf
    write(io_output_buff) table%ref_coeff
    write(io_output_buff) table%ref_csf

    ! deallocate the path
    call deallocate_orbital_path(the_path)

  end subroutine single_CSF_analysis

  !> \brief Writes out the post iteration analysis to the output file
  subroutine write_post_iteration_analysis(reference_analysis)
    ! inputs
    logical, optional, intent(in) :: reference_analysis
      
    ! local variables
    type(CSF_data_table) :: table 
    integer :: top_size, ref_size, num_ref_csf

    ! read in the analysis
    rewind(io_output_buff)
    read(io_output_buff) top_size
    read(io_output_bufF) ref_size
    read(io_output_buff) num_ref_csf

    ! allocate space to read in the data
    call allocate_data_table(table, top_size, ref_size)

    ! read in the rest of the table 
    read(io_output_buff) table%top_coeff
    read(io_output_buff) table%top_csf
    read(io_output_buff) table%ref_coeff
    read(io_output_buff) table%ref_csf
  
    ! produce the output 
    call write_wav_analysis(table%top_coeff, table%top_csf, "  Best guess wavefunction summary")
    if ( present(reference_analysis) ) then 
        if (reference_analysis) then 
            call write_wav_analysis(table%ref_coeff, table%ref_csf, "  Weight of all the reference CSFs")
        end if
    end if
            
    call deallocate_data_table(table)
  end subroutine write_post_iteration_analysis
  
!> \brief Writes out a neat little wavefunction analysis 
!> \param coeff  The coefficients of each csf I want to write out
!> \param configs  Some data on each csf I want to write out
!> \param title    A title for the analysis
subroutine write_wav_analysis(coeff, configs, title)
    implicit none
    
    ! inputs
    real(kind=real8), dimension(:) :: coeff
    integer, dimension(:,:) :: configs
    character(len=*) :: title

    ! local variables
    integer :: i,j,k
    integer :: n, count
    integer :: ierr, num_rows
    integer :: start, finish
    integer, parameter :: MAX_PER_ROW = 50
    integer, dimension(:), allocatable :: index    

    real(kind=real8) :: c
    integer :: weight, spin, virt_a, virt_b
    character(len=100) :: row_format, internal_space_format, virt_format, extra_format
    character(len=1024) :: internal_space

    ! figure out the size 
    n = size(coeff)
    num_rows = num_internal/MAX_PER_ROW + 1
    if (mod(num_internal,MAX_PER_ROW) == 0) num_rows = num_rows - 1
    
    ! sort the results !!!
    allocate(index(n),stat=ierr)
    call allocatecheck(ierr, "index in write_post_iter_analysis")
    do i = 1, n
        index(i) = i
    end do
    call sort_real_array_with_index(coeff, index)
    
    100 format("   Coefficient  ", 3x," Spin func. ", 3x,"  Arc Weight   ", 3x," Internal Occ | Virtuals ",/,'')
    write(*,*) title
    write(*,*) "  --------------------------------"
    write(*,100) 

    write(internal_space,*) min(MAX_PER_ROW, num_internal)
    row_format ="(1x,ES15.5,6x,i5,8x,i7,12x)"
    internal_space_format = "(" // trim(internal_space) // "i1)"
    extra_format = "(55x," // trim(internal_space) // "i1)"
    virt_format = "(1x,'|',i4,1x,i4)"

    count = 0 
    do i = n, 1, -1
        ! skip over any entries that don't have anything (coeff = 0.0 )
        if ( coeff(i) == 0.0 ) cycle
        count = count + 1

        ! Get the values I need for this line
        weight = configs(index(i),num_internal+1)
        spin   = configs(index(i),num_internal+2)
        c  = coeff(i)
        virt_a = configs(index(i),num_internal+3)
        virt_b = configs(index(i),num_internal+4)
        
        if ( num_rows == 1 ) then 
            ! simple case 
            write(*,row_format,advance='no') c, spin, weight
            write(*,internal_space_format,advance='no') (configs(index(i),k), k= 1,num_internal)
            write(*,virt_format) virt_a, virt_b
            write(*,*) 
        else
            ! this is the more complicated case
            start = 1 
            finish = min(MAX_PER_ROW, num_internal)

            write(*,row_format,advance='no') c, spin, weight
            write(*,internal_space_format) (configs(index(i),k), k= start,finish)
            do j = 2, num_rows
                start = finish + 1
                finish = min(start + MAX_PER_ROW-1, num_internal)
                write(*, extra_format,advance='no') (configs(index(i),k), k= start,finish)
                if (j /= num_rows) write(*,*) 
            end do
            write(*,virt_format) virt_a, virt_b
            write(*,*) 
        end if
    enddo     
    
    ! sometimes no references contribute at all (single ref - excited states) 
    ! produce something (otherwise I might mistake it for a bug, again)
    if (count == 0 ) write(*,*) "     No Configurations Contributed    "
   deallocate(index)
end subroutine write_wav_analysis



  subroutine analyze_CSF(wavefunc, the_path, table)
    ! inputs
    type(orbital_path)::the_path
    real(kind=real8), dimension(:) :: wavefunc
    type(CSF_data_table) :: table

    ! local variables
    integer::path_weight, path_address
    integer::path_address_aa, path_address_ab
    integer::start, spin_dim
    integer::i,j,k,jk_count, path_elec,p,q
    integer::size_virt
    integer, dimension(1)::min_location
    integer, dimension(:), allocatable :: allowed_virtuals
    real(real8)::min_coefficient

    ! START AUTOGENERATED INITIALIZATION 
    start = 0
    path_weight = 0
    path_address = 0
    k = 0
    i = 0
    j = 0
    q = 0
    path_address_ab = 0
    jk_count = 0
    min_coefficient = 0.0
    min_location = 0
    spin_dim = 0
    p = 0
    path_address_aa = 0
    size_virt = 0
    path_elec = 0
    ! END AUTOGENERATED INITIALIZATION 


!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    path_elec = sum(the_path%occupations(0:num_internal))
    path_weight = sum(the_path%arc_weights(0:num_internal))+1

    path_address = 0
    path_address_aa = 0
    path_address_ab = 0  
    if (path_elec == num_elec) then
       path_address = internal_index_vector0(path_weight)
    elseif (path_elec == num_elec-1) then
       path_address = internal_index_vector1(path_weight)      
    elseif (path_elec == num_elec-2) then
       path_address_ab = internal_index_vector2(path_weight)
       path_address_aa = internal_index_vector3(path_weight)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !// JUST LEAVE IF THIS ONE WAS THROWN OUT
    if (path_address < 0 .or. &
         path_address_aa < 0 .or.&
         path_address_ab < 0) return

    if (path_elec == num_elec) then

       spin_dim = fsn(v_singles(path_weight))
       do i = 1, spin_dim
          start = path_address + (i-1)
          j = start 
          !!check to see if it matches a reference
          if( is_reference(the_path%occupations(1:num_internal)) ) then 
             table%num_ref_csf = table%num_ref_csf+1
             if (table%num_ref_csf > size(table%ref_csf,1)) then      !size of the first dimension
                call grow_ref_arrays(table)  
             end if
             table%ref_csf(table%num_ref_csf,num_internal+1) = path_weight
             table%ref_csf(table%num_ref_csf,num_internal+2) = i
             table%ref_csf(table%num_ref_csf,num_internal+3) = 0
             table%ref_csf(table%num_ref_csf,num_internal+4) = 0
             table%ref_csf(table%num_ref_csf,1:num_internal) = the_path%occupations(1:num_internal)    
             table%ref_coeff(table%num_ref_csf) = abs(wavefunc(j)) 
          end if

          min_location = minloc(table%top_coeff)
          min_coefficient = minval(table%top_coeff)
          if (abs(wavefunc(j)) > min_coefficient) then
             table%top_coeff(min_location(1)) = abs(wavefunc(j))
             table%top_csf(min_location(1),num_internal+1) = path_weight
             table%top_csf(min_location(1),num_internal+2) = i
             table%top_csf(min_location(1),num_internal+3) = 0
             table%top_csf(min_location(1),num_internal+4) = 0
             table%top_csf(min_location(1),1:num_internal) = the_path%occupations(1:num_internal)
          endif
       enddo

    elseif (path_elec == num_elec-1) then

       size_virt = num_allowed_virtuals(path_weight,"S")
       allocate(allowed_virtuals(size_virt))
       allowed_virtuals = 0
       call  get_virtuals( path_weight,  "S", allowed_virtuals ) ! a list of virtual orbitals we excite to


       spin_dim = fsn(nm1_singles(path_weight)+1)
       do i = 1, spin_dim
          start = path_address + (i-1)*size_virt
          jk_count= 0
          do j = 1, size_virt
             jk_count = jk_count+1
             min_location = minloc(table%top_coeff)
             min_coefficient = minval(table%top_coeff)
             if (abs(wavefunc(start + jk_count - 1)) > min_coefficient) then
                table%top_coeff(min_location(1)) = abs(wavefunc(start + jk_count - 1))
                table%top_csf(min_location(1),num_internal+1) = path_weight
                table%top_csf(min_location(1),num_internal+2) = i
                table%top_csf(min_location(1),num_internal+3) = allowed_virtuals(j)
                table%top_csf(min_location(1),num_internal+4) = 0
                table%top_csf(min_location(1),1:num_internal) = the_path%occupations(1:num_internal)
             endif
          enddo
       enddo

    elseif (path_elec == num_elec-2) then

       size_virt = num_allowed_virtuals(path_weight,"D")
       allocate(allowed_virtuals(size_virt))
       allowed_virtuals = 0
       call get_virtuals( path_weight,  "D", allowed_virtuals ) ! a list of virtual orbitals we excite to

       !// ONE VIRTUAL DOUBLY OCCUPIED    
       spin_dim = fsn(nm2_singles(path_weight))
       do i = 1, spin_dim
          start = path_address_aa + (i-1)*size_virt
          jk_count = 0
          do j = 1, size_virt
             jk_count = jk_count + 1
             min_location = minloc(table%top_coeff)
             min_coefficient = minval(table%top_coeff)
             if (abs(wavefunc(start + jk_count - 1)) > min_coefficient) then
                table%top_coeff(min_location(1)) = abs(wavefunc(start + jk_count - 1))
                table%top_csf(min_location(1),num_internal+1) = path_weight
                table%top_csf(min_location(1),num_internal+2) = i
                table%top_csf(min_location(1),num_internal+3) = allowed_virtuals(j)
                table%top_csf(min_location(1),num_internal+4) = 0
                table%top_csf(min_location(1),1:num_internal) = the_path%occupations(1:num_internal)
             endif
          enddo
       enddo

       !// TWO VIRTUALS SINGLY OCCUPIED
       spin_dim = fsn(nm2_singles(path_weight)+2)
       do i = 1, spin_dim
          start = path_address_ab + (i-1)*size_virt*(size_virt-1)/2
          jk_count = 0
          do j = 1, size_virt
             do k = 1, j-1
                jk_count = jk_count + 1
                min_location = minloc(table%top_coeff)
                min_coefficient = minval(table%top_coeff)
                if (abs(wavefunc(start + jk_count - 1)) > min_coefficient) then
                   table%top_coeff(min_location(1)) = abs(wavefunc(start+jk_count-1))
                   table%top_csf(min_location(1),num_internal+1) = path_weight
                   table%top_csf(min_location(1),num_internal+2) = i
                   table%top_csf(min_location(1),num_internal+3) = allowed_virtuals(j)
                   table%top_csf(min_location(1),num_internal+4) = allowed_virtuals(k)
                   table%top_csf(min_location(1),1:num_internal) = the_path%occupations(1:num_internal)
                endif
             enddo
          enddo
       enddo

    endif

  end subroutine analyze_CSF


  !> If we haven't allocated enough space for the reference arrays 
  !> this will reallocate them for us (dynamically!). Note we are 
  !> bound by the f90/95 standards and can't use dynamically allocated
  !> arrays as dummy arguments. Hence the arrays are global to the 
  !> module.
  subroutine grow_ref_arrays(table)
    ! input
    type(CSF_data_table) :: table

    ! local variables 
    integer :: old_size, new_size 
    integer, dimension(:,:), allocatable::new_ref_csf
    real(kind=real8), dimension(:), allocatable::new_ref_coeff
    integer :: ierr


    ! usual algorithm, double the current size 
    old_size = size(table%ref_coeff)
    new_size = 2 * old_size 

    ! allocate the new arrays
    allocate(new_ref_csf(new_size,num_internal+4), new_ref_coeff(new_size), stat=ierr)
    call allocatecheck(ierr, "memory for reallocation of the ref data in post it")
    new_ref_csf = 0
    new_ref_coeff = 0.0

    ! copy the data to new allocation
    new_ref_coeff(1:old_size) = table%ref_coeff(:)
    new_ref_csf(1:old_size,:) = table%ref_csf(:,:)

    ! point our data structure to the new memory
    call move_alloc(new_ref_coeff, table%ref_coeff)
    call move_alloc(new_ref_csf,   table%ref_csf)
  end subroutine grow_ref_arrays

end module post_iteration_mod












