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
!**************************************************************
!>  \brief Tiger CI - A Cholesky Decomposed Local Symmetric Group Graphical Approach CI code
!>
!>  \author Derek Walter
!>  \author Andrew Szilva
!>  \author Arun Venkathanathan
!>  \author Jeremy Chwee
!>  \author David Krisiloff
!>  \author Leah Isseroff
!>
!>  THIS CODE WAS INITIALLY A RIPOFF OF THE CODE PRESENTED BY DUCH AND KARWOWSKI.
!>  INDEED, I LOOKED AT WLODEK DUCH'S CODE WHEN WRITING THIS CODE.
!>  ACTUALLY, I HAVE ADDED THE VECTORIZATION OF THE EXTERNAL 
!>  SPACE, WHICH IS A PRETTY MAJOR IMPROVEMENT.
!> 
!>  REFERENCE:
!>  W. DUCH; J. KARWOWSKI COMPUTER PHYSICS REPORTS 2 (1985) 93 - 170
!> 
!**************************************************************

!> \todo Need to find all the arrays still placed on the stack and move them to the heap
!> \todo Try and remove any remaining matmul calls of the form A= MATMUL(A,B) just in case

subroutine tiger_ci

  !// MAIN PROGRAM
  !// CALLS TO:
  !//    TIMESTAMP
  !//    CI_INPUT
  !//    INTEGRAL_SORT
  !//    GRAPH_DRIVER
  !//    SPIN
  !//    DIAGONALIZE


  use global_var_mod       !// GLOBAL VARIABLES
  use molecule_var_mod     !// MOLECULAR INFORMATION
  use utilities_mod        !// UTILITIES      
  use graph_mod            !// THE GRAPH MODULE
  use graph_var_mod        !// GRAPH VARIABLES
  use spin_mod             !// FOR THE SPIN INTEGRALS
  use natural_mod          !// FOR COMPUTING NATURAL ORBITALS
  use locist_mod
  use cholesky_driver_mod
  use molecular_orbital_mod
  use one_electron_integrals
  use sphere_cylinder 
  use wp_tov_mod
  use make_all_integrals
  use diag_element_mod
  use pseudo_seg_info_mod
  use get_integrals_mod
  use sgga_struct_mod
  use sdci_driver_mod
  use acpf_driver_mod
  use fortran_timing
  use get_basis_data
  use two_electron_integrals
  use density_fitting_interface
  use IOBuffer
  use number_virt_per_hole
  
  implicit none

  integer::i
  integer::deallocatestatus

  real(real8)::zero     !// 0.0
  real(real8)::nuclear_energy
  real(real8),allocatable::energy(:)
  
  type(clock)::timer
  type(sgga_struct) :: sgga_data
  zero = real(0.0,real8)
  !// OPEN INPUT AND OUTPUT FILES.  ADDITIONAL SCRATCH FILES ARE OPENED
  !// BELOW IN OPEN_FILES.  WE NEED THIS FOR NOW, HOWEVER, SO THAT
  !// WE CAN WRITE ERRORS OUT AND GET INFORMATION ABOUT THE SCRATCH
  !// FILES WE PLAN TO OPEN.

  call  flush(ioOutput)

  !// NOW WE'LL GET THE RELEVANT INPUT DATA.  NUMBER OF ELECTRONS, ORBITALS,
  !// SPIN, REFERENCE INFO, ETC.  THIS IS WHERE WE DO IT.  
  !call ci_input_read

  !// OPEN ALL FILES BUT INPUT,OUTPUT, AND ERROR FILES, WHICH WERE OPENED ABOVE
  call open_files

  !// TRIAL TO GET THE NUCLEAR REPULSION ENERGY
  nuclear_energy = zero
  call get_nuclear_repulsion(nuclear_energy)
  call flush(ioOutput)
  
  ! New interface to the MOs
  call get_molecular_orbitals()

  call get_num_atoms(num_atoms)
  call make_spheres_cylinders 
  call build_ignorable_pair_matrix

  ! Two-electron integrals
  write(*,*) "IN FORTRAN, CPP_DECOMPOSED_INTS HAS A VALUE OF ", CPP_DECOMPOSED_INTS
  write(*,*) "        AND CPP_TRANSFORMED_INTS HAS A VALUE OF ", CPP_TRANSFORMED_INTS
  if (CPP_TRANSFORMED_INTS) then 
      call setup_density_fitting(sgga_data%cho_data)
  else
    call cholesky_driver(sgga_data%cho_data,restart_from_old_CD, max_mem)
  end if
    
!write(*,*) "Done with cholesky"
!call flush(6)
!stop

#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: Entering pppp preparation"
  call flush(6)
#endif
  call getpppp_cho(sgga_data%cho_data)
#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: Entering prpr preparation"
  call flush(6)
#endif
  call getprpr_cho(sgga_data%cho_data)
#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: Entering pprr preparation"
  call flush(6)
#endif
  call getpprr_cho(sgga_data%cho_data)
#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: Finished pprr preparation"
  call flush(6)
#endif
 
  ! One-electron integral handling
  call one_electron_driver(molecular_orbitals)

  write(ioOutput,*) "Creating integral buffer of maximum size ",for_buf_maxmemInts
  write(ioOutput,*) "                   secondary AIJK buffer ",for_buf_maxmemAIJK
  write(*,*) "Number of threads: ",numberOfThreads
  flush(6)
  call for_double_buf_construct(for_buf_maxmemInts/(for_buf_blocksizeInts*real8),for_buf_blocksizeInts,0,numberOfThreads,for_buf_int_poolID)
  call for_int_buf_construct(for_buf_maxmemInts/(for_buf_blocksizeInts*real8),for_buf_blocksizeInts,0,numberOfThreads,for_buf_integer_int_poolID) 
 
!  if(.not. integralDirect .or. sphere_based_integral_truncations) then
    call for_double_buf_construct(for_buf_maxmemAIJK/(num_external*real8),num_external,0,numberOfThreads,for_buf_aijk_poolID)
    call for_double_buf_openfile(for_buf_aijk_poolID,ijka_test_no,&
       scratch_directory // 'ijka_test.dat',len(scratch_directory) + 13) ! open the file here to not run into problems with rebooting SGGA
!  endif
  write(ioOutput,*)
  write(ioOutput,*) "* Entering initial integral preparation"
  call flush(ioOutput)
  call start_clock(timer)
  call makeIntegrals(sgga_data%cho_data)
  call print_clock(timer, "  + integral time: ")
  write(ioOutput,*) " => Did integrals"
  call flush(ioOutput)
  if(for_buf_storeIntegrals) then
    write(ioOutput,*) "Storing buffered integrals to disk"
    call for_double_buf_syncpool(for_buf_int_poolID)
  endif
  
  ! a couple of files and pseudofiles
  call for_int_buf_openfile(for_buf_integer_int_poolID,twosingles_no,scratch_directory // '2s2singles.dat',len(scratch_directory) + 14)
  call setupPseudoInfoFiles()
  
  ! allocate out loc_scr
  call allocLocScratch(sgga_data%loc_scr,num_external)

  !// LETS MAKE THE GRAPH.  ONLY SINGLE AND DOUBLE
  !// EXCITATIONS WILL BE ALLOWED.  ALSO, WE WILL BUILD THE FULL GRAPH (INTERNAL AND
  !// EXTERNAL SPACE).  LATER WE WILL SEPARATE THE EXTERNAL AND INTERNAL SPACE FULLY.
  !// CALLING THIS ROUTINE WILL ALSO GENERATE THE INDEX VECTOR AND THE ARRAY FSN WHICH
  !// STORES THE NUMBER OF SPIN FUNCTIONS FOR A GIVEN NUMBER OF OPEN SHELLS.
  call start_clock(timer)
  call graph_driver(sgga_data%loc_scr)
  graph_time = get_clock_wall_time(timer)

  write(6,*) "Done with calling graph_driver"
  call count_virts_per_hole
  call flush(6)

  !// NOW WE ARE GOING TO DO THE COUPLING COEFFICIENTS (I.E. REPRESENTATION MATRICES OF
  !// THE SYMMETRIC GROUP).  WE GENERATE CYCLES AND TRANSPOSITIONS HERE.  THE COUPLING
  !// COEFFICIENTS WILL BE BUILT UP ON THE FLY FROM THESE MATRICES.
  call start_clock(timer)
  call spin
  spin_time = get_clock_wall_time(timer)

  write(6,*) "Done with calling spin"
  call flush(6)

  !// FLUSH THE OUTPUT FOR DEBUGGING
  call flush(ioOutput)
  
  ! Before we start the eigensolvers we construct the diagonal elements of the Hamiltonian and write them to disk
  ! Agreed this seems oddly specific to place in the main program but removing it from the solvers simplifies a lot of things ... so oh well
  call start_clock(timer)
  call diag_element(sgga_data%cho_data,sgga_data%loc_scr)
  diag_element_time = get_clock_wall_time(timer)
  rewind(iodiagstore)
  call write_big_vector(diagonal_elements,total_csfs,iodiagstore)
  deallocate(diagonal_elements, stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"diagonal")
  
#ifdef TIGER_USE_OMP
  ! the record stuff for threeiseg
  call for_int_buf_construct(for_buf_threeiseg_maxmem/(for_buf_threeiseg_block*real8),for_buf_threeiseg_block,0,numberOfThreads,for_buf_threeiseg_poolID)
  call for_int_buf_openfile(for_buf_threeiseg_poolID,threeiseg_no,scratch_directory // 'threeiseg_omp.dat',len(scratch_directory) + 17)
#endif
  
  !// NOW WE ARE READY TO DO THE DIAGONALIZATION STEP
  call start_clock(timer)
  allocate(energy(num_roots))
  if (acpf_flag  /= 0 ) then 
      call acpf_driver(num_roots, energy, sgga_data)
  else
      call sdci_driver(num_roots, energy, sgga_data)
  end if
  diag_tot_time = get_clock_wall_time(timer)
  call print_clock(timer, "eigenvalue solver")




  !// WRITE TO THE OUTPUT FILE THE TIMES
  write(ioOutput,400) sort_time, graph_time, spin_time, diag_tot_time,&
       diagtimesc_time,two_seg_time, &
       three_and_four_seg_time,ijkl_time,aijk_time,&
       abij_time,ikjk_time,abcd_time,iabc_time,ijaj_time

400 format(1x,/,1x,"Timing summary (sec)",/,&
       1x,"-----------------------------------------",/,&
       1x,"Integral sorting time.....................",f10.3,/,&
       1x,"Graph generation time.....................",f10.3,/,&
       1x,"Spin integral generation time.............",f10.3,/,&
       1x,"Total diagonalization time................",f10.3,/,&
       6x,     "* diagonal elements time.............",f10.3,/,&
       6x,     "* two segment loop time..............",f10.3,/,&
       6x,     "* three and four segment loop time...",f10.3,/,&
       8x,       "+ ijkl time........................",f10.3,/,&
       8x,       "+ aijk time........................",f10.3,/,&
       8x,       "+ abij time........................",f10.3,/,&
       8x,       "+ ikjk time........................",f10.3,/,&
       8x,       "+ abcd time........................",f10.3,/,&
       8x,       "+ iabc time........................",f10.3,/,&
       8x,       "+ ijaj time........................",f10.3,/)

  do i = 1, num_roots
      write(ioOutput,430) i,nuclear_energy+energy(i)
  430 format(/,1x,"Root " ,i5 ," total energy (electronic + nuclear) is : ",f24.15,/)
  end do

  !**********************************************************************************
  deallocate(ignorable_pair,energy,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"ignorable_pair  ")

#ifdef DEBUG_TIGER
  write(ioOutput,*) "DEBUG: Finished cleaning up arrays."
  write(ioOutput,*) "DEBUG: Cleaning up local scratch"
  flush(ioOutput)
#endif
  ! allocate out loc_scr
  call deallocLocScratch(sgga_data%loc_scr)

#ifdef DEBUG_TIGER
  write(ioOutput,*) "DEBUG: Taking care of buffers."
  flush(ioOutput)
#endif
  
  call deletePseudoInfoFile()
  
  ! close our integral buffer
  if(for_buf_storeIntegrals) then
#ifdef DEBUG_TIGER
  write(ioOutput,*) "DEBUG: Trying to close integral pools."
  flush(ioOutput)
#endif
    call for_double_buf_closepool(for_buf_int_poolID)
    if(.not. integralDirect .or. sphere_based_integral_truncations) then
      call for_double_buf_closepool(for_buf_aijk_poolID)
    endif
  else
#ifdef DEBUG_TIGER
  write(ioOutput,*) "DEBUG: Trying to remove integral pools."
  flush(ioOutput)
#endif
    call for_double_buf_removepool(for_buf_int_poolID)
    if(.not. integralDirect .or. sphere_based_integral_truncations) then
      call for_double_buf_removepool(for_buf_aijk_poolID)
    endif
  endif
  
!#ifdef DEBUG_TIGER
  write(ioOutput,*) "DEBUG: Took care of buffers, closing remaining files."
  flush(ioOutput)
!#endif
  
  ! close some leftover files
  close(53)
  close(72)
  close(75)
  close(76)
  close(77)
  close(80)
  close(96)

end subroutine tiger_ci
