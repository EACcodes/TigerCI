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
!//  GRAPH_MOD - THIS MODULE CONTAINS THE ROUTINES FOR CONSTRUCTING
!//  THE GRAPH.  ALSO, WE CONSTRUCT THE INDEX VECTOR AND  THE ARRAY FSN
!//  WHICH STORES THE NUMER OF SPIN FUNCTIONS FOR A GIVEN NUMBER OF OPEN 
!//  SHELLS.
!// 
!//  WRITTEN BY DEREK WALTER, 1999
!//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************

module graph_mod

  use io_unit_numbers

contains

  subroutine graph_driver(loc_scr)

    !// THIS SUBROUTINE IS THE DRIVER SUBROUTINE FOR THE CONSTUCTION OF
    !// THE GRAPH.  WE FOLLOW THE PRESCRIPTIONS GIVEN BY DUCH AND KARWOWSKI
    !// CALLED BY:           CALLS TO:
    !//    MR_CI                GET_FSN
    !//                         GRAPH_BORDER
    !//                         DRT_MAKER
    !//                         TREE_SEARCH_INTERNAL
    !//                         INDEX_VECTOR_INTERNAL
    use global_var_mod
    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use three_four_seg_var_mod
    use two_seg_var_mod

    implicit none
    type(locist_scratch)::loc_scr

    write(ioOutput,130)
130 format(/,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/&
         ,1x,"!//                              ",/&
         ,1x,"!// ENTERED DRT GENERATION       ",/&
         ,1X,"!// ROUTINE                      ",/&
         ,1x,"!//                              ",/&
         ,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/)

   
    !// FIRST, WE NEED TO CALL THIS SILLY ROUTINE TO GENERATE THE FSN ARRAY
    !// WHICH IS THE ARRAY OF THE NUMBER OF SPIN FUNCTIONS ASSOCIATED WITH
    !// A GIVEN NUMBER OF OPEN SHELLS AND A GIVEN MULTIPLICITY
    call get_fsn

    !// SECOND, WE NEED TO MAKE THE BORDERS OF THE GRAPH AND THE AUXILIARY VECTOR  
    call graph_border

    !// THIRD, WE GENERATE THE DRT (FULL SPACE, NOT JUST INTERNAL)
    call drt_maker

    !// FOURTH, WE HAVE TO DO THE TREE SEARCH ALGORITHM TO GET ALL THE PATHS
    call tree_search_internal

    !// FORM THE INDEX VECTORS; THIS IS A SUPER CRITICAL ROUTINE;
    !// ALL OF THE CONFIGURATION SELECTION YOU MIGHT WANT TO DO
    !// FOR LOCAL C.I. ETC. TAKES PLACE IN THIS ROUTINE, I.E.
    !// LOCAL_DRIVER AND ALLOCATE_PAO_QUANTITIES ARE CALLED FROM 
    !// HERE.  
    call index_vector_internal(loc_scr)


  end subroutine graph_driver

  !**************************************************************
  subroutine graph_border

    !// THIS SUBROUTINE CONSTRUCTS THE BORDERS OF THE GRAPH. 
    !// WE DO THE LEFT MOST BORDER, THE RIGHT MOST BORDER, AND
    !// THEN WE TRIM THE BORDERS DOWN TO MAKE SURE THEY ARE CONSISTENT
    !// WITH SINGLE AND DOUBLE EXCITATION FROM THE REFERENCES.  
    !// CALLED BY:
    !//    GRAPH_DRIVER

    use global_var_mod              !// GLOBAL VARIABLES
    use utilities_mod               !// VARIOUS UTILITIES
    use molecule_var_mod            !// INFO ABOUT THE MOLECULE
    use graph_var_mod               !// MODULE TO STORE GRAPH VARIABLES

    implicit none

    integer::i,j                           !// LOOP CONTROL VARIABLES
    integer::allocatestatus                !// STATUS VARIABLE FOR DYNAMIC MEMORY ALLOCATION
    integer::min_singles                   !// MINIMUM NUMBER OF SINGLES REQUIRED TO ACHIEVE DESIRED SPIN MUTIPLICITY
    integer::right_doubly, left_doubly     !// NUMBER OF DOUBLY OCCUPIED ORBITALS IN THE RIGHT AND LEFT PATHS
    integer::left_unocc                    !// NUMBER OF UNOCCUPIED ORBITALS IN LEFT BORDER
    integer::right_occ                     !// NUMBER OF OCCUPIED ORBITALS IN RIGT BORDER
    integer::left_singles, right_singles   !// NUMBER OF SINGLES IN LEFT, RIGHT BORDER BORDER
    integer::add                           !// NUMBER OF ELECTRONS WE ADD AT EACH STEP
    integer::min_elec, max_elec            !// MINIMUM AND MAXIMUM CUMMULATIVE OCCUPANCIES AT EACH LEVEL
    !//    CONSISTENT WITH SINGLES AND DOUBLES CI  
    integer::cummulative                   !// CUMULATIVE OCCUPANCY AT EACH LEVEL FOR A REFERENCE   
    integer,dimension(:),allocatable::temp      !// TEMPORARY WORKSPACE ARRAY

    min_singles = spinM - 1

    !// FOR THE RIGHTMOST BORDER WE NEED THE FOLLOWING INFO
    right_singles = min_singles  
    right_doubly = (num_elec - right_singles)/2
    right_occ = min_singles + right_doubly 

    !// FOR THE LEFT BORDER WE WILL NEED THE FOLLOWING
    left_singles = max(0,min_singles)
    left_doubly = (num_elec - left_singles)/2
    left_unocc = num_orbitals - left_singles - left_doubly

    !// ALLOCATE SOME ARRAYS TO STORE LEFT AND RIGHT BORDERS
    allocate(left_border(0:num_orbitals), stat = allocatestatus)
    call allocatecheck(allocatestatus,"left_bor")

    allocate(right_border(0:num_orbitals), stat = allocatestatus)
    call allocatecheck(allocatestatus,"right_bo")

    !// INITIALIZE OUR ARRAYS
    right_border = 0
    left_border = 0

    !// NOW LETS DO THE LEFT MOST BORDER.  IN THE ARRAYS WE ARE GOING TO
    !// STORE THE CUMULATIVE OCCUPATIONS
    add = 0
    do i = 1,num_orbitals  
       if (i <= left_unocc) then
          add =  0
       elseif ((i > left_unocc).and.(i <= (left_unocc + left_singles))) then
          add = 1
       elseif (i > (left_unocc + left_singles)) then
          add = 2
       endif

       if (i == 1) then
          left_border(i) = add
       else
          left_border(i) = left_border(i-1) + add
       endif
    enddo

    !// NOW THE RIGHT BORDER
    add = 0
    do i = 1,num_orbitals
       if (i <= right_doubly) then
          add = 2
       elseif ((i > right_doubly).and.(i <= (right_doubly + right_singles))) then
          add = 1
       elseif(i > (right_doubly + right_singles)) then
          add = 0
       endif

       if (i == 1) then
          right_border(i) = add
       else
          right_border(i) = right_border(i-1) + add
       endif
    enddo

    !// NOW WE HAVE TO CHECK AND MAKE SURE THESE BORDERS ARE CONSISTENT WITH SINGLE
    !// AND DOUBLE EXCITATIONS FROM THE REFERENCES.  WE'LL CHECK IT FOR EACH REFERENCE.

    !// WE'LL NEED A TEMPORARY ARRAY FOR THIS TASK
    allocate(temp(num_ref), stat = allocatestatus)
    call allocatecheck(allocatestatus, "temp    ")

    !// INITIALIZE TEMPORARY ARRAY
    temp(1:num_ref:1) = references(1,1:num_ref:1)
    orbitals: do i = 2, num_orbitals
       min_elec = 2*i   !// INITIALIZE WITH LARGEST POSSIBLE VALUE, AND SEE IF WE CAN FIND A SMALLER ONE
       max_elec = 0     !// INITIALIZE WITH SMALLEST POSSIBLE VALUE, AND SEE IF WE CAN FIND A LARGER ONE

       reference_loop: do j = 1, num_ref

          !// FOR EACH REFERENCE ADD UP THE TOTAL NUMBER OF ELECTRONS ORBITAL WISE
          cummulative = temp(j) + references(i,j)

          !// AT THIS PARTICULAR LEVEL (THE I VALUE), FIND THE SMALLEST AND LARGEST
          !// OCCUPATION IN ANY OF THE REFERENCES.
          if (cummulative > max_elec) max_elec = cummulative
          if (cummulative < min_elec) min_elec = cummulative

          !// STORE THE CUMMULATIVE OCCUPATION PATTERN FOR THIS REFERENCE
          temp(j) = cummulative
       enddo reference_loop

       !// NOW READJUST THE BORDERS SO THAT THEY ARE CONSISTENT WITH SINGLES
       !// AND DOUBLES CI
       left_border(i) = max(min_elec - exilevel,left_border(i))
       right_border(i) = min(max_elec + exilevel, right_border(i))

    enddo orbitals
    deallocate(temp)

    !// NOW WE'LL GET THE AUXILIARY ARRAY.  WE WANT THIS TO BE AN ARRAY SUCH THAT
    !//  AUXILIARY(I) + J GIVES THE NUMBER OF THE VERTEX WITH J ELECTRONS IN ROW I.   
    !// FIRST ALLOCATE THE ARRAY.
    allocate(auxiliary(0:num_orbitals), stat = allocatestatus)
    call allocatecheck(allocatestatus,"auxiliar")

    num_vert = 1
    auxiliary(0) = 1
    do i = 1, num_orbitals
       auxiliary(i) = num_vert + 1 - left_border(i)
       num_vert = right_border(i) + auxiliary(i)
       !     auxiliary(i) = right_border(i-1) - left_border(i) + auxiliary(i-1) + 1      
    enddo

    if (borderprint > 0) then
       write(ioOutput,*) "Graph border and auxiliary vector finished..... "
       write(ioOutput,*) "Graph border information: "
       write(ioOutput,200) 
200    format (1x," level  ", " left border  "," right border ",/,&
            1x,"------- ", "------------- ","--------------")
       do i = 1, num_orbitals
          write(ioOutput,210) i, left_border(i), right_border(i)
       enddo
210    format (1x,i4,6x,i6,8x,i6)
       write(ioOutput,220) num_vert
220    format(1x,/,&
            1x,"Total number of vertices.....",i5,/)
    endif
    if (auxprint > 0) then
       write(ioOutput,*)"Auxiliary vector:"
       write(ioOutput,*) (auxiliary(j), j = 1,num_orbitals) 
       write(ioOutput,*)
       call flush(ioOutput)
    endif

  end subroutine graph_border

  !**************************************************************
  subroutine drt_maker

    !// THIS SUBROUTINE GENERATES THE DRT USING THE BORDER INFORMATION AND AUXILIARY VECTOR
    !// THAT WE OBTAINED FROM GRAPH_BORDER.  I'LL DO THIS MORE OR LESS THE WAY DUCH DOES THIS,
    !// BUT I'LL ELIMINATE THE GOTOS AND MAKE IT MUCH PRETTIER  
    !// CALLED BY:
    !//    GRAPH_DRIVER

    use global_var_mod
    use graph_var_mod
    use molecule_var_mod
    use utilities_mod
    use ci_utilities_mod

    implicit none

    integer::allocatestatus            !// DYNAMIC MEMORY ALLOCATION STATUS VARIABLE
    integer::i,j                       !// LOOP CONTROL VARIABLES
    integer::vertex_above              !// VERTEX ABOVE THE CURRENT VERTEX
    integer::vert                      !// THE CURRENT VERTEX
    integer::weight                    !// STORES ARC WEIGHTS AND VERTEX WEIGHTS

    !// FIRST THING IS TO ALLOCATE THE Y0, Y1 AND Y2 ARRAYS
    allocate(y0(num_vert), stat = allocatestatus)
    call allocatecheck(allocatestatus, "y0      ")
    allocate(y1(num_vert), stat = allocatestatus)
    call allocatecheck(allocatestatus, "y1      ")
    allocate(y2(num_vert), stat = allocatestatus)
    call allocatecheck(allocatestatus, "y2      ")

    !// INITIALIZE THE GRAPH.  CHECK TO SEE IF THE RELEVANT VERTICES ARE CONTAINED WITHIN THE GRAPH.  IF 
    !// THEY ARE, THEN SET Y0 TO 1 AND THE Y1 AND Y2 TO 0.  
    if (belong(1,0)) then
       y0(1) = 1
    else
       y0(1) = -1
    endif

    if (belong(1,1)) then
       y1(1) = 0 
    else
       y1(1) = -1
    endif

    if (belong(1,2)) then
       y2(1) = 0
    else
       y2(1) = -1
    endif

    !// WE HAVE THINGS INITIALIZED.  NOW, WE WILL LOOP OVER LEVELS AND FOR EACH LEVEL WE WILL 
    !// LOOP OVER VERTICES, BUILDING UP THE GRAPH ALONG THE WAY.  
    orbitals: do i = 1, num_orbitals

       !WE ARE GOING TO LOOP OVER THE VERTICES FROM RIGHT TO LEFT!!
       vertices: do j = right_border(i), left_border(i), -1

          !// THIS IS THE VERTEX DIRECTLY ABOVE THE CURRENT VERTEX
          vertex_above = vertex(i-1,j)  

          !// AND THIS IS THE CURRENT VERTEX
          vert = vertex(i,j)            

          !// FIRST CALCULATE THE WEIGHT OF THE VERTEX.  WE ARE GOING TO CHECK AND
          !// SEE IF EACH VERTEX THAT CAN CONNECT TO THE CURRENT VERTEX IS THERE.  
          !// IF IT IS, WE ADD IN ITS CONTRIBUTION TO THE WEIGHT.
          weight = 0
          if (add0(i-1,j))   weight = weight + abs(y0(vertex_above))
          if (add1(i-1,j-1)) weight = weight + abs(y0(vertex_above - 1))
          if (add2(i-1,j-2)) weight = weight + abs(y0(vertex_above - 2))
          if (.not.belong(i+1,j)) then
             y0(vert) = -weight   !// THE VERTEX DIRECTLY BELOW THIS ONE IS OUT OF THE GRAPH
          else
             y0(vert) = weight
          endif

          !// NOW WE NEED TO LOOK AT THE LEVELS BELOW IN ORDER TO CALCULATE THE ARC WEIGHTS
          !// FIRST WE WILL DO THE 1 ARCS.  CHECK TO SEE IF THE ARC IS THERE, IF IT IS THEN
          !// COMPUTE THE WEIGHT.
          weight = 0
          if (.not.belong(i+1,j+1)) then  !// THE ARC IS OUTSIDE THE GRAPH BORDER
             weight = -1 
          else
             if (add0(i,j+1)) weight = abs(y0(vertex(i,j+1)))  !// FOR 1 ARCS ONLY ONE WEIGHT CAN CONTRIBUTE TO
             !// THE ARC WEIGHT
          endif
          y1(vert) = weight

          !// NOW LETS DO THE 2 ARCS.  MORE OF THE SAME HERE.   
          weight = 0
          if (.not.belong(i+1,j+2)) then
             weight = -1
          else 
             if (add1(i,j+1)) weight = weight + abs(y0(vertex(i,j+1)))
             if (add0(i,j+2)) weight = weight + abs(y0(vertex(i,j+2)))
          endif
          y2(vert) = weight

       enddo vertices
    enddo orbitals

    !// O.K., IF ALL WENT WELL, WE SHOULD HAVE THE DRT.  LETS PRINT IT TO THE OUTPUT FILE
    if (drtprint > 0) then
       write(ioOutput,*) "Finished generating DRT....."
       write(ioOutput,*) "DRT:"
       write(ioOutput,220)
220    format (1x, "  vertex  ",10x,"    y0     ",10x,"     y1     ",10x,"     y2     ",/,&
            1x, "----------",10x,"-----------",10x,"------------",10x,"------------")
       do i = 1,num_vert
          write(ioOutput,230)i,y0(i),y1(i),y2(i)
       enddo
230    format(1x,i5,6x,i15,7x,i15,7x,i15)
       write(ioOutput,*)
    endif

  end subroutine drt_maker

  !**************************************************************
  subroutine get_fsn

    !// THIS SUBROUTINE DETERMINES THE NUMBER OF SPIN FUNCTIONS FOR A GIVEN
    !// NUMBER OF OPEN SHELLS AND THE GIVEN SPIN MULTIPLICITY.  THE METHOD THE
    !// ALGORITHM USES IS BASED ON SPIN BRANCHING DIAGRAMS.  IF YOU STARE AT ONE
    !// OF THESE ALGORITHMS LONG ENOUGH, YOU CAN FIGURE OUT EXACTLY WHAT IS
    !// GOING ON.
    !// CALLED BY:
    !//    GRAPH_DRIVER

    use global_var_mod
    use utilities_mod
    use molecule_var_mod
    use graph_var_mod

    implicit none

    integer::i,j            !// LOOP CONTROL VARIABLES
    integer::allocatestatus    !// DYNAMIC MEMORY ALLOCATION VARIABLE
    integer::num_functions     !// NUMBER OF SPIN FUNCTIONS
    integer :: ref_singles
    
    ! Adding the following logic for when we need to run the graph driver more than once
    ! (We need to recalculate the original number of open shells)
    ! this logic is copied right from ci_input.f90
      open_shells = 0
      do i = 1, num_ref

          ref_singles = 0
          do j = num_inactive+1, num_inactive+num_active
              if (references(j,i)==1) ref_singles = ref_singles+1
          enddo

          open_shells = max(ref_singles+4,open_shells)

      enddo
      if (num_elec < open_shells) open_shells = num_elec
    
    

    !// IF THE SPIN MULTIPLICITY IS AN EVEN MULTIPLE OF TWO, IT MEANS
    !// THAT WE HAVE AN ODD NUMBER OF OPEN SHELLS.  WE ARE GOING TO REDIMENSION
    !// OPEN_SHELLS, WHICH IS THE MAXIMUM NUMBER OF OPEN SHELLS, SO THAT IT IS EITHER
    !// ODD OR EVEN TO CONFORM TO THE PARTICULAR SITUATION.
    if ((mod(spinM,2) == 0).and.(mod(open_shells,2)== 0)) open_shells = open_shells + 1
    write(ioOutput,*) "Maximum number of open shells allowed: ",open_shells

    !// NOW WE CAN ALLOCATE OUR ARRAYS.  REMEMBER, THAT THE ARRAYS MUST RUN
    !// FROM 0 TO OPEN_SHELLS BECAUSE 0 IS CERTAINLY A VALID VALUE FOR THE 
    !// NUMBER OF SINGLES. 

    if (open_shells > fsn_size) then

       write(ioOutput,*) "You are requesting a calculation with ", open_shells,"open shells." 
       write(ioOutput,*) "This is alot.  I can do this many if you really want to."
       write(ioOutput,*) "Right now I am set to stop at ", fsn_size, " open shells."
       write(ioOutput,*) "If you really want to do more, you need to go into global_var_mod"
       write(ioOutput,*) "and change the variable fsn_size."
       stop

    endif


    allocate(fsn(0:fsn_size),stat = allocatestatus)
    call allocatecheck(allocatestatus, "fsn     ")
    allocate(fsn3(0:open_shells+1),stat = allocatestatus)
    call allocatecheck(allocatestatus, "fsn3    ")

    !// INITIALIZE THE ARRAYS 
    fsn = 0
    fsn3 = 0
    fsn(0) = 1
    fsn3(0) = 1
    fsn3(1) = 1

    !// WE NOW LOOP OVER ALL THE POSSIBLE VALUES FOR THE OPEN SHELLS
    do i = 1, open_shells

       !// FOR EACH VALUE OF THE NUMBER OF OPEN SHELLS, WE LOOP OVER THE ALLOWED VALUES OF THE SPIN MULTIPLICITY
       !// CONSISTENT WITH THE CURRENT NUMBER OF OPEN SHELLS
       do j = mod(i,2) + 1, i+1, 2 
          if (j == 1) then

             !// THE NUMBER OF ALLOWED SPIN FUNCTIONS EQUALS THE NUMBER OF SPIN FUNCTIONS WHICH
             !// WERE DOUBLETS AT THE PREVIOUS LEVEL, BECAUSE WE CAN COUPLE 1 ELECTRON TO MAKE
             !// THEM ALL SINGLETS.
             num_functions = fsn3(2)

          elseif (j == i+1) then

             !// NOW ALL OUR ELECTRONS ARE HIGH-SPIN COUPLED SO WE HAVE ONLY
             !// ONE ALLOWED SPIN FUNCTION
             num_functions = 1

          else
             !// NOW, WE CAN GENERATE THE DESIRED MULTIPLICITY TWO WAYS.
             !// WE CAN COUPLE AN ELECTRON WITH SPIN UP, OR WE CAN COUPLE AN
             !// ELECTRIN WITH SPIN DOWN.  THUS WE NEED TO ADD THE TWO VALUES FROM 
             !// THE PREVIOUS LEVEL
             num_functions = fsn3(j-1) + fsn3(j+1)

          endif

          !// NOW UPDATE FSN3 WITH THE APPROPRIATE NUMBER OF SPIN FUNCTIONS
          fsn3(j) = num_functions
       enddo

       fsn(i) = fsn3(spinM)
    enddo

  end subroutine get_fsn

!   !*****************************************************************
!   subroutine get_csfs_off_disk(iounit)

!     !// THIS ROUTINE READS CSFS OFF OF A LOGICAL UNIT.  ALSO ITS SETS
!     !// THE INTERNAL INDEX VECTORS IF NECESSARY

!     use graph_var_mod
!     use molecule_var_mod
!     use utilities_mod
!     use global_var_mod
!     use locist_mod

!     implicit none

! !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
!     integer::iounit                     !// UNIT WE ARE READING FROM
!     integer::i                          !// LOOP VARIABLES
!     integer::count                      !// FOR KEEPING TRACK OF INTERNAL CSFS THROW OUT
!     integer::allocatestatus             !// FOR DUNAMIC MEMORY ALLOCATION


! !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     if (.not.allocated(internal_index_vector0)) then
!        allocate(internal_index_vector0(abs(y0(vertex(num_internal,num_elec)))),stat=allocatestatus)
!        allocate(internal_index_vector1(abs(y0(vertex(num_internal,num_elec-1)))),stat=allocatestatus)
!        allocate(internal_index_vector2(abs(y0(vertex(num_internal,num_elec-2)))),stat=allocatestatus)
!        allocate(internal_index_vector3(abs(y0(vertex(num_internal,num_elec-2)))),stat=allocatestatus)
!        call allocatecheck(allocatestatus,"internaL")
!     endif

! !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     write(ioOutput,*)        
!     write(ioOutput,*) "Determining CSFs off disk ....."
!     write(ioOutput,*)

!     write(ioOutput,*) "ioUnit ", ioUnit

!     rewind(ioUnit)  
!     read(ioUnit) internal_index_vector0
!     read(ioUnit) internal_index_vector1
!     read(ioUnit) internal_index_vector2
!     read(ioUnit) internal_index_vector3

!     count = 0
!     do i =1, size(internal_index_vector0)
!        if (internal_index_vector0(i) == kill_number) then
!           count = count + 1
!        endif
!     enddo
!     write(ioOutput,*) "# Valence configuratons thrown out: ",count

!     count = 0
!     do i =1, size(internal_index_vector1)
!        if (internal_index_vector1(i) == kill_number) then
!           count = count + 1
!        endif
!     enddo
!     write(ioOutput,*) "# N-1 configurations thrown out: ",count

!     count = 0
!     do i =1, size(internal_index_vector2)
!        if (internal_index_vector2(i) == kill_number) then
!           count = count + 1
!        endif
!     enddo
!     write(ioOutput,*) "# N-2 AB configurations thrown out: ",count

!     count = 0
!     do i =1, size(internal_index_vector3)
!        if (internal_index_vector3(i) == kill_number) then
!           count = count + 1
!        endif
!     enddo
!     write(ioOutput,*) "# N-2 AA configurations thrown out: ",count

!   end subroutine get_csfs_off_disk

!   !**************************************************************
  subroutine tree_search_internal 

    !// THIS SUBROUTINE LOOKS FOR PATHS IN THE SPATIAL GRAPH AND COMPUTES 
    !// THEIR LEXICAL INDICES.  THIS IS DONE IN PREPARATION FOR THE 
    !// INDEX VECTOR SUBROUTINE. THIS ROUTINE ONLY LOOKS FOR INTERNAL
    !// SPACE CSFS.  IT IS FOR USE WHEN OPERATING IN VECTORIZED MODE.
    !// CONFIGS WHICH NEED TO BE ELIMINATED SHOULD BE THROWN OUT BY SETTING
    !// THEIR ELEMENTS OF THE INDEX VECTOR TO SOME NEGATIVE NUMBER.
    !// ONE OTHER BYPRODUCT OF THIS ROUTINE ARE THE VECTORS WHICH KEEP
    !// TRACK OF THE TOTAL NUMBER OF SINGLES IN EACH INTERNAL 
    !// PATH
    !//
    !// CALLED BY:
    !//    GRAPH_DRIVER

    use global_var_mod
    use utilities_mod
    use ci_utilities_mod
    use graph_var_mod
    use molecule_var_mod

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::max_cfgs                            !// MAXIMUM NUMBER OF CONFIGURATIONS
    integer::allocatestatus                      !// DYNAMIC MEMORY STATUS VARIABLE
    integer::vert                                !// KEEPS TRACK OF THE CURRENT VERTEX
    integer::orbital                             !// KEEPS TRACK OF CURRENT ORBITAL NUMBER
    integer::electrons                           !// KEEPS TRACK OF NUMBER OF ELECTRONS
    integer::paths                               !// NUMBER OF COMPLETE PATHS
    integer::arc                                 !// ARC OR VERTEX WEIGHT

    integer, dimension(:), allocatable::singles  !// KEEPS TRACK OF NUMBERS OF SINGLES WE HAVE IN EACH PATH
    integer, dimension(:), allocatable::change   !// KEEPS TRACK OF STEPS WE'VE TAKEN DOWN THE DRT
    integer, dimension(:), allocatable::lex      !// KEEPS TRACK OF LEXICAL INDEX AT EACH LEVEL


!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    max_cfgs = abs(y0(num_vert))

    !// ALLOCATE THE SINGLES VECTOR
    allocate(singles(0:num_internal), stat = allocatestatus)
    call allocatecheck(allocatestatus,"singles ")

    !// ALLOCATE THE CHANGE VECTOR
    allocate(change(0:num_orbitals), stat = allocatestatus)
    call allocatecheck(allocatestatus,"change  ")

    !// ALLOCATE THE LEXICAL INDEX VECTOR
    allocate(lex(0:num_internal), stat = allocatestatus)
    call allocatecheck(allocatestatus,"lex     ")

    !// LETS ALLOCATE SOME SPACE FOR THE INDEX VECTORS AND NUMBERS OF 
    !// SINGLES
    if (.not.allocated(internal_index_vector0)) then

       allocate(internal_index_vector0(abs(y0(vertex(num_internal,num_elec)))),stat=allocatestatus)
       allocate(internal_index_vector1(abs(y0(vertex(num_internal,num_elec-1)))),stat=allocatestatus)
       allocate(internal_index_vector2(abs(y0(vertex(num_internal,num_elec-2)))),stat=allocatestatus)
       allocate(internal_index_vector3(abs(y0(vertex(num_internal,num_elec-2)))),stat=allocatestatus)
       call allocatecheck(allocatestatus,"interna2")

       allocate(v_singles(abs(y0(vertex(num_internal,num_elec)))),stat=allocatestatus)
       allocate(nm1_singles(abs(y0(vertex(num_internal,num_elec-1)))),stat=allocatestatus)
       allocate(nm2_singles(abs(y0(vertex(num_internal,num_elec-2)))),stat=allocatestatus)
       call allocatecheck(allocatestatus,"singles ")

    endif

    max_cfgs = max(abs(y0(vertex(num_internal,num_elec-2))),&
         abs(y0(vertex(num_internal,num_elec-1))),&
         abs(y0(vertex(num_internal,num_elec-0))))


!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// NOW WE ARE GOING TO GO INTO OUR TREE SEARCH ALGORITHM TO FIND ALL THE PATHS
    !// FROM THE TOP TO THE BOTTOM AND THEIR LEXICAL INDICES.  FIRST, INITIALIZE THINGS
    !// AND THEN HAVE A LOOP WHICH KEEPS LOOPING TILL WE FINISH
    vert = 1
    change = 0
    change(0) = 2
    orbital = 0
    singles = 0
    electrons = 0
    paths = 0
    lex = 0
    lex(0) = 1

    internal_index_vector0 = 0
    internal_index_vector1 = 0
    internal_index_vector2 = 0
    internal_index_vector3 = 0

    repeater: do 

       !// STARTING FROM THE CURRENT VERTEX, SEE IF WE CAN STEP DOWN.  WE
       !// TRY STEP 2 FIRST, THEN STEP 1, AND FINALLY STEP 0
       if (change(orbital) == 2) then

          arc = y2(vert)

          if (arc >= 0) then  !// THE STEP IS ALLOWED, WORK WITH IT.

             !// INCREMENT NUMBER OF ELECTRONS
             electrons = electrons + 2 

             !// INCREMENT NUMBER OF ORBITALS
             orbital = orbital + 1

             !// GET NEW VERTEX FOR NEXT ITERATION
             vert = auxiliary(orbital) + electrons

             !// ADD IN RELEVANT NUMBER OF SINGLES
             singles(orbital) = singles(orbital - 1)

             !// ADJUST LEXICAL INDEX
             lex(orbital) = lex(orbital - 1) + arc
             if (lex(orbital) > max_cfgs) then
                write(ioError,*) "Maximum number of configurations is:", max_cfgs
                write(ioError,*) "Current lexical index is:", lex(orbital)
                write(ioError,*) "Paths processed:", paths
                write(ioError,*) "Current level:",orbital
                write(ioError,*) "Trying change:",change(orbital)
                write(ioError,*) "Good luck debugging this one"
                stop
             endif

             !// PROCESS THE PATH IF WE ARE AT THE BOTTOM, AND STEP BACK
             if (orbital == num_internal) then 

                if (check_excitation().and.&
                     check_frozen()) then

                   if (electrons == num_elec) then
                      internal_index_vector0(lex(orbital)) = singles(orbital)
                      v_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-1) then
                      internal_index_vector1(lex(orbital)) = singles(orbital)
                      nm1_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-2) then
                      internal_index_vector2(lex(orbital)) = singles(orbital)
                      internal_index_vector3(lex(orbital)) = singles(orbital)
                      nm2_singles(lex(orbital)) = singles(orbital)
                   else
                      write(ioError,*) "I'm not at the internal external space border."
                      stop
                   endif

                   paths = paths + 1
                   orbital = orbital - 1
                   change(orbital) = change(orbital) - 1
                   electrons = electrons - 2
                   vert = auxiliary(orbital) + electrons
                   cycle repeater

                else

                   if (electrons == num_elec) then
                      internal_index_vector0(lex(orbital)) =  kill_number
                      v_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-1) then
                      internal_index_vector1(lex(orbital)) =  kill_number
                      nm1_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-2) then
                      internal_index_vector2(lex(orbital)) = kill_number
                      internal_index_vector3(lex(orbital)) = kill_number
                      nm2_singles(lex(orbital)) = singles(orbital)
                   else
                      write(ioError,*) "I'm not at the internal external space border."
                      stop
                   endif

                   paths = paths + 1
                   orbital = orbital - 1
                   change(orbital) = change(orbital) - 1
                   electrons = electrons - 2
                   vert = auxiliary(orbital) + electrons
                   cycle repeater

                endif
             endif

             !// IF WE ARE NOT AT THE END, JUST CYCLE REPEATER STARTING AT
             !// THE NEXT LEVEL WITH CHANGE = 2
             change(orbital) = 2
             cycle repeater

          else  !// THE STEP IS NOT ALLOWED, TRY CHANGE =1.  DON'T CHANGE THE VERTEX
             change(orbital) = 1
             cycle repeater
          endif

       elseif (change(orbital) ==1) then

          arc = y1(vert)

          if (arc >= 0) then  !// THE STEP IS ALLOWED, WORK WITH IT.

             !// INCREMENT NUMBER OF ELECTRONS
             electrons = electrons + 1

             !// INCREMENT NUMBER OF ORBITALS
             orbital = orbital + 1

             !// GET NEW VERTEX FOR NEXT ITERATION
             vert = auxiliary(orbital) + electrons

             !// ADD IN RELEVANT NUMBER OF SINGLES
             singles(orbital) = singles(orbital - 1) + 1

             !// ADJUST LEXICAL INDEX
             lex(orbital) = lex(orbital - 1) + arc

             if (lex(orbital) > max_cfgs) then
                write(ioError,*) "Maximum number of configurations is:", max_cfgs
                write(ioError,*) "Current lexical index is:", lex(orbital)
                write(ioError,*) "Paths processed:", paths
                write(ioError,*) "Current level:",orbital
                write(ioError,*) "Trying change:",change(orbital)
                write(ioError,*) "Good luck debugging this one"
                stop
             endif

             !// PROCESS THE PATH IF WE ARE AT THE BOTTOM, AND STEP BACK
             if(orbital == num_internal) then 

                if (check_excitation().and.&
                     check_frozen()) then

                   if (electrons == num_elec) then
                      internal_index_vector0(lex(orbital)) = singles(orbital)
                      v_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-1) then
                      internal_index_vector1(lex(orbital)) = singles(orbital)
                      nm1_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-2) then
                      internal_index_vector2(lex(orbital)) = singles(orbital)
                      internal_index_vector3(lex(orbital)) = singles(orbital)
                      nm2_singles(lex(orbital)) = singles(orbital)
                   else
                      write(ioError,*) "I'm not at the internal external space border."
                      stop
                   endif

                   paths = paths + 1
                   orbital = orbital - 1
                   change(orbital) = change(orbital) - 1
                   electrons = electrons - 1
                   vert = auxiliary(orbital) + electrons
                   cycle repeater

                else

                   if (electrons == num_elec) then
                      internal_index_vector0(lex(orbital)) = kill_number 
                      v_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-1) then
                      internal_index_vector1(lex(orbital)) = kill_number
                      nm1_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-2) then
                      internal_index_vector2(lex(orbital)) = kill_number
                      internal_index_vector3(lex(orbital)) = kill_number
                      nm2_singles(lex(orbital)) = singles(orbital)
                   else
                      write(ioError,*) "I'm not at the internal external space border."
                      stop
                   endif

                   paths = paths + 1
                   orbital = orbital - 1
                   change(orbital) = change(orbital) - 1 
                   electrons = electrons - 1
                   vert = auxiliary(orbital) + electrons
                   cycle repeater

                endif

             endif

             !// IF WE ARE NOT AT THE END, JUST CYCLE REPEATER STARTING AT
             !// THE NEXT LEVEL WITH CHANGE = 0
             change(orbital) = 2 

             cycle repeater

          else  !// THE STEP IS NOT ALLOWED, TRY CHANGE = 0
             change(orbital) = 0 
             cycle repeater
          endif

       elseif (change(orbital) == 0) then

          arc = y0(vert)
          if (arc >= 0) then  !// THE STEP IS ALLOWED, WORK WITH IT.

             !// INCREMENT NUMBER OF ORBITALS
             orbital = orbital + 1

             !// GET NEW VERTEX FOR NEXT ITERATION
             vert = auxiliary(orbital) + electrons

             !// ADD IN RELEVANT NUMBER OF SINGLES
             singles(orbital) = singles(orbital - 1) 

             !// ADJUST THE LEXICAL INDEX
             lex(orbital) = lex(orbital -1)

             !// PROCESS THE PATH IF WE ARE AT THE BOTTOM, AND STEP BACK
             if(orbital == num_internal) then 

                if (check_excitation().and.&
                     check_frozen()) then
                   if (electrons == num_elec) then
                      internal_index_vector0(lex(orbital)) = singles(orbital)
                      v_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-1) then
                      internal_index_vector1(lex(orbital)) = singles(orbital)
                      nm1_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-2) then
                      internal_index_vector2(lex(orbital)) = singles(orbital)
                      internal_index_vector3(lex(orbital)) = singles(orbital)
                      nm2_singles(lex(orbital)) = singles(orbital)
                   else
                      write(ioError,*) "I'm not at the internal external space border."
                      stop
                   endif

                   paths = paths + 1
                   orbital = orbital - 1
                   change(orbital) = change(orbital) - 1
                   vert = auxiliary(orbital) + electrons
                   cycle repeater
                else
                   if (electrons == num_elec) then
                      internal_index_vector0(lex(orbital)) = kill_number
                      v_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-1) then
                      internal_index_vector1(lex(orbital)) = kill_number
                      nm1_singles(lex(orbital)) = singles(orbital)
                   elseif (electrons == num_elec-2) then
                      internal_index_vector2(lex(orbital)) = kill_number
                      internal_index_vector3(lex(orbital)) = kill_number
                      nm2_singles(lex(orbital)) = singles(orbital)
                   else
                      write(ioError,*) "I'm not at the internal external space border."
                      stop
                   endif

                   paths = paths + 1
                   orbital = orbital - 1
                   change(orbital) = change(orbital) - 1
                   vert = auxiliary(orbital) + electrons
                   cycle repeater        
                endif
             endif

             !// IF WE ARE NOT AT THE END, JUST CYCLE REPEATER STARTING AT
             !// THE NEXT LEVEL WITH CHANGE = 2 
             change(orbital) =  2        
             cycle repeater

          else  !// THE STEP IS NOT ALLOWED, BACK OFF A LEVEL AND TRY THE NEXT STEP FROM THAT LEVEL
             change(orbital) = change(orbital) -1
             cycle repeater
          endif

       elseif (change(orbital) == -1) then   !// WE RAN OUT OF STEPS FOR THIS LEVEL.  IF WE ARE BACK AT THE TOP,
          !// GET OUT OF HERE.  OTHERWISE, GO BACK ONE MORE LEVEL
          if (orbital /= 0) then
             orbital = orbital - 1
             electrons = electrons - change(orbital)
             change(orbital) = change(orbital) - 1
             vert = auxiliary(orbital) + electrons
          else
             exit repeater  !WE FINISHED, YIPEE!!
          endif

       endif

    enddo repeater

    !// LETS WRITE SOME STATISTICS OUT TO THE OUTPUTFILE
    write(ioOutput,*)"Tree search algorithm completed....."
    write(ioOutput,*)"Paths processed: ", paths
    write(ioOutput,*)

    deallocate(lex,singles,change)


  contains

    function check_excitation()

      !// THIS FUNCTION RETURNS TRUE IF THE CONFIGURATION IS NOT MORE
      !// THAN DOUBLY EXCITED FROM THE REFERENCE(S)

      logical::check_excitation
      integer::i,j,difference, electrons

      check_excitation = .false.
      electrons = sum(change(0:num_internal-1))

      !// LOOP OVER REFERENCES AND TEST THEM ALL
      do i = 1, size(references,2)

         !// LOOP OVER ORBITALS IN REFERENCES
         difference=0
         do j = 0, num_internal-1
            difference = difference + abs(change(j) - references(j+1,i))
         enddo

         if (electrons == num_elec) then
            difference = difference
         elseif (electrons == num_elec-1) then
            difference = difference + 1
         elseif (electrons == num_elec-2) then
            difference = difference + 2
         else
            write(ioError,*) "I am not at the internal/external space border."
            write(ioError,*) "Number of electrons: ", electrons
            stop
         endif

         if (difference <= 4) then
            check_excitation =.true.
            return
         endif

      enddo

    end function check_excitation

    function check_frozen()

      !// THIS FUNCTION RETURNS TRUE IF THE CONFIGURATION DOES
      !// NOT HAVE EXCITATIONS OUT OF FROZEN ORBITALS

      logical::check_frozen
      integer::j,difference

      check_frozen = .false.


      !// LOOP OVER ORBITALS IN REFERENCES
      difference=0
      do j = 0, num_frozen-1
         difference = difference + abs(change(j) - 2)
      enddo

      if (difference == 0) then
         check_frozen =.true.
         return
      endif

    end function check_frozen

  end subroutine tree_search_internal


  !**************************************************************
  subroutine index_vector_internal(loc_scr)

    !// THIS SUBROUTINE GENERATES THE INDEX VECTOR FROM THE SINGLES INFORMATION
    !// AND THE FSN ARRAY GENERATED IN TREE_SEARCH AND GET_FSN, RESPECTIVELY
    !// CALLED BY:
    !//    GRAPH_DRIVER
    !// THIS ROUTINE IS JUST FOR THE VECTORIZED MODE, WHERE WE ONLY NEED THE 
    !// INTERNAL CSF VECTORS. 
    !// 
    !// THE MAP OF THE ENTIRE CI VECTOR IS SCHEMATICALLY LAYED OUT BELOW.
    !//
    !//                 _
    !//                 |
    !//     VALENCE     | INTERNAL CSF 1 => INTERNAL_INDEX_VECTOR_0(1)
    !//                 |    SPIN 1
    !//                 |    SPIN 2
    !//                 |      :
    !//                 | INTERNAL CSF 2 => INTERNAL_INDEX_VECTOR_0(2)
    !//                 |       :
    !//                 |
    !//                 -
    !//                 |
    !//                 | INTERNAL CSF 1 => INTERNAL_INDEX_VECTOR_1(1)
    !//     N-1 STATES  |    SPIN 1
    !//                 |       VIRT 1   
    !//                 |       VIRT 2
    !//                 |         :
    !//                 |    SPIN 2
    !//                 |       VIRT 1
    !//                 |         :
    !//                 | INTERNAL CSF 2 => INTERNAL_INDEX_VECTOR_1(2)
    !//                 |    SPIN 1
    !//                 |       :
    !//                 | INTERNAL CSF 3 => INTERNAL_INDEX_VECTOR_1(3)
    !//                 |       :
    !//                 |
    !//                 -
    !//                 |
    !//     N-2 STATES  |  INTERNAL CSF 1
    !//                 !     STATES HAVING ONE VIRT DOUBLY OCCUPIED => INTERNAL_INDEX_VECTOR_3(1)
    !//                 |       SPIN 1
    !//                 |           VIRT 1 DOUBLY OCCUPIED
    !//                 |           VIRT 2 DOUBLY OCCUPIED
    !//                 |                 :
    !//                 |       SPIN 2
    !//                 |           VIRT 1 DOUBLY OCCUPIED
    !//                 |                 :
    !//                 |     STATES HAVING TWO VIRTS DOUBLE OCCUPIED => INTERNAL_INDEX_VECTOR_2(1)
    !//                 |       SPIN 1
    !//                 |           VIRT 1, VIRT 2
    !//                 |           VIRT 1, VIRT 3
    !//                 |           VIRT 1, VIRT 4
    !//                 |                 :
    !//                 |           VIRT 2, VIRT 3
    !//                 |           VIRT 2, VIRT 4
    !//                 |                 :
    !//                 |       SPIN 2
    !//                 |           :
    !//                 |   INTERNAL CSF 2
    !//                 |         :
    !//                 |
    !//                 -
    !//                  

    use graph_var_mod
    use molecule_var_mod
    use utilities_mod
    use global_var_mod
    use locist_mod
    use cholesky_structs
    use locist_var_mod,only:locist_scratch
    use three_four_seg_var_mod
    use two_seg_var_mod

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::allocatestatus                    !// DYNAMIC MEMORY ALLOCATION VARIABLE
    integer::shell1                            !// FIRST # OF OPEN SHELLS CONSISTENT WITH SPIN 
    integer::shellU                            !// FIRST UNALLOWED # OF OPEN SHELLS 
    integer::i,j,k                           !// LOOP VARIABLES
    integer::removed                           !// NUMBER OF REMOVED CONFIGURATIONS
    integer::kept                              !// THE NUMBER OF PATHS WE KEEP
    integer::singles                           !// NUMBER OF SINGLES         
    integer::num_csfs                          !// NUMBER OF CSFS, TOTAL OR ASSOCIATED WITH A GIVEN CONFIGURATION
    integer::cfgs_w_i                          !// CONFIGURATIONS WITH I OPEN SHELLS
    integer::dim                               !// NUMBER OF SPIN FUNCTIONS ASSOCIATED WITH 
    !// A GIVEN NUMBER OF SINGLES - SHORT FOR DIMENSION
    integer::num_valence                       !// NUMBER OF VALENCE STATES
    integer::num_nm1                           !// NUMBER OF N-1 STATES
    integer::num_nm2                           !// NUMBER OF N-2 STATES
    integer::index_start                       !// FOR SETTING UP INDEX VECTOR IN N-1 AND N-2 CASES
    integer, dimension(:), allocatable::temp   !// TEMPORARY WORK ARRAY
    integer, dimension(:), allocatable::last   !// LAST CSF WITH I OPEN SHELLS, CORE(XINFOF + 1) OF DUCH'S CODE

    integer::num_virtC2

    type(locist_scratch)::loc_scr


!!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// WE ARE GOING TO NEED A TEMPORARY ARRAY.  THE ARRAY IS DIMENSIONED TO BE OF
    !// SIZE EQUAL TO THE MAX NUMBER OF OPEN SHELLS, OR 2 WHICH EVER IS LARGER
    allocate(temp(0:max(fsn_size,2)),stat = allocatestatus)
    call allocatecheck(allocatestatus,"temp    ")
    allocate(last(0:max(open_shells,2)),stat = allocatestatus)
    call allocatecheck(allocatestatus,"last    ")

    temp = 0

    !// NOW WE DEFINE SHELL1 AND SHELLU.  AS AN EXAMPLE, IF WE HAVE SINGLET SPIN
    !// THEN THE FIRST ALLOWED VALUE OF THE NUMBER OF OPEN SHELLS IS EQUAL TO 0.  
    !// PREVIOUSLY WE ADJUSTED OPEN_SHELLS SO THAT IT WAS EVEN FOR ODD SPIN MULTIPLICITIES
    !// AND ODD FOR EVEN SPIN MULTIPLICITIES.  
    shell1 = mod(open_shells,2)
    shellU = shell1 + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! We need to make sure that the index vectors are allocated
    if (.not.allocated(internal_index_vector0)) then
       allocate(internal_index_vector0(abs(y0(vertex(num_internal,num_elec)))),stat=allocatestatus)
       call allocatecheck(allocatestatus,"interna0")
       allocate(internal_index_vector1(abs(y0(vertex(num_internal,num_elec-1)))),stat=allocatestatus)
       call allocatecheck(allocatestatus,"interna1")
       allocate(internal_index_vector2(abs(y0(vertex(num_internal,num_elec-2)))),stat=allocatestatus)
       call allocatecheck(allocatestatus,"interna2")
       allocate(internal_index_vector3(abs(y0(vertex(num_internal,num_elec-2)))),stat=allocatestatus)
       call allocatecheck(allocatestatus,"interna2")
       internal_index_vector0 = 0
       internal_index_vector1 = 0
       internal_index_vector2 = 0
       internal_index_vector3 = 0
    endif

    !// TAKE CARE OF THE LOCAL CI DATA.  IN THIS ROUTINE WE ARE GOING TO MODIFY 
    !// THE INDEX VECTOR.  IF THERE ARE CONFIGURATIONS THAT SHOULD BE THROWN OUT
    !// ON THE BASIS OF THE SPATIAL CHARACTERSTICS OF THE ORBITALS, THE ELEMENT OF 
    !// THE INDEX VECTOR GETS SET TO -1.  BELOW, THIS WILL CAUSE OT TO BE THROWN OUT.  
    !// THE FOLLOWING ROUTINE CAN ALSO BE USED TO TAKE CARE OF THROWING OUT 
    !// CONFIGURATIONS FOR DOING A VALENCE ONLY CI
    call local_driver(loc_scr)      
!!!!!!!!!!!!!!!!!!!!!
    !//
    !// VALENCE STATES
    !//
!!!!!!!!!!!!!!!!!!!!!

    write(ioOutput,*) "Forming index vector for valence states....."

    !// THROW OUT PATHS
    do i = 1, abs(y0(vertex(num_internal,num_elec)))

       if (internal_index_vector0(i) > fsn_size) then
          write(ioOutput,*) "You are requesting a calculation with ", internal_index_vector0(i),"open shells." 
          write(ioOutput,*) "This is alot.  I can do this many if you really want to."
          write(ioOutput,*) "Right now I am set to stop at ", fsn_size, " open shells."
          write(ioOutput,*) "Ha.If you really want to do more, you need to go into global_var_mod"
          write(ioOutput,*) "and change the variable fsn_size."
          stop
       endif

       !// THROW OUT OFFENDING CONFIGURATIONS
       if (internal_index_vector0(i) < (spinM - 1).or.&
            internal_index_vector0(i) > open_shells) then
          internal_index_vector0(i) = shellU
       endif

       temp(internal_index_vector0(i)) = temp(internal_index_vector0(i)) + 1

    enddo

    removed = temp(shellU)

    !// LOOP OVER THE NUMBERS OF OPEN SHELLS CONSISTENT WITH THE SPIN.  
    num_csfs = 0
    do i = spinM-1, open_shells, 2

       dim = fsn(i)
       num_csfs  = num_csfs + temp(i)*dim
       last(i)   = num_csfs

       if (temp(i) > 0) then
          write(ioOutput,270) i,dim,temp(i)
       endif

    enddo

270 format(5x,"For",i3," singles there are ",i5," spinfunctions and ",i10," configurations.")
    total_csfs = num_csfs
    num_valence = num_csfs

    !// NOW, WE DO THE FOLLOWING SO THAT WE WILL HAVE A FLAG FOR THE THROWN OUT
    !// CONFIGURATIONS.  YOU WILL SEE HOW IT ALL WORKS SHORTLY.
    temp(shellU) = -1
    fsn(shellU)  =  0

    !// FOR THE VALENCE STATES THE CSFS ARE ORDERED EXACTLY AS THEY ARE ORDERED 
    !// IN THE CASE THAT WE ARE NOT OPERATING IN VECTORIZED MODE
    num_csfs = 1 
    do i = shell1, open_shells, 2
       cfgs_w_i = temp(i)
       temp(i) = num_csfs
       num_csfs = num_csfs + cfgs_w_i*fsn(i)
    enddo

    !// OFFSETS ARE SET UP, GET THE INDEX VECTOR
    do i = 1, abs(y0(vertex(num_internal,num_elec)))

       singles = internal_index_vector0(i)

       if (singles == shellU) then
          internal_index_vector0(i) = kill_number
       else
          internal_index_vector0(i) = temp(singles)
          temp(singles) = temp(singles) + fsn(singles)
       endif

    enddo

    write(ioOutput,280) num_valence
280 format(5x,"Number of valence CSFs: ",i10) 

!!!!!!!!!!!!!!!!!!!!!
    !//
    !// N-1 STATES
    !//
!!!!!!!!!!!!!!!!!!!!!

    write(ioOutput,*) "Forming index vector for N-1 states....."

    temp = 0
    !// THROW OUT PATHS
    do i = 1, abs(y0(vertex(num_internal,num_elec-1)))

       if (internal_index_vector1(i)+1 > fsn_size) then
          write(ioOutput,*) "You are requesting a calculation with ", internal_index_vector1(i),"open shells." 
          write(ioOutput,*) "This is alot.  I can do this many if you really want to."
          write(ioOutput,*) "Right now I am set to stop at ", fsn_size, " open shells."
          write(ioOutput,*) "Ha.If you really want to do more, you need to go into global_var_mod"
          write(ioOutput,*) "and change the variable fsn_size."
          stop
       endif

       !// THROW OUT OFFENDING CONFIGURATIONS
       if (internal_index_vector1(i)+1 < (spinM - 1).or.&
            internal_index_vector1(i)+1 > open_shells) then
          internal_index_vector1(i) = shellU-1
       endif
          loc_scr%num_virt = num_allowed_virtuals(i,"S")
          temp(internal_index_vector1(i)+1) = &
          temp(internal_index_vector1(i)+1) + loc_scr%num_virt
    enddo

    removed = removed + temp(shellU)

    !// LOOP OVER THE NUMBERS OF OPEN SHELLS CONSISTENT WITH THE SPIN.  
    num_csfs = 0
    do i = spinM-1, open_shells, 2

       dim = fsn(i)
       num_csfs  = num_csfs + temp(i)*dim
       last(i)   = num_csfs

       if (temp(i) > 0) then
          write(ioOutput,270) i,dim,temp(i)
       endif

    enddo

    total_csfs = total_csfs + num_csfs

    !// NOW, WE DO THE FOLLOWING SO THAT WE WILL HAVE A FLAG FOR THE THROWN OUT
    !// CONFIGURATIONS.  YOU WILL SEE HOW IT ALL WORKS SHORTLY.
    temp(shellU) = -1
    fsn(shellU)  =  0

    !// THE INDEX VECTOR HERE IS SET UP HERE DIFFERENTLY THAN IN THE
    !// VECTORIZED MODE.  HERE THE CSFS ARE LISTED IN MANIFOLDS DEFINED BY
    !// THE SPIN FUNCTION.  WITHIN THE MANIFOLD THE CSFS DIFFER ACCORDING
    !// TO THE VIRTUAL WHICH IS OCCUPIED

    !// GET THE INDEX VECTOR
    index_start = num_valence + 1
    do i = 1, abs(y0(vertex(num_internal,num_elec-1)))

       singles = internal_index_vector1(i)+1

       if (singles == shellU) then
          internal_index_vector1(i) = kill_number
       else
          internal_index_vector1(i) = index_start
             loc_scr%num_virt = num_allowed_virtuals(i,"S")
             index_start = index_start + fsn(singles)*loc_scr%num_virt
       endif

    enddo

    num_nm1 = index_start - num_valence - 1 
    write(ioOutput,290)  num_nm1
290 format(5x,"Number of N-1 states: ",i10)

!!!!!!!!!!!!!!!!!!!!!
    !//
    !// N-2 STATES; MIGHT BE KIND OF A BITCH HERE.  
    !// WE NEED TO TAKE INTO ACCOUNT THE CASES WHERE
    !// TWO VIRTUALS ARE OCCUPIED AND ONLY ONE
    !// VIRTUAL IS OCCUPIED
    !//
!!!!!!!!!!!!!!!!!!!!!

    write(ioOutput,*) "Forming index vector for N-2 states....."

    !// WE NEED NM2CSFS FOR SETTING UP THE PSEUDOSPECTRAL IJAB COUPLING
    !// COEFFICIENTS
    nm2csfs = 0
    temp = 0

    !// THROW OUT PATHS
    do i = 1, abs(y0(vertex(num_internal,num_elec-2)))

       if (internal_index_vector2(i)+2 > fsn_size) then
          write(ioOutput,*) "You are requesting a calculation with ", internal_index_vector2(i)+2,"open shells." 
          write(ioOutput,*) "This is alot.  I can do this many if you really want to."
          write(ioOutput,*) "Right now I am set to stop at ", fsn_size, " open shells."
          write(ioOutput,*) "Ha.If you really want to do more, you need to go into global_var_mod"
          write(ioOutput,*) "and change the variable fsn_size."
          stop
       endif

       if (internal_index_vector3(i) > fsn_size) then
          write(ioOutput,*) "You are requesting a calculation with ", internal_index_vector2(i),"open shells." 
          write(ioOutput,*) "This is alot.  I can do this many if you really want to."
          write(ioOutput,*) "Right now I am set to stop at ", fsn_size, " open shells."
          write(ioOutput,*) "Ha.If you really want to do more, you need to go into global_var_mod"
          write(ioOutput,*) "and change the variable fsn_size."
          stop
       endif

       !// THROW OUT OFFENDING CONFIGURATIONS; DO IT FOR BOTH 
       !// POSSIBLE CASES OF TWO ELECTRONS IN THE VIRTUAL SPACE
       if (internal_index_vector2(i)+2 < (spinM - 1).or.&
            internal_index_vector2(i)+2 > open_shells) then
          !            write(6,*) "Removed,spinM - 1,open_shells",internal_index_vector2(i),spinM - 1,open_shells
          internal_index_vector2(i) = shellU-2
       endif

       if (internal_index_vector2(i) /= shellU-2) nm2csfs = nm2csfs + fsn(internal_index_vector2(i)+2)

          loc_scr%num_virt = num_allowed_virtuals(i,"D")
          num_virtC2 = loc_scr%num_virt*(loc_scr%num_virt-1)/2
          temp(internal_index_vector2(i)+2) = temp(internal_index_vector2(i)+2) + num_virtC2

       if (internal_index_vector3(i) < (spinM - 1).or.&
            internal_index_vector3(i) > open_shells) then
          internal_index_vector3(i) = shellU
       endif
          loc_scr%num_virt = num_allowed_virtuals(i,"D")
          temp(internal_index_vector3(i)) = temp(internal_index_vector3(i)) &
               + loc_scr%num_virt
    enddo

    removed = removed + temp(shellU)
    !    write(6,*) "Stage 2 Removed",removed


    !// LOOP OVER THE NUMBERS OF OPEN SHELLS CONSISTENT WITH THE SPIN.  
    num_csfs = 0
    do i = spinM-1, open_shells, 2

       dim = fsn(i)
       num_csfs  = num_csfs + temp(i)*dim
       last(i)   = num_csfs

       if (temp(i) > 0) then
          write(ioOutput,270) i,dim,temp(i)
       endif

    enddo



    total_csfs = total_csfs + num_csfs

    !// NOW, WE DO THE FOLLOWING SO THAT WE WILL HAVE A FLAG FOR THE THROWN OUT
    !// CONFIGURATIONS.  YOU WILL SEE HOW IT ALL WORKS SHORTLY.
    temp(shellU) = -1
    fsn(shellU)  =  0


    !// THE INDEX VECTOR HERE IS SET UP HERE DIFFERENTLY THAN IN THE
    !// NON VECTORIZED MODE.  HERE THE CSFS ARE LISTED IN MANIFOLDS DEFINED BY
    !// THE SPIN FUNCTION.  WITHIN THE MANIFOLD THE CSFS DIFFER ACCORDING
    !// TO THE VIRTUAL WHICH IS OCCUPIED
    !// GET THE INDEX VECTOR

    !    write(6,*) "ShellU",shellU
    index_start = num_valence + num_nm1 + 1
    do i = 1, abs(y0(vertex(num_internal,num_elec-2)))
       singles = internal_index_vector3(i)
       if (singles == shellU) then
          internal_index_vector3(i) = kill_number
       else
          internal_index_vector3(i) = index_start
             loc_scr%num_virt = num_allowed_virtuals(i,"D")
             index_start = index_start + fsn(singles)*loc_scr%num_virt
       endif

       singles = internal_index_vector2(i)+2
       !        write(6,*) "singles",singles
       if (singles == shellU) then
          internal_index_vector2(i) = kill_number
       else
          internal_index_vector2(i) = index_start
             loc_scr%num_virt = num_allowed_virtuals(i,"D")
             num_virtC2 = loc_scr%num_virt*(loc_scr%num_virt-1)/2
             index_start = index_start + fsn(singles)*num_virtC2
       endif

    enddo

    num_nm2 = index_start - num_valence - num_nm1 - 1 
    write(ioOutput,300) num_nm2
300 format(5x,"Number of N-2 states: ",i10)

    !// TOTAL NUMBER OF CONFIGURATIONS KEPT
    kept = abs(y0(num_vert)) - removed

    !// O.K., ALL DONE.  LETS WRITE INFO TO THE OUTPUT.
    write(ioOutput,*)"Index vector finished....." 
    write(ioOutput,240)abs(y0(num_vert)),kept,removed,total_csfs,real(total_csfs,real8)/real(kept,real8)
240 format(1x,/,&
         1x,"Total number of configurations.......",i14,/,&
         1x,"Number of configurations kept........",i14,/,&
         1x,"Number of configurations removed.....",i14,/,&
         1x,"Total number of full CSFs ...........",i14,/,&
         1x,"Average CSFs/configuration...........",f8.2,/,&
         1x,/)

    !// PRINT OUT THE INDEX VECTOR IF REQUESTED
    if (indexprint > 0) then

!!!!!!!
       !//
       !// VALENCE VECTOR
       !//
!!!!!!!
       write(ioOutput,*)
       write(ioOutput,*) "Valence index vector printed below: "
       k = 1
       do i = 1,abs(y0(vertex(num_internal,num_elec)))/15
          write(ioOutput,250)(j,j = k, k + 14),(internal_index_vector0(j),j = k,k + 14)
          k = k + 15 
       enddo
       write(ioOutput,260)(j,j = k, abs(y0(vertex(num_internal,num_elec))))
       write(ioOutput,260)(internal_index_vector0(j), j = k,  abs(y0(vertex(num_internal,num_elec))))

!!!!!!!
       !//
       !// N-1 VECTOR
       !//
!!!!!!!
       write(ioOutput,*)
       write(ioOutput,*) "N-1 index vector printed below: "
       k = 1
       do i = 1,abs(y0(vertex(num_internal,num_elec-1)))/15
          write(ioOutput,250)(j,j = k, k + 14),(internal_index_vector1(j),j = k,k + 14)
          k = k + 15 
       enddo
       write(ioOutput,260)(j,j = k, abs(y0(vertex(num_internal,num_elec-1))))
       write(ioOutput,260)(internal_index_vector1(j), j = k,  abs(y0(vertex(num_internal,num_elec-1))))

!!!!!!!
       !//
       !// N-2 VECTOR FOR TWO OCCUPIED VIRTUALS
       !//
!!!!!!!
       write(ioOutput,*)
       write(ioOutput,*) "N-2/two virtuals occupied index vector printed below: "
       k = 1
       do i = 1,abs(y0(vertex(num_internal,num_elec-2)))/15
          write(ioOutput,250)(j,j = k, k + 14),(internal_index_vector2(j),j = k,k + 14)
          k = k + 15 
       enddo
       write(ioOutput,260)(j,j = k, abs(y0(vertex(num_internal,num_elec-2))))
       write(ioOutput,260)(internal_index_vector2(j), j = k,  abs(y0(vertex(num_internal,num_elec-2))))

!!!!!!!
       !//
       !// N-2 VECTOR FOR ONE OCCUPIED VIRTUAL
       !//
!!!!!!!
       write(ioOutput,*)
       write(ioOutput,*) "N-2/one virtual occupied index vector printed below: "
       k = 1
       do i = 1,abs(y0(vertex(num_internal,num_elec-2)))/15
          write(ioOutput,250)(j,j = k, k + 14),(internal_index_vector3(j),j = k,k + 14)
          k = k + 15 
       enddo
       write(ioOutput,260)(j,j = k, abs(y0(vertex(num_internal,num_elec-2))))
       write(ioOutput,260)(internal_index_vector3(j), j = k,  abs(y0(vertex(num_internal,num_elec-2))))


250    format(1x,15i7,/,&
            1x,15i7,/,&
            /) 
260    format(1x,15i7)

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// 
    !//  WE NEED TO DO THIS SO THAT WE KNOW THAT FOR 0 OPEN
    !//  SHELLS WE HAVE 0 SPIN FUNCTIONS.  THESE NEXT FEW LINES
    !//  OF CODE ARE OF CRITICAL IMPORTANCE.  WE DO THIS HERE, AND
    !//  ONLY HERE BECAUSE WE ARE TRULY DONE USING FSN
    !//  FOR ANYTHING ELSE BUT DETERMINING THE DIMENSION
    !//  OF THE SPIN SPACE.  IN THIS ROUTINE, BEFORE
    !//  THIS POINT, WE WERE USING SOME ELEMENTS
    !//  OF FSN FOR TEMPORARY STORAGE....BUT NOT ANYMORE.
    !// 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 0,spinM-2
       fsn(i) = 0
    enddo

    do i = spinM, open_shells-1,2
       fsn(i) = 0
    enddo

    do i = open_shells+1, fsn_size
       fsn(i) = 0
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !// 
    !// CALL A ROUTINE TO ALLOCATE THE CI
    !// CI VECTOR AND SIGMA VECTOR IN 
    !// THE PAO BASIS
    !//
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !    stop

       allocate(loc_scr%sigma_nm2_lambda(num_external,num_external),loc_scr%sigma_nm2_mu(num_external,num_external),stat=allocatestatus)
       call allocatecheck(allocatestatus,"scepsinm2")

       allocate(loc_scr%sigma_nm1_mu(num_external),stat=allocatestatus)
       call allocatecheck(allocatestatus,"scepsinm1")

       allocate(loc_scr%virt_lam_allow(num_external),loc_scr%virt_mu_allow(num_external),stat=allocatestatus)
       call allocatecheck(allocatestatus,"virt_allow")


  end subroutine index_vector_internal


end module graph_mod

