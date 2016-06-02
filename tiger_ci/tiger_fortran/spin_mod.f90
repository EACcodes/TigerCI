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
  !//  SPIN_MOD - THIS MODULE CONTAINS ALL THE ROUTINES NECESSARY FOR
  !//  GENERATING THE SYMMETRIC GROUP REPRESENTATION MATRICES.  WE DO THE
  !//  FOLLOWING THINGS HERE:
  !//       1) MAKE BRANCHING DIAGRAM
  !//       2) MAKE SFT GRAPH
  !//       3) MAKE MATRICES FOR ALL CYCLES
  !//       4) MAKE MATRICES FOR ALL ELEMENTARY TRANSPOSITIONS
  !//       5) MAKE ALL CONTRACTED MATRICES (I..J)|J
  !// 
  !//  A FEW ADDITIONAL REFERENCES ARE IN ORDER HERE
  !//       1) EFFICIENT METHOD FOR THE COMPUTATION OF REPRESENTATION
  !//          MATRICES OF THE UNITARY GROUP GENERATORS. W. DUCH, 
  !//          INT. J. QUANTUM CHEM. 27, 59 (1985)
  !//       2) DUCH W (1986) GRMS OR GRAPHICAL REPRESENTATION OF
  !//          MODEL SPACES, VOL I. LECT NOTES CHEM 42. SPRINGER, BERLIN
  !//          HEIDELBERG NEW YORK
  !// 
  !//  WRITTEN BY DEREK WALTER, 1999
  !//  WARNING: THIS CODE DOES NOT CONFORM TO Y2K STANDARDS!! 
!**************************************************************
  
module spin_mod

  use io_unit_numbers
  
  contains
  
subroutine spin
  
  !// THIS IS THE MAIN DRIVER ROUTINE FOR MY SPIN PACKAGE.  IT IS BASED
  !// ON THE SPIN PACKAGE IN DUCH'S SGGA.  
  !// CALLED BY:             CALLS TO:
  !//    MR_CI                  BRANCHING_DIAGRAM
  !//                           SFT_GRAPH 
  !//                           MAKE_CYCLES
  !//                           MAKE_TRANSPOSITIONS
  
  use global_var_mod
  
  implicit none
  
  
  write(ioOutput,130) 
  130 format(/,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/&
              ,1x,"!//                              ",/&
              ,1x,"!// ENTERED SPIN PACKAGE         ",/&
              ,1x,"!//                              ",/&
              ,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/)
  
  !// LETS MAKE THE BRANCHING DIAGRAM FIRST
  call branching_diagram
  
  !// LETS MAKE THE SFT GRAPH
  call sft_graph
  
  !// NOW GENERATE THE MATRICES FOR THE CYCLES 
  call make_cycles
 
  !// NOW MAKE THE TRANSPOSITIONS
  call make_transpositions
  
  !// MAKE THE CONTRACTED CYCLES AND PRINT THEM
  if (contractedprint > 0) call make_all_contracted
  
end subroutine spin
!**************************************************************
subroutine branching_diagram
  
  !// THIS SUBROUTINE MAKES THE BRANCHING DIAGRAM.  
  !// CALLED BY:
  !//    SPIN
  
  use global_var_mod
  use graph_var_mod
  use spin_var_mod
  use molecule_var_mod
  use utilities_mod
  
  implicit none
  
  integer::allocatestatus              !// FOR DYNAMIC MEMORY ALLOCATION
  integer::up_steps                    !// IN ANY PATH, THIS IS THE NUMBER OF ARCS THAT MUST POINT UP (FROM RIGHT TO LEFT)
  integer::down_steps                  !// LIKE UP_STEPS, EXCEPT POINT DOWN FROM RIGHT TO LEFT.  
  integer::a                           !// LABELS THE NUMBER OF UP STEPS NEEDED TO REACH CURRENT VERTEX
  integer::b                           !// LABELS THE NUMBER OF DOWN STEPS NEEDED TO REACH CURRENT VERTEX 
  integer::amax                        !// NUMBER OF UP STEPS NEEDED TO REACH HIGHEST VERTEX IN CURRENT COLUMN
  integer::bmax                        !// NUMBER OF DOWN STEPS NEEDED TO REACH HIGHEST VERTEX IN CURRENT COLUMN
  integer::i,j                         !// LOOP CONTROL VARIABLES
  integer::finish                      !// AN INDEX TO TELL IS WHEN TO STOP PRINTING THE BRANCHING DIAGRAM  
  
  !// WHAT WE DO HERE IS THINK OF THE BRANCHING DIAGRAM AS A MATRIX.  TAKE THE HIGHEST POINT
  !// ON THE BRANCHING DIAGRAM (IT HELPS TO LOOK AT THE FIGURE ON PG. 123 OF THE COMPUTER
  !// PHYSICS REPORTS PAPER WHERE THIS IS THE POINT S = 5/2, K = 3) AND TREAT THIS IS THE CORNER OF A 
  !// MATRIX.  THIS POINT WILL HAVE THE LARGEST ROW DIMENSION AND THE SMALLEST COLUMN DIMENSION.  LETS
  !// ALLOCATE THE SPACE FOR THE BRANCHING DIAGRAM.
  up_steps = (open_shells + spinM - 1)/2
  down_steps = open_shells - up_steps 
  allocate(branch(up_steps,down_steps + 1), stat = allocatestatus)
  call allocatecheck(allocatestatus,"branch  ")
  
  !// THE BASIC ALGORITHM IS AS FOLLOWS.  ONE CAN THINK OF THE VERTICES ON THE BRANCHING DIAGRAM AS BEING
  !// ARRANGED INTO COLUMNS WHERE ALL THE VERTICES IN A GIVEN COLUMN ARE LABELED WITH THE SAME NUMBER OF 
  !// ELECTRONS.  WE START AT THE LEFTMOST EDGE OF THE BRANCHING DIAGRAM AND BUILD UP THE DIAGRAM COLUMN BY COLUMN 
  !// STARTING AT THE HIGHEST SPIN VALUE FOR EACH COLUMN AND MOVING DOWN.  THEN, GO ON TO THE NEXT COLUMN, ETC.
  !// NOTE THAT WE DON'T HAVE THE VERTEX WITH SPIN = 0 AND MAXIMUM NUMBER OF OPEN SHELLS (INCIDENTALLY, THIS IS WHY
  !// WE DON'T DIMENSION BRANCH HAE NUMBER OF ROWS = UP_STEPS + 1).
  !// START AT LEFT MOST VERTEX
  branch = 0
  a = up_steps
  b = down_steps + 1
  amax = up_steps
  bmax = down_steps + 1
  
  !// THE FOLLOWING LOOP WILL ASSIGN THE NEW VERTEX WEIGHT AND THEN CHOOSE THE NEW
  !// VERTEX OVER AND OVER UNTIL WE REACH THE RIGHTMOST VERTEX ON THE BRANCHING DIAGRAM.
  !// THEN IT WILL EXIT.
  repeater: do
      
      !// FIRST, LETS ASSIGN THE NEW VERTEX WEIGHT.
      !// FIRST CASE, WE ARE AT LEFTMOST VERTEX
      if ((a == up_steps).and.(b == down_steps + 1)) then
          branch(a,b) = 1
      
      !// SECOND CASE, BOTTOM VERTEX OF A COLUMN
      elseif ((b == down_steps + 1).or.&
          (a + 1 == b)) then
          branch(a,b) = branch(a+1,b)
      
      !// THIRD CASE, TOP VERTEX IN A COLUMN TO THE RIGHT OF THE APEX
      elseif ((b == 1).and.(a < up_steps)) then
          branch(a,b) = branch(a+1,b) + branch(a, b+1)
      
      !// FOURTH CASE, TOP VERTEX IN COLUMN LEFT OF APEX AND THE APEX
      elseif (a == up_steps) then
          branch(a,b) = 1
      !// FIFTH AND FINAL CASE, INTERIOR VERTEX
      else
          branch(a,b) = branch(a+1,b) + branch(a,b+1)
      
      endif
  
      !// O.K., NOW LETS GET OUR NEW POINT OR EXIT IF WE ARE DONE.
      !// EXIT IF WE ARE AT LAST VERTEX.
      if ((a==1).and.(b==1)) then
          exit repeater
      
      !// NEXT COLUMN OVER 
      elseif ((b == down_steps + 1).or.&   !// RUN OUT OF DOWN STEPS, ONLY HAPPENS WHEN SPINM > 1
          (a + 1 == b).or.&    !// AT THE BOTTOM WITH A \/ SEGMENT
          (a == b)) then       !// AT THE BOTTOM WITH A /\ SEGMENT
          if ((bmax > 1).and.(amax == up_steps)) then        !// WE ARE TO THE LEFT OF THE APEX
              amax = amax
              bmax = bmax -1
              a = amax
              b = bmax
              cycle repeater
          elseif ((bmax == 1).and.(amax <= up_steps)) then  !// RIGHT OF APEX
              amax = amax - 1
              bmax = bmax
              a = amax
              b = bmax
              cycle repeater
          endif
      
      !// ONLY REMAINING POSSIBILITY, NEXT VERTEX DOWN
      else 
          a = a - 1
          b = b + 1
          cycle repeater
      endif  
  
  enddo repeater
  
  !// WE FINISHED THE BRANCHING DIAGRAM.  LETS PRINT IT OUT TO THE OUTPUT FILE. 
  !// ACTUALLY, DON'T PRINT IT IF IT IS TOO BIG
  write(ioOutput,*) "Finished branching diagram....."
  if (branchprint > 0) then
      260 format(1x,"Branching diagram for",i3," open shells and spin multplicity",i3,":")
      write(ioOutput,260) open_shells, spinM
      finish = 1 
      do i = 1, down_steps + 1
          if (i == 1) then
              write(ioOutput,250) (branch(j,i),j = up_steps,1,-1),branch(1,1)
          else 
              write(ioOutput,250) (branch(j,i),j = up_steps,finish,-1)
              !// NEXT TIME THROUGH WE PRINT ONE LESS ELEMENT
              finish = finish + 1
          endif
      enddo
      250 format(//1x,20i7)
  endif
  
  !// FLUSH OUTPUT FOR DEBUGGING
  call flush(ioOutput)
  
end subroutine branching_diagram
      
  
!**************************************************************
subroutine sft_graph
  
  !// THIS SUBROUTINE MAKES THE SFT GRAPH. ONCE THIS ROUTINE IS FINISHED, WE
  !// DON'T NEED THE BRANCHING DIAGRAM ANYMORE AND WE WILL DEALLOCATE THE MEMORY.
  !// THE MAIN ALGORITHM IS AS FOLLOWS:
  !//  1) GENERATE A SPIN FUNCTIONS IN REVERSE LEXICAL ORDER
  !//  2) LOOP OVER ALL ELEMENTARY TRANSPOSITIONS, BUILDING THE COLUMN 
  !//     OF THE SFT GRAPH CORRESPONDING TO THE CURRENT SPIN FUNCTION
  !//  3) GET THE NEXT SPIN FUNCTION, GOTO (2). 
  !// CALLED BY:
  !//    SPIN
  
  use graph_var_mod           
  use spin_var_mod          
  use utilities_mod
  use molecule_var_mod
  use global_var_mod
  
  implicit none
  
  integer::i,j,k                    !// LOOP CONTROL VARIABLES
  integer::a,b                      !// NUMBER OF UP, DOWN STEPS REQUIRED TO REACH CURRENT VERTEX
  integer::allocatestatus           !// FOR DYNAMIC MEMORY
  integer::shells                   !// THE CURRENT NUMBER OF OPEN SHELLS
  integer::weight                   !// THE WEIGHT OF THE CURRENT SPIN FUNCTION
  integer::new_weight               !// THE WEIGHT OF DAUGHTER FUNCTIONS 
  integer::up_steps                 !// NUMBER OF UP_STEPS IN 1 COMPLETE PATH (SPIN FUNCTION)
  integer::down_steps               !// NUMBER OF DOWN STEPS IN 1 COMPLETE PATH (SPIN FUNCTION)
  integer::num_print_cycles         !// NUMBER OF TIMES WE NEED TO PRINT A GROUP OF COLUMNS IN A MATRIX
  integer::funct                    !// LABELS THE SPIN FUNCTION WHEN WE PRINT THINGS OUT
  integer::permutation              !// LABELS THE PERMUTATION WHEN WE PRINT THINGS OUT
  logical::previous_step            !// TRUE IF THE PREVIOUS STEP WAS AN UP STEP (\), FALSE OTHERWISE (/)
  logical::current_step             !// TRUE IF CURRENT STEP IS UP, FALSE OTHERWISE
  logical, dimension(:), allocatable::path           !// THE FULL CURRENT SPIN FUNCTION
  integer, dimension(:,:), allocatable::pathAB       !// THE PATH AS A LIST OF (A,B) VALUES.  FIRST NUMBER IS A, SECOND B.
! START AUTOGENERATED INITIALIZATION 
permutation = 0
current_step = .false.
previous_step = .false.
shells = 0
b = 0
a = 0
j = 0
k = 0
i = 0
num_print_cycles = 0
weight = 0
up_steps = 0
! END AUTOGENERATED INITIALIZATION 
  
  !// ALLOCATE THE SFT GRAPH.  TO GET THE INFO FOR ET (I,I+1) LOOK AT ROW I OF SFT.
  allocate(sft(open_shells-1,fsn(open_shells)), stat = allocatestatus)
  call allocatecheck(allocatestatus,"sft     ")
  
  !// ALLOCATE THE SPIN MATRIX.  HERE WE STORE THE INTERMEDIATE SPINS FOR ALL THE
  !// FUNCTIONS AT ALL OF THE VERTICES.  SOME OF THE INFORMATION IS USELESS BECAUSE
  !// WE KNOW THAT THE SPIN AT SHELLS=0 SHOULD BE SPINM.  ALSO, AT SHELLS =OPEN_SHELLS
  !// AND SHELLS=OPEN_SHELLS-1 THE SPIN MULTIPLICITY SHOULD BE 2 AND 1 RESPECTIVELY.  I DO
  !// IT THIS WAY TO KEEP THE INDEXING SIMPLE.  
  allocate(spin_matrix(0:open_shells,fsn(open_shells)),stat=allocatestatus)
  call allocatecheck(allocatestatus,"spin_mat")
  
  !// ALLOCATE PATH. THIS ARRAY IS LABELLING ARCS, SO IT RUNS FROM 1 - OPEN_SHELLS
  allocate(path(open_shells), stat = allocatestatus)
! AUTOGENERATED INITALIZATION
path = .false.
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"path    ")
  
  !// ALLOCATE PATHAB.  WE ARE STORING VERTEX LABELS SO WE NEED TO GO FROM 0 TO OPEN_SHELLS.
  !// AGAIN SOME OF THE DATA IS USELESS, BUT I DO IT FOR THE SAKE OF SIMPLICITY IN INDEXING.  
  allocate(pathAB(0:open_shells,2),stat = allocatestatus)
! AUTOGENERATED INITALIZATION
pathab = 0
! END AUTOGENERATED INITIALIZATION 
  call allocatecheck(allocatestatus,"pathAB  ")
  
  up_steps = (open_shells + spinM - 1)/2
  down_steps = open_shells - up_steps 
  
  !// INITIALIZE THINGS
  path(open_shells) = .true.
  shells = open_shells - 1
  pathAB(open_shells,1) = 0
  pathAB(open_shells,2) = 1
  pathAB(shells,1) =  1
  pathAB(shells,2) =  1
  weight = 1
  
  !// FILL IN THE FIRST AND LAST ROW MANUALLY 
  spin_matrix(0,1:fsn(open_shells):1) = spinM
  spin_matrix(open_shells,1:fsn(open_shells):1) = 1
  
  repeater: do
      
      !// LETS GET A NEW SPIN FUCTION.  ANOTHER REPEATING LOOP
      !// THAT QUITS WHEN WE GET TO THE END
    get_spin_function: do
  
          !// FIRST TEST TO SEE IF WE GOT A NEW FUNCTION
          if ((pathAB(shells,1) == up_steps).and.&
              (pathAB(shells,2) == down_steps + 1).and.&
              (shells == 0)) then
            exit get_spin_function
              
          !// NOW, WE WANT TO GENERATE THE SPIN FUNCTIONS IN REVERSE 
          !// LEXICAL ORDER.  THIS MEANS THAT WE GO DOWN EVERYCHANCE
          !// WE GET.  WE ONLY GO UP IF WE HAVE TO.  WHICH IS WHAT WE
          !// DO NEXT
          elseif ((pathAB(shells,2) == down_steps + 1).or.&
              (pathAB(shells,2) == pathAB(shells,1) + 1)) then
              path(shells) = .true.
              !// NOW UPDATE OUR VALUES FOR THE NEXT CYCLE  
              shells = shells - 1
              pathAB(shells,1) = pathAB(shells + 1,1) + 1
              pathAB(shells,2) = pathAB(shells + 1,2)
            cycle get_spin_function
          else
              !// OTHERWISE WE GO DOWN.
              path(shells) = .false.
              shells = shells - 1
              pathAB(shells,1) = pathAB(shells + 1,1)
              pathAB(shells,2) = pathAB(shells + 1,2) + 1
            cycle get_spin_function
          endif
      
    enddo get_spin_function
  
      !// O.K., NOW WE HAVE OUR SPIN FUNCTION.  LETS LOOP OVER ELEMENTARY TRANSPOSITIONS
      !// AND BUILD UP THE COLUMN OF THE SFT MATRIX FOR THIS SPIN FUNCTION.  THE ELEMENTARY 
      !// TRANSPOSITIONS WILL BE LABELED BY (I,I+1).  SO SFT(I,X) DESCRIBES WHAT
      !// YOU GET WHEN YOU OPERATE ON SPIN FUNCTION "X" WITH TRANSPOSITION (I,I+1).
  
      !// IN EVERYTHING THAT FOLLOWS IT IS HELPFUL TO KEEP IN MIND THAT IN EACH PATH
      !// THE ARC LABELLED BY "I" POINTS TO THE VERTEX LABELLED BY "I" WHERE VERTEX 0
      !// IS AT THE LEFTMOST EDGE OF THE GRAPH.  SO PATH(I) LABELS THE ARC WHICH POINTS
      !// (FROM LEFT TO RIGHT) TO THE VERTEX LABELED BY PATHAB(I,X) X = 1,2.
  
      do i = open_shells - 1,1,-1
      
          previous_step = path(i+1)  
          current_step = path(i)  
  
          !// THE VALUES OF A AND B FOR THE VERTEX THAT ARC I IN THE
          !// CURRENT SPIN FUNCTION POINTS TO.
          a = pathAB(i,1)
          b = pathAB(i,2)
   
          !// LETS FIRST FILL UP OUR SPIN MATRIX
          if (.not.previous_step) then   !// WE STEPPED DOWN
              spin_matrix(i,weight) = spin_matrix(i+1,weight) - 1
          else
              spin_matrix(i,weight) = spin_matrix(i+1,weight) + 1 
          endif
  
          !// NOW, LETS DO THE SFT GRAPH.  WE ARE GOING TO DO THIS A LITTLE
          !// DIFFERENTLY THAN WHATS DONE IN THE ARTICLES, SO PAY ATTENTION.
          !//                   /
          !//  1) IF WE HAVE A / OR A \ SEGMENT THEN WE STORE 0.  SUCH A SEGMENT DOES NOT
          !//                          \
          !//     CHANGE UNDER THE CURRENT ET (ELEMENTARY TRANSPOSITION).
          !//  2) IF WE HAVE A /\ SEGMENT WE STORE THE LEXICAL INDEX OF THE NEW DAUGHTER
          !//     FUNCTION THE CURRENT ET GENERATES WHEN IT OPERATES ON THE CURRENT
          !//     SPIN FUNCTION.  THIS NEW SPIN FUNCTION IS THE SAME AS THE CURRENT SPIN FUNCTION
          !//     EXCEPT THE SEGMENT IN QUESTION BECOMES \/. 
          !//  3) IF WE HAVE A \/ SEGMENT, WE DO THE SAME AS IN THE /\ CASE EXCEPT WE STORE
          !//     THE NEGATIVE OF THE LEXICAL INDEX OF THE NEW DAUGHTER FUNCTION
          !//  4) IF WE HAVE A /\ SEGMENT WHICH IS ON THE "GROUND" LEVEL OF THE BRANCHING 
          !//     DIAGRAM, THEN WE STORE THE LEXICAL INDEX OF THE CURRENT SPIN FUNCTION.  
          
          !// FIRST,  TREAT CASE (4) IN THE DESCRIPTION ABOVE
          if ((previous_step).and.(.not.current_step).and.&
  
              (a==b)) then
              sft(i,weight) = weight 
  
          !// NOW DO CASE (2) 
          elseif ((previous_step).and.(.not.current_step)) then
  
              !// WE NEED TO SUBTRACT OFF THE WEIGHT OF THE \ ARC IN THE CURRENT 
              !// FUNCTION AND ADD IN THE WEIGHT OF THE \ ARC IN THE DAUGHTER
              !// FUNCTION.
              new_weight = weight - branch(a-1,b+1)
              if (b + 2 <= down_steps + 1) new_weight = new_weight + branch(a-1,b+2)
                  sft(i,weight) = new_weight
  
          !// NOW DO CASE (3)
          elseif ((.not.previous_step).and.(current_step)) then
  
              !// WE NEED TO DO THE SAME THING AS IN CASE(2), BUT THEN REMEMBER
              !// TO CHANGE THE SIGN.  
              new_weight = weight + branch(a,b) 
              if ((b <= a).and.(b < down_steps + 1)) new_weight = new_weight - branch(a,b+1)
              new_weight = -new_weight
              sft(i,weight) = new_weight
  
          !// NOW DO CASE (1)
                !          elseif (previous_step == current_step) then
          elseif ((previous_step .and. current_step) .or. (.not. previous_step &
                 .and. .not. current_step)) then
  
              sft(i,weight) = 0
      
          !// A BIT OF ERROR CHECKING
          else
              write(ioError,*) "This case not covered."
              write(ioError,*) "current_step = ", current_step
              write(ioError,*) "previous_step = ", previous_step
            write(ioError,*) "subroutine sft_graph"
              stop
          endif
  
      enddo
  
      !// O.K., NOW WE HAVE TO GO BACK AND SEARCH FOR A BRANCH POINT SO WE CAN GET
      !// OUR NEW FUNCTION.  THIS IS WHERE ALL THE ACTION IS AS FAR AS GETTING THE
      !// WEIGHTS OF THE NEW SPIN FUNCTIONS.  WHAT WE DO IS TRACE OUR CURRENT SPIN
      !// FUNCTION BACK LOOKING FOR A BRANCH POINT.
      
      !// SET THIS TO 1 BECAUSE WE WIDDLED IT DOWN TO ZERO WHEN WE WORKED BACK.  
      shells = 1
  
      look_for_branch_point: do
          if (shells == open_shells) then
          
              exit repeater
          
          elseif (.not.path(shells).and.&
              (pathAB(shells,1) < up_steps)) then  !// WE FOUND A BRANCH POINT
       
              path(shells)=.true.
              weight = weight + branch(pathAB(shells,1),pathAB(shells,2)+1)
              shells = shells - 1
              pathAB(shells,1) = pathAB(shells + 1,1) + 1
              pathAB(shells,2) = pathAB(shells + 1,2)
              cycle repeater
          
          else   !// JUST STEP LEFT TO RIGHT LOOKING FOR ANOTHER BRANCH POINT
       
              shells = shells + 1
              if ((pathAB(shells,2) < down_steps + 1).and.&
              (pathAB(shells,2) <= pathAB(shells,1)).and.&
              (path(shells)))&
!              (path(shells) == .true.))&
              weight = weight - branch(pathAB(shells,1),pathAB(shells,2) + 1)
       
          endif
      enddo look_for_branch_point
  
  enddo repeater
  
  !// NOW THE SFT GRAPH SHOULD BE DONE.  WE SHOULD ALSO HAVE OUR SPIN MATRIX DONE.  LETS WRITE
  !// THIS STUFF OUT TO THE OUTPUT FILE.  WE WILL WRITE THE MATRICES OUT IN COLUMNS 8 COLUMNS AT 
  !// A TIME.
  write(ioOutput,*)
  write(ioOutput,*) "SFT matrix finished....."
  if (sftprint > 0) then
      write(ioOutput,*) "SFT matrix"
      num_print_cycles= fsn(open_shells)/8 
    funct = 1
      permutation = 1
      do i = 1,num_print_cycles
        write(ioOutput,280) (j,j = funct, funct + 7)
          do k = 1, open_shells -1
            write(ioOutput,290)k,k+1,(sft(k,j), j = funct, funct + 7) 
          enddo
        funct = funct + 8
      enddo
  
      !// NOW DO THE REMAINDERS THAT WE DIDN'T GET
    write(ioOutput,280)(j,j = funct, fsn(open_shells))
    write(ioOutput,*)"----------",(" --------",j = funct, fsn(open_shells)) 
      do k = 1,open_shells-1
        write(ioOutput,290)k,k+1,(sft(k,j), j = funct, fsn(open_shells))
      enddo
  
  
      280 format(1x,/,&
             1x,"    ET    ",1x,8(i5,4x),/,&
             1x,"----------",8(1x,"--------"))
      290 format(1x,"(",i3,",",i3,")",2x,8(i5,4x)) 
      write(ioOutput,*)
  endif
  
  !// BEFORE LEAVING DEALLOCATE UNUSED ARRAYS
  deallocate(pathAB, path, branch)
  
end subroutine sft_graph
  
  
!**************************************************************
subroutine make_cycles
  
  !// THIS SUBROUTINE CREATES THE REPRESENTATION MATRICES FOR ALL THE
  !// CYCLES (K..OSHELL) K =1,2,...,OSHELL-1. BASICALLY, WE BUILD UP THE
  !// MATRICES ONE COLUMN AT A TIME USING THE ALGORITHM DESCRIBED BY DUCH
  !// (REFERENCE 1 IN THE MODULE INTRODUCTION).  WE DECOMPOSE EACH CYCLE
  !// INTO A PRODUCT OF ETS AND THEN APPLY THE ETS ONE AT A TIME TO 
  !// A GIVEN SPIN FUNCTION FILLING IN MATRIX ELEMENTS ALONG THE WAY.
  !// CALLED BY:
  !//    SPIN
  
  use global_var_mod
  use utilities_mod
  use molecule_var_mod
  use spin_var_mod
  use graph_var_mod
  
  implicit none
  
  integer::num_cycles                  !// THE TOTAL NUMBER OF CYCLES
  integer::allocatestatus              !// FOR DYNAMIC MEMORY ALLOCATION
  integer::i,j,k,l                     !// LOOP CONTROL VARIABLES
  integer::cycle_label                 !// TAKES L AND K AND GETS AN INTEGER FOR STORING
  integer::spin_mult                   !// SPIN MULTIPLICITY OF THE CURRENT FUNCTION AT THE CURRENT K VALUE
  integer::num_states                  !// KEEPS TRACK OF NUMBER OF DAUGHTER STATES
  integer::sft_value                   !// THIS IS THE VALUE WE GET FROM THE SFT MATRIX
  integer::non_zero                    !// NUMBER OF NON ZERO ELEMENTS IN THE CYCLE MATRICES
  integer::num_multiplies              !// LETS COUNT THE TOTAL NUMBER OF MULTIPLICATIONS NEEDED
  
  integer, dimension(:), allocatable::reorder      !// A VECTOR WHICH ALLOWS ME TO KEEP TRACK OF THE STATES I AM
                                                   !// PROCESSING
  real(real8)::a_kl                    !// 1/SPIN_MULT
  real(real8)::one                     !// 1.0
  real(real8)::zero                    !// 0.0
  real(real8)::small_number            !// 1.0D-15
  real(real8)::net_result              !// WE HAVE TO KEEP A RUNNING PRODUCT OF THE MATRIX ELEMENT VALUE
  real(real8), dimension(:), allocatable::column    !// TEMPORARILY HOLDS A COLUMN OF THE MATRIX
  
  num_cycles = open_shells*(open_shells - 1)/2 
  non_zero = 0
  num_multiplies = 0
  
  one = real(1.0,real8)
  zero = real(0.0,real8)
  small_number = real(1.0e-9,real8)
  
  !// FIRST, WE NEED TO ALLOCATE THE MEMORY TO STORE THE CYCLES.  THE FIRST INDEX WILL
  !// LABEL THE CYCLE AND TH LAT TWO INDICES WILL LABEL THE INDIVIDUAL MATRIX ELEMENTS IN
  !// THE CURRENT CYCLE.  
  allocate (cycles(fsn(open_shells), fsn(open_shells),num_cycles), stat = allocatestatus)
  call allocatecheck(allocatestatus,"cycles  ")
  
  !// ALLOCATE TEMPORARY ARRAY
  allocate(column(fsn(open_shells)), stat=allocatestatus)
  call allocatecheck(allocatestatus,"column  ")
  
  !// ALLOCATE REORDERING VECTOR
  allocate(reorder(fsn(open_shells)), stat=allocatestatus)
  call allocatecheck(allocatestatus,"reorder ")

  cycles = 0
  column = 0
  reorder = 0


  !// NOW WE'RE READY TO LOOP OVER ALL THE COLUMNS IN THE MATRIX
  columns: do i = 1, fsn(open_shells)
      reorder(1) = i    !// THE FIRST STATE WE OPERATE ON IS THE CURRENT COLUMN WE ARE CALCULATING 
      
      !// NOW LOOP OVER ALL THE POSSIBLE CYCLES (L..K), L < K
      first_index: do k = 2, open_shells
      !// A NEW GROUP OF CYCLES.  START WITH 1 STATE.   
      num_states = 1
      !// REINITIALIZE THE TEMPORARY VECTOR HOLDING THE CURRENT COLUMN 
      column = 0
      column(1) = one
      
      !// NOW COMPUTE ALL THE CYCLES (L..K), L = K-1, K-2, K-3 ...., 1 SIMULTANEOUSLY
      second_index: do l = k-1,1,-1
              
          !// THIS NEXT VARIABLE IS IMPORTANT BECAUSE IT TELLS US HOW WE
          !// OUR INDEXING OUR CYCLES
          cycle_label = (k-1)*(k-2)/2 + l
  
          !// NOW, COMPUTE THE CYCLE (L..K).  IN EVERYTHING THAT FOLLOWS, WE TACK ON EXTRA MINUS
          !// SIGNS TO TAKE INTO ACCOUNT THE (-1)^P FACTOR. THE NEXT LOOP IS OVER ALL THE STATES
          !// THAT HAVE BEEN GENERATED BY THE CURRENT CYCLE.  
          states: do j = 1, num_states
  
              !// LETS GET THE SPIN MULTIPLICITY FOR THE STATE WE ARE GOING TO OPERATE ON.  
                  spin_mult = spin_matrix(l+1,reorder(j))   !// THE SPIN OF THE ITH SPIN FUNCTION AT VERTEX L+1 
                  a_kl = one/real(spin_mult,real8) 
  
          net_result = column(j)        
          !// SEE HOW THE CURRENT ET TRANSFORMS THE CURRENT SPIN FUNCTION
          sft_value = sft(l,reorder(j))
  
          if (sft_value == 0) then       !// (K,L)|I> = (-1)|I>
          
              column(j) = -column(j)
          
          elseif (sft_value > 0) then    !// (K,L)|I> = AK|I> - BK|SFT_VALUE>
          
              column(j) = a_kl*net_result
              num_multiplies = num_multiplies + 1
          
              !// IF THIS ET GENERATES A NEW DAUGHTER STATE, TAKE CARE OF IT NOW.  NOTE
              !// THAT THIS IS THE ONLY PLACE WHERE WE HAVE TO APPLY THE IF-THEN
              !// STATEMENT BECAUSE FOR SFT_VALUE < 0 BY DEFINITION IT MEANS
              !// WE GENERATE A DAUGHTER STATE. 
              if (sft_value /= reorder(j)) then
                  num_states = num_states + 1
                  reorder(num_states) = abs(sft_value)
                  column(num_states) = -sqrt(one - a_kl**2)*net_result
                  num_multiplies = num_multiplies + 1
              endif
          
          elseif (sft_value < 0) then    !// (K,L)|I> = -AK|I> - BK|SFT_VALUE>
          
              column(j) = -a_kl*net_result
              num_multiplies = num_multiplies + 1
              num_states = num_states + 1
              reorder(num_states) = abs(sft_value)
              column(num_states) = -sqrt(one - a_kl**2)*net_result
          
          endif
          enddo states
  
          !// NOW WE HAVE COMPUTED THE ITH COLUMN OF THE MATRIX FOR CYCLE (L..K).
          !// LETS STORE IT.
          do  j = 1, num_states
              if (abs(column(j)) < small_number) column(j) = zero
              cycles(reorder(j),i,cycle_label) = column(j)
          enddo
      
      enddo second_index 
      enddo first_index
  
  enddo columns
  
  !// LETS PRINT OUT THE CYCLES IF REQUESTED
  if (cycleprint > 0) write(ioOutput,*)"Cycles printed below:"
  do k = 2, open_shells
      do l = k-1, 1, -1
         cycle_label = (k-1)*(k-2)/2 + l
         if (cycleprint > 0) write(ioOutput,300)l,k
         do i = 1,fsn(open_shells)
              if (cycleprint > 0) write(ioOutput,310) (cycles(i,j,cycle_label), j = 1,fsn(open_shells))
              do j = 1, fsn(open_shells)
                  if(abs(cycles(i,j,cycle_label)) > small_number) non_zero = non_zero + 1
              enddo
          enddo
      enddo
  enddo
  
  write(ioOutput,320)num_cycles, fsn(open_shells), non_zero, &
        real(non_zero,real8)/real(num_cycles*(fsn(open_shells)**2),real8),&
        real(num_multiplies,real8)/real(non_zero,real8) 
  
  300 format(1x,/,&
         1x,"(",i3,"..",i3,")",&
         /,1x,"--------------")
  310 format(1x,20(1x,f6.3))
  320 format(1x,/,&
         1x,"Statistics for cycle matrices:",/,&
         1x,"Total number of cycles......................", i10,/,&
       1x,"Total number of spin functions..............", i10,/,&
         1x,"Total number of non zero matrix elements....", i10,/,&
         1x,"Percentage of elements which are non zero...", f10.3,/,&
         1x,"Number of multiplications/matrix element....", f10.3)
  
  !// DEALLOCATE ARRAYS
  deallocate(column, reorder)
  
end subroutine make_cycles
  
!**************************************************************
subroutine make_transpositions
  
  !// THIS SUBROUTINE USES THE MATRICES FOR THE CYCLIC PREMUTATIONS TO MAKE 
  !// MATRICES FOR ALL THE TRANSPOSITIONS.  A FEW NOTES HERE.  FIRST, WE ONLY NEED
  !// TO STORE THE UPPER OR LOWER HALF OF THE MATRIX.  THIS IS BECAUSE THE MATRICES 
  !// FOR TRANSPOSITIONS ARE SYMMETRIC (THIS IS BECAUSE THE INVERSE OF THE TRANPOSITION
  !// IS THE TRANSPOSITION ITSELF).  THE METHOD WE ARE GOING TO USE TO COMPUTE THE 
  !// TRANSPOSITIONS IS TO DO A MATRIX MULTIPLICATION :
  !//  (I,J)  =(J+1..I)  (J..I)          I > J
  !//       KL         MK      ML
  !// DON'T FORGET I > J, OTHERWISE THIS WON'T WORK.  
  !// CALLED BY
  !//    SPIN
  
  use global_var_mod
  use utilities_mod
  use molecule_var_mod
  use graph_var_mod
  use spin_var_mod
  
  implicit none
  
  integer::allocatestatus            !// FOR DYNAMIC MEMORY
  integer::num_transps               !// NUMBER OF TRANSPOSITIONS
  integer::size_triangle             !// SIZE OF EACH TRIANGULAR MATRIX
  integer::i,j,k,l                   !// LOOP CONTROL VARIABLES
  integer::transp_label              !// LABEL FOR THE TRANSPOSITION 
  integer::cycle_label               !// LABEL FOR THE CYCLE
  integer::num_sf                    !// MAXIMUM UMBER OF SPIN FUNCTIONS
  integer::start                     !// STARTING POINT OF ROW COPY
  integer::finish                    !// ENDING POINT OF ROW COPY
  
  
  real(real8)::zero                  !// 0.0
  real(real8), dimension(:,:), allocatable::temp_matrix    !// MATRIX TO TEMPORARILY HOLD A MATRIX PRODUCT
  real(real8), dimension(:,:), allocatable::cycle1t        !// MATRIX TO TEMPORARILY HOLD A TRANSPOSED CYCLE MATRIX
  real(real8), dimension(:,:), allocatable::cycle2         !// MATRIX TO TEMPORARILY HOLD A CYCLE MATRIX 
  
  num_transps = (open_shells)*(open_shells -1)/2
  num_sf = fsn(open_shells)
  size_triangle = num_sf*(num_sf+1)/2
  zero = real(0.0, real8)
  
  !// FIRST LETS ALLOCATE THE MEMORY FOR THE TRIANGULAR MATRICES FOR TRANSPOSITION STORAGE.
  !// IN ORDER TO SAVE SPACE ON JUST STORING THE TRIANGULAR MATRICES WE NEED TO USE A 
  !// CANONICAL INDEX TO LABEL THE INDIVIDUAL ELEMENTS.  THIS IS KIND OF A PAIN IN THE
  !// ASS, BUT WE HAVE TO DO IT.
  allocate(transpositions(num_transps, size_triangle), stat = allocatestatus)
  call allocatecheck(allocatestatus,"transpos")
  
  allocate(temp_matrix(num_sf,num_sf),cycle1t(num_sf, num_sf),cycle2(num_sf,num_sf),&
       stat=allocatestatus)
  call allocatecheck(allocatestatus,"temp_mat")
  
  transpositions = zero
  temp_matrix = zero
  cycle1t = zero
  cycle2 = zero
  
  !// WE ARE GOING TO LABEL ALL THE TRANSPOSITIONS USING THE FOLLOWING SCHEME:
  !// FOR TRANSPOSITION (I,J) I > J USE TRANS_LABEL = (I-1)*(I-2)/2 + J 
  
  !// NOW, WE'LL DO ALL THE (I,I+1) TRANSPOSITIONS
  do i = 1, open_shells - 1
      cycle_label = i*(i-1)/2 + i 
      transp_label = cycle_label
  
      !// NOW COPY THINGS OVER
      do j = 1, num_sf 
          start = j*(j-1)/2 + 1
          finish = start + j - 1
          transpositions(transp_label,start:finish:1) = cycles(j,1:j:1,cycle_label)
      enddo   
  
  enddo
  
  !// NOW DO ALL THE OTHER TRANSPOSITIONS
  do i = 3, open_shells
      do j = i-2, 1, -1
          cycle_label = (i-1)*(i-2)/2 + j+1
          cycle1t = transpose(cycles(1:num_sf:1,1:num_sf:1,cycle_label))
          cycle_label = cycle_label-1 
          cycle2 = cycles(1:num_sf:1,1:num_sf:1,cycle_label)
          temp_matrix = zero 
          temp_matrix = matmul(cycle1t,cycle2)
  
          !// NOW, COPY TEMP_MATRIX INTO TRANSPOSITIONS JUST KEEPING THE LOWER HALF
          transp_label = (i-1)*(i-2)/2 + j
          do k = 1,num_sf
              start = k*(k-1)/2 + 1
              finish = start + k - 1
              transpositions(transp_label,start:finish:1) = temp_matrix(k,1:k:1)
          enddo
  
      enddo
  enddo
  
  !// THE FOLLOWING CODE HAS BEEN COMMENTED OUT BECAUSE THE ABOVE PIECE OF CODE IS FASTER.
  !// IF YOU EVER DO A SYSTEM WHICH IS REALLY, REALLY LARGE MAYBE YOU MIGHT WANT TO TRY
  !// THE BOTTOM PIECE OF CODE.
  !// LOOP OVER TRANSPOSITIONS 
  !// DO I = 3, OPEN_SHELLS
  !//     DO J = I-2,1,-1
  !//     TRANSP_LABEL = (I-1)*(I-2)/2 + J
  !//     CYCLE_LABEL = TRANSP_LABEL + 1  !FOR THE (J+1..I) CYCLE 
  !// 
      !// LOOP OVER ELEMENTS IN THE MATRIX TO COMPUTE THE K,L ELEMENT
      !// OF THE TRANSPOSITION MATRIX
  !//     DO K = 1, NUM_SF
  !//     DO L = 1,K  
  !//        ELEMENT_LABEL = K*(K-1)/2 + L    
  
          !// COMPUTE THE NEW ELEMENT
  !//         DO M = 1,NUM_SF
  !//         TRANSPOSITIONS(TRANSP_LABEL,ELEMENT_LABEL) = TRANSPOSITIONS(TRANSP_LABEL, ELEMENT_LABEL) +&
  !//                                 CYCLES(CYCLE_LABEL,M,K)*CYCLES(CYCLE_LABEL-1,M,L)
  !//         ENDDO
  
  !//     ENDDO
  !//     ENDDO
  
  !//     ENDDO
  
  !// ENDDO
  deallocate(temp_matrix, cycle1t, cycle2)
  
  !// NOW LETS PRINT IT OUT AND SEE WHAT WE HAVE.  
  write(ioOutput,330)
  330 format (1x,/,"Finished computing transpositions.....")
  
  !// PRINT OUT TRANSPOSITIONS IF REQUESTED
  if (transprint > 0) then
  
      write(ioOutput,360)   
      
      !// LOOP OVER TRANSPOSITIONS
      do i = 1, open_shells 
          do j = 1, i-1     
              transp_label = (i-1)*(i-2)/2 + j 
              write(ioOutput,340)i,j
              
              !// LOOP OVER ROWS OF MATRIX
              do k = 1, num_sf
                  start = k*(k-1)/2 + 1     
                  finish = start + k -1
                  write(ioOutput,350)(transpositions(transp_label,l),l = start,finish)
              enddo
          enddo
      enddo
  
  endif
  
  
  340 format(1x,/,"(",i3,",",i3,")",&
         1x,/,"--------------------")
  350 format(1x,20(1x,f6.3))
  360 format(1x,"Transposition matrices printed below.",/)
  
end subroutine make_transpositions
  
!*****************************************************************
subroutine contracted_cycle_stateless(i, j, k, l, d1, d2, open_shells_left, open_shells_right, contracted, cycles, sft, spin_matrix)
  
  !// IN THIS SUBROUTINE WE GENERATE THE CONTRACTED CYCLES (L..K)(J..I)|D1,D2.
  !// WE CAN GENERATE EITHER A SINGLE OR DOUBLE CONTRACTION.  THE RESULTING CYCLE 
  !// IS HELD IN THE "CONTRACTED" MATRIX WHICH IS IN THE SPIN_VAR_MOD MODULE.   THIS
  !// ARRAY WILL BE ALLOCATED HERE ONE TIME.  AFTER THIS, IT WILL JUST BE WRITTEN OVER.  
  !// FOR THOSE OF YOU HAVING A TOUGH TIME UNDERSTANDING THIS ROUTINE, IT WILL BE HELPFUL
  !// FOR YOU TO TAKE A LOOK AT PG 151 OF WLODEK DUCH'S GRMS BOOK.  
  
  !// ADDED NOTES:
  !//  TO GENERATE (J..I)|J CALL:
  !//   CONTRACTED_CYCLE(I,J,0,0,J,0,OPEN_SHELLS_LAMBDA,OPEN_SHELLS_MU)
  !// 
  !//  TO GENERATE (L..K)(J..I)|L,J CALL:
  !//   CONTRACTED_CYCLE(I,J,K,L,L,J,OPEN_SHELLS_LAMBDA,OPEN_SHELLS_MU)
  !// 
  !//  THE RESULTS WILL BE IN THE MATRIX "CONTRACTED"
  !// FURTHER NOTES: THIS ROUTINE ASSUMES THAT ALL MATRICES (IN AND OUTPUT) HAVE BEEN ALLOCATED
 
  use global_var_mod,only:real8
  use utilities_mod,only:allocatecheck,deallocatecheck,index2m1
  use graph_var_mod,only:open_shells,fsn
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(real8),dimension(:,:),intent(inout)::contracted
  real(real8),dimension(:,:,:),intent(in)::cycles
  integer,dimension(:,:),intent(in)::sft
  integer,dimension(:,:),allocatable,intent(in)::spin_matrix !KEEP THE ALLOCATABLE AS THE IFORT V2013 OTHERWISE CAUSES WRONG (!!!) RESULTS
  
  integer,intent(in)::d1,d2              !// SEE THE INTRODUCTION TO THE SUBROUTINE
  integer::d1_a,d2_b
  integer,intent(in)::i,j,k,l            !// LABELS THE CYCLES - NOT LOOP CONTROL VARIABLES
  integer::temp,temp1,temp2   !// TEMPORARY VARIABLE
  integer::allocatestatus     !// FOR DYNAMIC MEMORY ALLOCATION
  integer::deallocatestatus   !// ALSO FOR DYNAMIC MEMORY ALLOCATION
  integer,intent(in)::open_shells_right  !// NUMBER OF OPEN SHELLS IN RIGHT CONFIGURATION
  integer,intent(in)::open_shells_left   !// NUMBER OF OPEN SHELLS IN LEFT CONFIGURATION
  integer::max_shells         !// LARGER OF LEFT OR RIGHT OPEN SHELLS
  integer::max_spin_functions !// MAXIMUM NUMBER OF SPIN FUNCTIONS
  integer::nsfl               !// NUMBER OF SPIN FUNCTIONS FOR BRA
  integer::nsfr               !// NUMBER OF SPIN FUNCTIONS FOR KET
  integer::nsfri              !// NUMBER OF SPIN FUNCTIONS NEEDED IN INITIAL CYCLE MULTIPLICATION
  integer::cycle_label        !// LABEL FOR CYCLE 
  integer::spin_function      !// LOOP CONTROL VARIABLE
  integer::column_counter     !// COUNTS THE NONZERO COLUMNS IN THE SFT MATRIX
  integer::spin_mult          !// SPIN_MULTIPLICITY REQUIRED TO COMPUTE A_K
  integer::sft_value          !// VALUE FROM SFT MATRIX
  
  integer::sft_value1         !// VALUES FROM SFT MATRIX INDEXED BY 1 AND 2
  integer::sft_value2         !// FOR USE IN THE DOUBLE CONTRACTION.  
  
  integer::spin_mult1         !// FOR USE IN THE DOUBLE CONTRACTION
  integer::spin_mult2
  integer::o,p,q              !// DUMMY INTEGERS FOR MATRIX MULTIPLICATION
  
  real(real8),parameter::zero = real(0.0,real8),one = real(1.0,real8),two = real(2.0,real8)
  real(real8),parameter::sqrt2 = real(sqrt(real(2.0,real8)),real8)
  real(real8)::a_k            !// NEEDED TO CONSTRUCT COEFFICIENT FOR COMBINING COLUMNS WHEN
                              !// CONSTRUCTING CONTRACTED MATRICES
  real(real8)::a_k1, a_k2     !// USED IN THE DOUBLE CONTRACTIONS
  real(real8)::m              !// +1 OR -1 DEPENDING ON SEGMENT TYPE AT D1+1    
  real(real8)::m1,m2          !// USED ALSO FOR DOUBLE CONTRACTIONS
  real(real8)::coeff_spinfunc !// COEFFICIENT FOR THE COLUMN LABELLED BY SPIN_FUNCTION
  real(real8)::coeff_sftval   !// COEFFICIENT FOR THE COLUMN LABELLED BY SFT_VALUE
  
  real(real8)::coeff_sftval1  !// FOR USE IN THE DOUBLE CONTRACTION
  real(real8)::coeff_sftval2
  real(real8)::coeff_sftval12
  
  real(real8)::scale
  
  real(real8),dimension(:,:), allocatable::cycle_space  !// SPACE TO HOLD A CYCLE
  integer::maxD
  
  !write(*,*) "DEBUG: i j k l d1 d2 in real ",i,j,k,l,d1,d2
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  d1_a = d1
  d2_b = d2
  
  !// DO SOME QUICK ERROR CHECKING
  if (open_shells_left < open_shells_right) then
  
      write(ioOutput,*) "Error in contracted_cycle! "
      write(ioOutput,*) "Open shells in lambda function: ", open_shells_left
      write(ioOutput,*) "Open shells in mu function: ", open_shells_right
      write(ioOutput,*) "Open shells in lambda should be greater than "
      write(ioOutput,*) "or equal to open shells in mu!!"
      stop
  
  endif
  max_spin_functions = fsn(open_shells)
  
  !// FIRST SET THINGS UP SO THAT D1 > D2
  temp1 = min(d1_a,d2_b)
  temp2 = max(d1_a,d2_b)
  d1_a = temp2
  d2_b = temp1
  
  !// IF NECESSARY, ALLOCATE THE ARRAY TO HOLD THE CONTRACTED CYCLE
  !if (.not.allocated(contracted)) then
  !  allocate(contracted(max_spin_functions, max_spin_functions), stat=allocatestatus)
  !  call allocatecheck(allocatestatus,"contract")
  !endif
    
  !// ALLOCATE THE TEMPORARY MATRIX, WHICH WE SHOULD NEED ALWAYS
  if (.not.allocated(cycle_space)) then
    allocate(cycle_space(max_spin_functions,max_spin_functions), stat=allocatestatus)
    call allocatecheck(allocatestatus,"cycle_sp")
  endif
  
  !// GET THE NUMBERS OF SPIN FUNCTIONS
  nsfl = fsn(open_shells_left)
  nsfr = fsn(open_shells_right)
  maxD = max(nsfl,nsfr)
  temp = open_shells_right
  
  !// MAX_SHELLS HERE IS THE LARGEST NUMBER OF OPEN SHELLS
  !// IN EITHER OF THE FUNCTIONS, WHICH I'VE RESTRICED
  !// TO BE OPEN_SHELLS_LEFT.   
  max_shells = open_shells_left
  
  !// DECIDE IF WE NEED TO GET THE ENTIRE CYCLE MATRIX.  BASICALLY,
  !// IF THE POSITION OF THE SINGLET COUPLED PAIR WE ARE SHIFTING IS 
  !// ALREADY AT THE END OF THE SPIN FUNCTION, THEN THE MATRIX FOR SHIFT
  !// IS DIAGONAL IN ITS UPPER NSFR X NSFR SECTION.  BELOW THIS,
  !// IN THE (NSFL-NSFR) X NSFR SECTION, IT IS EQUAL TO ZERO.  THUS, WE
  !// ONLY NEED TO WORK WITH A SMALLER SECTION OF THE MATRIX.  HOWEVER,
  !// IF IT IS NOT THEN WE NEED THE ENTIRE CYCLE MATRIX SO WE ADD IN
  !// EXTRA PAIRS OF OPEN SHELLS AS NECESSARY TO MAKE THE DIMENSION RIGHT.  
  if ( (d1_a > 0).and.(d1_a < max_shells -1) ) temp=temp+2
  if ( (d2_b > 0).and.(d2_b < max_shells -3) ) temp=temp+2
  
  
  !// NOW THIS IS THE SMALLER DIMENSION OF THE CYCLE MATRIX WE NEED TO WORK 
  !// WITH.  
  nsfri  = fsn(temp)
  
  !// REINITIALIZE THE MATRIX WHICH WILL HOLD THE CONTRACTED CYCLE
  contracted(1:maxD,1:maxD) = zero
  cycle_space(1:maxD,1:maxD) = zero
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !// FIRST, WE WILL GET THE FIRST CYCLE.  THIS CYCLE WILL HAVE TO BE NSFL BY NSFRI.  
  !// WE WANT TO MAKE THE CYCLE (J..I)|D1, OR SOME OTHER CONTRACTION.  
  !// WE HAVE ALL THE CYCLES (K..L) WHERE K < L
  !// ALREADY BUILT AND STORED.  THUS, IF J < I WE JUST NEED TO GET THE CYCLE
  !// MATRIX DIRECTLY FROM STORAGE.  HOWEVER, IF J > I WE NEED TO TRANSPOSE THE 
  !// MATRIX THAT WE HAVE STORED.  
  
  !// write(ioOutput,*) "start: ", i,j,k,l,d1,d2
  
  if ( (j-i) < 0) then
  
       !// write(ioOutput,*) "straight"
      
      !// HERE WE JUST GET THE MATRIX WE NEED
      cycle_label = index2m1(i,j)
      contracted(1:nsfl,1:nsfri) = cycles(1:nsfl,1:nsfri,cycle_label)
  
  elseif ( (j-i) > 0) then
  
       !// write(ioOutput,*) "transpose"
  
      !// HERE WE NEED TO TRANSPOSE THE MATRIX
      cycle_label = index2m1(i,j)
      contracted(1:nsfl,1:nsfl) =&
        transpose(cycles(1:nsfl,1:nsfl,cycle_label))
      
      !// I HAVE TAKEN OUT THE CODE BELOW AND REPLACED IT SIMPLY WITH A CALL 
      !// TO TRANSPOSE TWO LINES ABOVE.
      
      !// NOW TRANSPOSE THE MATRIX STORING IT IN CONTRACTED
      !// DO K = 1, MAX_SPIN_FUNCTIONS
      !//     DO L = 1, K
      !//         CYCKE_SPACE(K,L) = CONTRACTED(L,K)
      !//         CYCLE_SPACE(L,K) = CONTRACTED(K,L)
      !//     ENDDO
      !// ENDDO
      !// 
      !// CONTRACTED = ZERO
      !// CONTRACTED(1:NSFL:1,1:NSFRI:1) = CYCLE_SPACE(1:NSFL:1,1:NSFRI:1)
      
  else
  
       !// write(ioOutput,*) "identity: ", nsfri
  
      !// HERE WE JUST NEED AN IDENTITY MATRIX FOR THE CYCLE.  FOR THE
      !// TOP HALF OF THE CYCLE.
      !// AND WE DO NOT USE THE FUNKY IDENTITY ROUTINE AS IT SEEMS VERY "UNKOSHER" (QUOTE DAVID)
      do o=1,nsfri
         contracted(o,o) = one
      enddo
      
  endif
  
  !// NOW WE HAVE OUR FIRST CYCLE.  WE MUST EITHER MAKE A CONTRACTION OF THE
  !// FIRST CYCLE NOW OR WE MUST GO ON TO THE NEXT CYCLE.  HERE WE MAKE EITHER
  !// A DOUBLE OR SINGLE CONTRACTION BASED ON THE VALUES OF D1 AND D2.
  if (d1_a /= 0.and.d2_b == 0) then
  
    !// NOW WE MUST MAKE A SINGLE CONTRACTION.  WE WILL LOOP OVER ALL OF THE
    !// NFSRI SPIN FUNCTIONS (THE COLUMNS OF THE MATRIX) AND CHECK TO SEE 
    !// IF THE SEGMENT (D1,D1+1) IS EITHER A /\ OR \/ SEGMENT.  IF IT 
    !// IS NOT ONE OF THESE SEGMENTS, THEN THE CORRESPONDING COLUMN OF THE
    !// MATRIX FOR THE SHIFT OF THE SINGLET COUPLED PAIR IS POPULATED 
    !// ENTIRELY BY ZEROS.  
      
    column_counter = 0
    cycle_space(1:maxD,1:maxD) = zero
      
    do spin_function = 1, nsfri
          
        !// CALL ARRAY_DEBUG(SIZE(SFT,1),D1,"SFT     ")        
        sft_value = sft(d1_a, spin_function)
          
        !// IF THE SPINS ARE TRIPLET COUPLED AT THIS LEVEL, 
        !// DON'T DO ANYTHING.  
        if (sft_value <= 0) cycle
          
        !// NOW, WE HAVE THE CASE WHERE THERE IS A /\ SEGMENT ON
        !// THE GROUND LEVEL, WE HAVE TO MULTIPLY THE COLUMN BY 
        !// SQRT2.
        if (sft_value == spin_function) then
          
              column_counter = column_counter + 1
              cycle_space(1:nsfl,column_counter) =&
                sqrt2*contracted(1:nsfl,spin_function)
              
        !// NOW, IF WE HAVE A CASE WHERE THERE IS A /\ OR A \/
        !// WE NEED TO FORM A LINEAR COMBINATION OF TWO COLUMNS.  
        else 
          
              column_counter = column_counter + 1
              
              !// GET THE COEFFICIENT TO GO IN FRONT OF THE CURRENT COLUMN.
              !// FIRST GET A_K.
              spin_mult = spin_matrix(d1_a+1,spin_function)
              a_k = one/(real(spin_mult,real8))
              
              !// NOW GET M, AND CORRECT SFT_VALUE IF IT IS LESS THAN 0
              m = one
              if (sft_value > 0) m = -m
              if (sft_value < 0) sft_value = -sft_value
              
              !// NOW, COMPUTE THE COEFFICIENTS
              coeff_spinfunc = -m*sqrt(1 - m*a_k) 
              coeff_sftval =  m*sqrt(1 + m*a_k)
  
              !// DO THE MATRIX MULTIPLICATION
              cycle_space(1:nsfl,column_counter) =&
                coeff_spinfunc*contracted(1:nsfl,spin_function) +&
                coeff_sftval*contracted(1:nsfl,sft_value)
              
        endif
          
    enddo
      
    !// NOW, WE NEED TO COPY THINGS FROM TEMP_MATRIX BACK TO THE CONTRACTED
    !// MATRIX
    contracted(1:maxD,1:maxD) = cycle_space(1:maxD,1:maxD)
      
  elseif (d1_a /= 0.and.d2_b /= 0) then
  
    !// HERE WE NEED TO MAKE A DOUBLE CONTRACTION.  ONCE THE DOUBLE CONTRACTION
    !// IS MADE, IT NEEDS TO BE MUTPLIED BY THE (K..L) CYCLE.  THE TREATMENT HERE
    !// OF THE DOUBLE CONTRACTION IS VERY SIMILAR TO THE TREATMENT OF THE SINGLE
    !// CONTRACTION.  
      
    column_counter = 0
    cycle_space(1:maxD,1:maxD) = zero
      
    do spin_function = 1, nsfri
      
        if (column_counter == nsfr) exit
          
        !// CALL ARRAY_DEBUG(SIZE(SFT,1),D1,"SFT     ")        
        sft_value1 = sft(d1_a, spin_function)
        sft_value2 = sft(d2_b, spin_function)
          
        !// IF THE SPINS ARE TRIPLET COUPLED AT THIS LEVEL, 
        !// DON'T DO ANYTHING.  
        if (sft_value1 <= 0.or.sft_value2 <= 0) cycle
          
        if (sft_value1 == spin_function.and.&
            sft_value2 == spin_function) then
              
              !// NOW, WE HAVE THE CASE WHERE THERE IS A /\ SEGMENT ON
              !// THE GROUND LEVEL, WE HAVE TO MULTIPLY THE COLUMN BY 
              !// SQRT2.
          
              column_counter = column_counter + 1
              cycle_space(1:nsfl,column_counter) =&
                two*contracted(1:nsfl,spin_function)
              
              
        else 
          
              !// NOW, IF WE HAVE A CASE WHERE THERE IS A /\ OR A \/
              !// WE NEED TO FORM A LINEAR COMBINATION OF FOUR COLUMNS.  
          
              column_counter = column_counter + 1
              
              !// GET THE COEFFICIENTS.  HERE THERE ARE GOING TO BE FOUR
              !// COEFFICIENTS BECAUSE WE ARE GOING TO MAKE A LINEAR COMBINATION
              !// OF FOUR ROWS.  
              
              spin_mult1 = spin_matrix(d1_a+1,spin_function)
              a_k1 = one/(real(spin_mult1,real8))
              
              spin_mult2 = spin_matrix(d2_b+1,spin_function)
              a_k2 = one/(real(spin_mult2,real8))
              
              !// NOW GET MS, AND CORRECT SFT_VALUE IF IT IS LESS THAN 0
              m1 = one
              if (sft_value1 > 0) m1 = -m1
              if (sft_value1 < 0) sft_value1 = -sft_value1
              
              m2 = one
              if (sft_value2 > 0) m2 = -m2
              if (sft_value2 < 0) sft_value2 = -sft_value2
              
              !// NOW, COMPUTE THE COEFFICIENTS
              coeff_spinfunc  = (-m1*sqrt(1-m1*a_k1))*(-m2*sqrt(1-m2*a_k2)) 
              coeff_sftval2   = (-m1*sqrt(1-m1*a_k1))*( m2*sqrt(1+m2*a_k2))  
              coeff_sftval1   = ( m1*sqrt(1+m1*a_k1))*(-m2*sqrt(1-m2*a_k2))  
              coeff_sftval12  = ( m1*sqrt(1+m1*a_k1))*( m2*sqrt(1+m2*a_k2)) 
              
              !// HERE IS THE BASIC IDEA IN GRAPHICAL FORM WITH SOME VERY SIMPLE
              !// SPIN FUNCTIONS THAT ILLUSTRATE THE BASIC IDEA.  BASICALLY, FROM
              !// THIS WE CAN GET THE MATRIX REPRESENTATIONS OF THE DOUBLE SHIFT
              !// OPERATORS.  ALSO, FOR COEFF_SFTVAL12 WE SEE THAT THE APPROPRIATE
              !// SPIN FUNCTION TO USE IS THE ONE THAT YOU GET BY TAKING FUNCTION 
              !// SFT_VALUE2 AND FINDING OUT WHAT YOU GET WHEN YOU "INVERT" THE
              !// SPIN COUPLING AT D1.  
              !// 
              !//                       
              !//  (D2..N-2)^2(D1..N)^2{/\/\  }                  {/\/\  }
              !//                      {    \ } = COEFF_SPINFUNC*{    \ }
              !//                      {     \}                  {     \}
              !//                    
              !//                                                {/\    }
              !//                               + COEFF_SFTVAL1 *{  \/\ }
              !//                                                {     \}
              !// 
              !//                                                {  /\  }
              !//                               + COEFF_SFTVAL2 *{\/  \ }
              !//                                                {     \}
              !// 
              !//                                                {      }
              !//                               + COEFF_SFTVAL12*{\/\/\ }
              !//                                                {     \}
              !// 
              !// 
              
              !// DO THE MATRIX MULTIPLICATION
              cycle_space(1:nsfl,column_counter) =&
                coeff_spinfunc*contracted(1:nsfl,spin_function) +&
                coeff_sftval1*contracted(1:nsfl,sft_value1) +&
                coeff_sftval2*contracted(1:nsfl,sft_value2) +&
                coeff_sftval12*contracted(1:nsfl,sft(d1_a,sft_value2))
              
          endif
          
      enddo
      
      !// NOW, WE NEED TO COPY THINGS FROM CYCLE_SPACE BACK TO THE CONTRACTED
      !// MATRIX
      contracted(1:maxD,1:maxD) = cycle_space(1:maxD,1:maxD)
      
  endif
  
  !// NOW WE NEED TO MULTIPLY BY THE (L..K) CYCLE.  SHOULD BE STRAIGHTFORWARD.  
  !// NOTE THAT WE DO THIS MULTIPLICATION IN THE CASE OF NO CONTRACTIONS, A SINGLE 
  !// CONTRACTION, OR A DOUBLE CONTRACTION.  WE CONVENIENTLY USE THIS SAME SECTION
  !// OF CODE IN ALL CASES.  
  if ((l-k) < 0 ) then
  
      !// write(ioOutput,*) "straight"
  
      !// WE CAN USE THE CYCLE STRAIGHT (NO TRANSPOSE NEEDED!)
      cycle_label = index2m1(k,l)
      
      do o = 1, nsfr
          do p = 1, nsfl
              cycle_space(p,o) = zero
          enddo
      enddo        
  
      do q = 1, nsfl
          do p = 1,nsfr
              
              scale = contracted(q,p)
              
              do o = 1,nsfl
                  cycle_space(o,p) = cycle_space(o,p) + cycles(o,q,cycle_label)*scale
              enddo
           
          enddo
      enddo
  
      do o = 1, nsfr
          do p = 1, nsfl
              contracted(p,o) = cycle_space(p,o)
          enddo
      enddo                    
                      
      return
  
  elseif ((l-k) > 0 ) then
  
      !// write(ioOutput,*) "tranpose"
      
      !// WE NEED TO USE THE CYCLE TRANSPOSED
      cycle_label = index2m1(k,l)
      
      do o = 1, nsfl
          do p = 1,nsfr
              
              scale = zero
              do q = 1,nsfl
                  scale = scale + cycles(q,o,cycle_label)*contracted(q,p)
              enddo
              cycle_space(o,p) = scale
              
          enddo
      enddo
      
      do o = 1, nsfr
          do p = 1, nsfl
              contracted(p,o) = cycle_space(p,o)
          enddo
      enddo                
      
      
      return
  endif                                
  
!  110 format(1x,100f10.3)           
!  120 format(1x,"(",i3,"..",i3,")(",i3,"..",i3,")|",i3,",",i3,/,&
!             1x,"--------------------------")           
  
  !// DEALLOCATE THE TEMPORARY MATRIX
  deallocate(cycle_space,stat=deallocatestatus)
  call deallocatecheck(deallocatestatus,"cycle_sp")
  
end subroutine contracted_cycle_stateless
  
  
!*****************************************************************
subroutine make_all_contracted()
  
  !// THIS IS JUST A DUMMY DRIVER ROUTINE THAT COMPUTES ALL THE CONTRACTED CYCLES
  !// AND PRINTS THEM FOR REFERENCE
  
  use global_var_mod
  use spin_var_mod
  use graph_var_mod
  use utilities_mod
  
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!  VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::i,j,k,l         !// LOOP VARIABLES
  integer::allocatestatus,max_spin_functions
  real(real8),dimension(:,:),allocatable::contracted
  
  !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  max_spin_functions = fsn(open_shells)
  allocate(contracted(max_spin_functions, max_spin_functions), stat=allocatestatus)
  call allocatecheck(allocatestatus,"contract")
  
  write(ioOutput,*) "Contracted cycles printed below:"
  write(ioOutput,*)
  do i = 2,open_shells
      do j = 1,i-1
          call contracted_cycle_stateless(j,i,0,0,j,0,open_shells,open_shells-2,contracted,cycles,sft,spin_matrix)
          
          !// NOW WE HAVE THE CYCLE, WRITE IT OUT
          write(ioOutput,*)
          write(ioOutput,100)i,j,j
          do k = 1, fsn(open_shells)
              write(ioOutput,110) (contracted(k,l),l=1,fsn(open_shells-2))
          enddo
          
      enddo
  enddo
  
  call flush(ioOutput)
  !// WRITE(IOOUTPUT,*) "DOUBLY CONTRACTED CYCLES"
  !// I=2
  !// J=4
  !// K=0
  !// L=0
  !// D1=2
  !// D2=0
  
  !// CALL CONTRACTED_CYCLE(I,J,K,L,D1,D2,4,2)
  !// WRITE(IOOUTPUT,120) L,K,J,I,D1,D2
  !// DO K = 1, FSN(4)
  !//     WRITE(IOOUTPUT,110) (CONTRACTED(K,L),L=1,FSN(2))
  !// ENDDO
  !// STOP
              
  !// CALL CONTRACTED_CYCLE(1,4,0,0,2,0,6,2)
  !// WRITE(IOOUTPUT,100)4,1,2
  !// DO K = 1, FSN(OPEN_SHELLS)
  !//     WRITE(IOOUTPUT,110) (CONTRACTED(K,L),L=1,FSN(OPEN_SHELLS-2))
  !// ENDDO
  
  100 format(1x,"(",i3,"..",i3,")|",i3,/,&
             1x,"--------------------------")
             
             
!  120 format(1x,"(",i3,"..",i3,")(",i3,"..",i3,")|",i3,",",i3,/,&
!             1x,"--------------------------")           
             
  110 format(1x,100f10.3)           
             
end subroutine
  
end module spin_mod
