! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module new_tree_search_mod
    use global_var_mod
    use graph_var_mod
    use new_tree_search_structs
    use ci_utilities_mod
    implicit none

contains
    subroutine init_tree_search(G, top_level, top_occ, bottom_level, bottom_occ)
        implicit none

        integer, intent(in) :: top_level, top_occ, bottom_level, bottom_occ
        type(graph_search_state) :: G
        G % top_level = top_level
        G % top_occ = top_occ
        G % bottom_level = bottom_level
        G % bottom_occ = bottom_occ
        G % ending_vertex = vertex(G % bottom_level, G % bottom_occ)
        G % path_elecs = 0

    end subroutine init_tree_search


    logical function get_internal_path(the_path, G)
        implicit none
        ! inputs
        type(orbital_path) :: the_path
        type(graph_search_state) :: G

        ! local variable
        logical :: finished

        call tree_searchX(the_path, G, finished)
        ! if the tree search finished it did not find another path
        if (finished) then
            get_internal_path = .false.
        else
            get_internal_path = .true.
        end if
    end function get_internal_path

    subroutine tree_searchX(the_path, G, finished)
        implicit none

        ! inputs
        type(orbital_path) :: the_path
        type(graph_search_state) :: G
        logical, intent(out) :: finished

        ! local variables
        real(real8) :: zero

        !!!!!!!!!!!!!!!!!!!!!!! VARIABLE INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zero = real(0.0, real8)
        finished = .false.
        if (G % path_elecs == 0) then
            G % level = G % top_level
            G % step_type = 2
            G % path_elecs = G % top_occ
        end if
        G % current_vertex = vertex(G % level, G % path_elecs)


        ! debug output
        !        write(*,*) "Searching the internal occupation tree"
        !        write(*,*) "--------------------------------------"
        !        write(*,*) "Current level   = " , G % level 
        !        write(*,*) "Current step    = " , G % step_type
        !        write(*,*) "Number of elec  = " , G % path_elecs
        !        write(*,*) "Current vertex  = " , G % current_vertex, " -> (" , G % level , ",", G % path_elecs , ")"
        !        

        !!!!!!!!!!!!!!!!!!!!!!!!!!! START OF ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !// NOW, WE HAVE THE TREE SEARCH ALGORITHM.  WE ARE GOING TO GO FROM THE TOP 
        !// VERTEX TO THE BOTTOM VERTEX.  ALL OF THESE VARIABLES SHOULD BE INITIALIZED IN
        !// THE CALLING ROUTINE.    

        ! If the top level and bottom_level are equivalent then we aren't going to be able to find 
        ! a path between them
        if (G % top_level == G % bottom_level) then
            write(*, *) "Top level and bottom level are equivalent in the treesearch"
            write(*, *) "in the new treesearch this raises an error conditon"
            write(*, *) "whoops"
            write(*, *) "FATAL ERROR"
            stop
        end if


        ! If the current vertex is the ending vertex we already found one path from the starting vertex
        ! to the ending vertex. Step back one level and start searching again
        if (G % current_vertex == G % ending_vertex) then
            G % step_type = the_path % occupations(G % level) - 1
            G % path_elecs = G % path_elecs - the_path % occupations(G % level)

            the_path % arc_weights(G % level) = 0
            if (the_path % occupations(G % level) == 1) then

                the_path % singles(the_path % num_singles) = 0
                the_path % num_singles = the_path % num_singles - 1

            endif
            G % level = G % level - 1
        end if


levels: do

        if (G % level >= 0) G % current_vertex = vertex(G % level, G % path_elecs)

        !// FIRST CHECK TO SEE IF WE SHOULD EXIT
        if ((G % level <= G % top_level).and.(G % step_type == -1)) then
            exit levels

            !// NOW CHECK TO SEE IF WE HAVE A PATH TO PROCESS
        elseif (G % current_vertex == G % ending_vertex) then
            return

            !// NOW TRY A TWO STEP
        elseif (G % step_type == 2) then

            !// TRY TO SEE IF WE CAN MAKE A 2 STEP
            if (add2(G % level, G % path_elecs).and.&
                ((G % path_elecs + 2) <= G % bottom_occ).and.&
                (G % level < G % bottom_level)) then

                G % level = G % level + 1

                the_path % occupations(G % level) = G % step_type
                the_path % arc_weights(G % level) = abs(y2(G % current_vertex))
                G % path_elecs = G % path_elecs + 2
                G % step_type = 2
            else
                G % step_type = 1
                cycle levels
            endif

            !// NOW TRY A 1 STEP    
        elseif (G % step_type == 1) then

            !// TRY TO SEE IF WE CAN MAKE A 1 STEP
            if (add1(G % level, G % path_elecs).and.&
                ((G % path_elecs + 1) <= G % bottom_occ).and.&
                (G % level < G % bottom_level).and.&
                (the_path % num_singles + 1 <= open_shells)) then

                G % level = G % level + 1

                the_path % occupations(G % level) = G % step_type
                the_path % arc_weights(G % level) = abs(y1(G % current_vertex))
                G % path_elecs = G % path_elecs + 1
                the_path % num_singles = the_path % num_singles + 1
                the_path % singles(the_path % num_singles) = G % level
                G % step_type = 2
            else
                G % step_type = 0
                cycle levels
            endif

            !// NOW TRY A 0 STEP
        elseif (G % step_type == 0) then

            !// TRY TO SEE IF WE CAN MAKE A 0 STEP
            if (add0(G % level, G % path_elecs).and.&
                (G % level < G % bottom_level)) then

                G % level = G % level + 1
                the_path % occupations(G % level) = G % step_type
                G % step_type = 2
            else
                G % step_type = -1
                cycle levels
            endif

            !// UH OH, NOW WE NEED TO STEP BACK.  REMEMBER, THE ARCS ARE INDEXED
            !// ACCORDING TO THE LEVELS THEY POINT TO, AND "LEVEL" IS LABELLING
            !// THE ARC WE JUST TRIED TO JUMP OFF FROM.  BASICALLY, ALL I'M TRYING TO
            !// SAY HERE IS PAY ATTENTION AND KEEP THE INDEXING STRAIGHT.   
        elseif (G % step_type == -1) then

            G % step_type = the_path % occupations(G % level) - 1
            G % path_elecs = G % path_elecs - the_path % occupations(G % level)

            the_path % arc_weights(G % level) = 0

            if (the_path % occupations(G % level) == 1) then
                the_path % singles(the_path % num_singles) = 0
                the_path % num_singles = the_path % num_singles - 1
            endif

            G % level = G % level - 1
            cycle levels
        endif

    enddo levels

    ! exit flag
    finished = .true.

end subroutine tree_searchX




end module new_tree_search_mod
