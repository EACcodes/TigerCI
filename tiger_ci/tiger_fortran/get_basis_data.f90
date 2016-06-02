! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!     
! File:   get_basis_data.f90
! Author: Francis Ricci
!
! Created on June 20, 2014, 11:04 AM
!

module get_basis_data
    use iso_c_binding
    implicit none
        
    interface
        !get nuclear repulsion energy, returning value into nuclear_energy
        subroutine c_get_nuclear_repulsion(nuclear_energy) bind(c)
            use iso_c_binding
            implicit none
            real(kind=C_DOUBLE), intent(out)::nuclear_energy
        end subroutine c_get_nuclear_repulsion
    end interface
    
    interface
        !get number of atoms, returning value through num_atoms
        subroutine c_get_num_atoms(num_atoms) bind(c)
            use iso_c_binding
            implicit none
            integer(kind=C_INT64_T), intent(out)::num_atoms
        end subroutine c_get_num_atoms
    end interface
    
    interface
        !get number of shells, returning value through nShell
        subroutine c_nr_shells(nShell) bind(c)
            use iso_c_binding
            implicit none
            integer(kind=C_INT64_T), intent(out)::nShell
        end subroutine c_nr_shells
    end interface
    
    interface
        !get nuclear coordinates, returning value through coordinates
        subroutine c_get_coordinates(n_atoms, coordinates) bind(c)
            use iso_c_binding
            implicit none
        
            integer(kind=C_INT64_T), intent(in)::n_atoms
            real(kind=C_DOUBLE), intent(out)::coordinates(1,1)
        end subroutine c_get_coordinates
    end interface
    
    interface
        !generate mapping of each basis function to its atomic center
        subroutine c_get_basis_map(number_bas, basismap) bind(c)
            use iso_c_binding
            implicit none
            
            integer(kind=C_INT64_T), intent(in)::number_bas
            integer(kind=C_INT64_T), intent(out)::basismap(1)
        end subroutine c_get_basis_map
    end interface

    interface
        !get some datarmation about the basis functions
        subroutine c_basis_func_data(nFuncInShell, number_bas, nShell) bind(c)
            use iso_c_binding
            implicit none
            
            integer(kind=C_INT64_T), intent(in)::number_bas
            integer(kind=C_INT64_T), intent(out)::nFuncInShell(1)
            integer(kind=C_INT64_T), intent(out)::nShell
            
        end subroutine c_basis_func_data
    end interface
        
contains
    !get nuclear repulsion energy, returning value into nuclear_energy
    subroutine get_nuclear_repulsion(nuclear_energy)
        use iso_c_binding
        implicit none
        
        real(kind=8), intent(out)::nuclear_energy
        
        call c_get_nuclear_repulsion(nuclear_energy)        
    end subroutine get_nuclear_repulsion   
    
    !get number of atoms, returning value through num_atoms
    subroutine get_num_atoms(num_atoms)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(out)::num_atoms
        
        call c_get_num_atoms(num_atoms)        
    end subroutine get_num_atoms  
    
    !get number of shells, returning value through nshell
    subroutine nr_shells(nShell)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(out)::nShell
        
        call c_nr_shells(nShell)
    end subroutine nr_shells
        
    !get nuclear coordinates, returning value through coordinates
    subroutine get_coordinates(n_atoms, coordinates)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::n_atoms
        real(kind=8), intent(out)::coordinates(n_atoms,3)
        
        call c_get_coordinates(n_atoms, coordinates)
    end subroutine get_coordinates
    
    subroutine get_basis_map(number_bas, basisfunc2atom)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::number_bas
        integer(kind=8), intent(out)::basisfunc2atom(number_bas)
        
        call c_get_basis_map(number_bas, basisfunc2atom)
    end subroutine get_basis_map
    
    !get some datarmation about the basis functions
    subroutine basis_func_data(nFuncInShell, number_bas, nShell)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::number_bas
        integer(kind=8), intent(out)::nFuncInShell(number_bas)
        integer(kind=8), intent(out)::nShell
        
        call c_basis_func_data(nFuncInShell, number_bas, nShell)
    end subroutine basis_func_data
    
end module get_basis_data
