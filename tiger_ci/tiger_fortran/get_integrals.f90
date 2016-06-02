! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!     
! File:   get_integrals.f90
! Author: Francis Ricci
!
! Created on June 20, 2014, 11:04 AM
!

module get_integrals
    use iso_c_binding
    implicit none
    
    interface
        subroutine c_get_cholesky_data(nCho) bind(c)
            use iso_c_binding
            implicit none
            integer(kind=C_INT64_T), intent(out)::nCho
        end subroutine c_get_cholesky_data
    end interface
    
    interface
        !get two electron integrals for shells ijkl, returning values into TInt
        subroutine c_eval_ijkl(iShell, jShell, kShell, lShell, TInt, nTInt) bind(c)
            use iso_c_binding
            implicit none
            integer(kind=C_INT64_T), intent(in)::iShell
            integer(kind=C_INT64_T), intent(in)::jShell
            integer(kind=C_INT64_T), intent(in)::kShell
            integer(kind=C_INT64_T), intent(in)::lShell
            real(kind=C_DOUBLE), intent(out)::TInt(1)
            integer(kind=C_INT64_T), intent(in)::nTInt
        end subroutine c_eval_ijkl
    end interface
    
    interface
        !get one electron integrals, returning values into one_ints_AO
        subroutine c_get_one_ints(one_ints_AO, number_bas, nao) bind(c)
            use iso_c_binding
            implicit none
            real(kind=C_DOUBLE), intent(out)::one_ints_AO(1)
            integer(kind=C_INT64_T), intent(in)::number_bas
            integer(kind=C_INT64_T), intent(in)::nao
        end subroutine c_get_one_ints
    end interface
    
    interface
        !get overlap integrals, returning values into AOoverlap
        subroutine c_get_overlap(nBas, AOoverlap) bind(c)
            use iso_c_binding
            implicit none
            integer(kind=C_INT64_T), intent(in)::nBas
            real(kind=C_DOUBLE), intent(out)::AOoverlap(1,1)
        end subroutine c_get_overlap
    end interface

    interface
        !get molecular orbitals, returning value through molecular_orbitals
        subroutine c_get_coefficients(nBas, nOrb, molecular_orbitals) bind(c)
            use iso_c_binding
            implicit none
            integer(kind=C_INT64_T), intent(in)::nBas
            integer(kind=C_INT64_T), intent(in)::nOrb
            real(kind=C_DOUBLE), intent(out)::molecular_orbitals(1,1)
        end subroutine c_get_coefficients
    end interface
    
contains
    subroutine get_cholesky_data(nCho)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(out)::nCho
        
        call c_get_cholesky_data(nCho)        
    end subroutine get_cholesky_data

    !get two electron integrals for shells ijkl, returning values into TInt
    subroutine eval_ijkl(iShell, jShell, kShell, lShell, TInt, nTInt)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::iShell
        integer(kind=8), intent(in)::jShell
        integer(kind=8), intent(in)::kShell
        integer(kind=8), intent(in)::lShell
        integer(kind=8), intent(in)::nTInt
        real(kind=8), intent(out)::TInt(nTInt)
            
        !Decrement for c indices
        call c_eval_ijkl(iShell-1, jShell-1, kShell-1, lShell-1, TInt, nTInt)
    end subroutine eval_ijkl
    
    !get one electron integrals, returning values into one_ints_AO
    subroutine get_one_ints(one_ints_AO, number_bas, nao)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::number_bas
        integer(kind=8), intent(in)::nao
        real(kind=8), intent(out)::one_ints_AO(nao)
        
        call c_get_one_ints(one_ints_AO, number_bas, nao)
    end subroutine get_one_ints
    
    !get overlap integrals, returning values into AOoverlap
    subroutine get_overlap(nBas, AOoverlap)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::nBas
        real(kind=8), intent(out)::AOoverlap(nBas,nBas)
        
        call c_get_overlap(nBas, AOoverlap)
    end subroutine get_overlap
    
    !get molecular orbitals, returning value through molecular_orbitals
    subroutine get_coefficients(nBas, nOrb, molecular_orbitals)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::nBas, nOrb
        real(kind=8), intent(out)::molecular_orbitals(nBas,nOrb)
        
        call c_get_coefficients(nBas, nOrb, molecular_orbitals)
    end subroutine get_coefficients
    
end module get_integrals
