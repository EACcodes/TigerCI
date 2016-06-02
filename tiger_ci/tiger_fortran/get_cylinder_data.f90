! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!     
! File:   get_cylinder_data.f90
! Author: Francis Ricci
!
! Created on June 20, 2014, 11:04 AM
!

module get_cylinder_data
    use iso_c_binding
    use molecule_var_mod
    implicit none
        
    interface
        !get nuclear repulsion energy, returning value into nuclear_energy
        subroutine c_num_cylinders(orb_num, num_cyl) bind(c)
            use iso_c_binding
            implicit none
            
            integer(kind=C_INT64_T), intent(in)::orb_num
            integer(kind=C_INT64_T), intent(out)::num_cyl
        end subroutine c_num_cylinders
    end interface
    
    interface
        ! get manual cylinder endpoints
        subroutine c_get_cylinders(iOrb, cyls) bind(c)
            use iso_c_binding
            implicit none

            integer(kind=C_INT64_T), intent(in)::iOrb
            integer(kind=C_INT64_T), intent(out)::cyls(1)
        end subroutine c_get_cylinders
    end interface

contains
    !get nuclear repulsion energy, returning value into nuclear_energy
    subroutine num_cylinders(orb_num, num_cyl)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::orb_num
        integer(kind=8), intent(out)::num_cyl
        
        call c_num_cylinders(orb_num - 1, num_cyl)        
    end subroutine num_cylinders
    
    ! get manual cylinder endpoints
    subroutine get_cylinders(iOrb, cyls)
        use iso_c_binding
        implicit none
        
        integer(kind=8), intent(in)::iOrb
        integer(kind=8), dimension(:), intent(out)::cyls
        
        call c_get_cylinders(iOrb - 1, cyls)
    end subroutine get_cylinders
    
end module get_cylinder_data

