! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!     
! File:   fetch_globals.f90
! Author: Francis Ricci
!
! Created on June 10, 2014, 3:45 PM
!
!*********************************************

! Module for the extraction of global variables from c++ data structures
module fetch_globals
    
    use molecule_var_mod
    use utilities_mod
    use graph_var_mod
    use global_var_mod
    use get_basis_data
    
#ifdef TIGER_USE_OMP
    use omp_lib
#endif
    
    use iso_c_binding
    implicit none
    
    interface
        ! return int variable name to val
        subroutine c_get_global_int(name, val) bind(c)
            use iso_c_binding
            implicit none
            character(kind=c_char), intent(in)::name(*)
            integer(kind=C_INT64_T), intent(out)::val
        end subroutine c_get_global_int
    end interface
    
    interface
        ! return double variable name to val
        subroutine c_get_global_double(name, val) bind(c)
            use iso_c_binding
            implicit none
            character(kind=c_char), intent(in)::name(*)
            real(kind=C_DOUBLE), intent(out)::val
        end subroutine c_get_global_double
    end interface
    
    interface
        ! return bool variable name to val
        subroutine c_get_global_bool(name, val) bind(c)
            use iso_c_binding
            implicit none
            character(kind=c_char), intent(in)::name(*)
            integer(kind=C_INT64_T), intent(out)::val
        end subroutine c_get_global_bool
    end interface
    
    interface
        ! return string variable name to val
        subroutine c_get_global_string(name, val, length) bind(c)
            use iso_c_binding
            implicit none
            character(kind=c_char), intent(in)::name(*)
            character(kind=c_char,len=1), intent(out)::val
            integer(kind=C_INT64_T), intent(out)::length
        end subroutine c_get_global_string
    end interface
        
    interface
        ! return integer vector variable name, vector number num, to val
        subroutine c_get_global_int_vec(name, num, val) bind(c)
            use iso_c_binding
            implicit none
            character(kind=c_char), intent(in)::name(*)
            integer(kind=C_INT64_T), intent(in)::num
            integer(kind=C_INT64_T), intent(out)::val(1)
        end subroutine c_get_global_int_vec
    end interface
    
contains
    !extract the global variables
    subroutine grab_globals() bind(c)
        use iso_c_binding
        implicit none
        
        integer(kind=8)::i,j,n2
        integer(kind=8)::ref_singles = 0
        
        ! set global variable values
        call get_global_int("num_orbitals", num_orbitals)
        call get_global_int("num_orbitalsC2", num_orbitalsC2)
        call get_global_double("energy_tol", energy_tol)
        call get_global_double("norm_tol", norm_tol)
        call get_global_int("valence_ci_flag", valence_ci_flag)
        call get_global_bool("restart_flag", restart_flag)
        call get_global_bool("integralDirect", integralDirect)
        call get_global_bool("fullyIntegralDirect", fullyIntegralDirect)
        call get_global_bool("cdVecsInMemory", cdVecsInMemory)
        call get_global_bool("directFourInternal", directFourInternal)
        call get_global_bool("directLowMem", directLowMem)
        call get_global_bool("directSuperLowMem", directSuperLowMem)
        call get_global_int("number_bas", number_bas)
        call get_global_int("number_basC2", number_basC2)
        call get_global_int("num_frozen", num_frozen)
        call get_global_int("num_inactive", num_inactive)
        call get_global_int("num_active", num_active)
        call get_global_int("spinM", spinM)
        call get_global_int("num_elec", num_elec)
        call get_global_int("num_ref", num_ref)
        call get_global_int("nat_orb_flag", nat_orb_flag)
        call get_global_string("scratch_directory", scratch_directory)
        call get_global_double("integral_threshold", integral_threshold)
        call get_global_double("ao_integral_threshold", ao_integral_threshold)
        call get_global_double("internal_threshold", internal_threshold)
        call get_global_double("wp_default_radius", wp_default_radius)
        call get_global_double("tov_occupied_default_radius", tov_occupied_default_radius)
        call get_global_double("tov_occupied_multiplier", tov_occupied_multiplier)
        call get_global_double("wp_multiplier", wp_multiplier)
        call get_global_double("virtual_threshold", virtual_threshold)
        call get_global_double("tov_virtual_default_radius", tov_virtual_default_radius)
        call get_global_bool("sphereprint", sphereprint)
        call get_global_double("tov_virtual_multiplier", tov_virtual_multiplier)
        call get_global_double("tov_cylinder_radius", tov_cylinder_radius)
        call get_global_int("reference_ci_flag", reference_ci_flag)
        call get_global_double("wp_cylinder_radius", wp_cylinder_radius)
        call get_global_int("acpf_flag", acpf_flag)
        call get_global_double("cd_thresh", cd_thresh)
        call get_global_mem("max_mem", max_mem)
        max_mem = max_mem / 8
        call get_global_mem("max_mem_ints", max_mem_ints)
        max_mem_ints = max_mem_ints / 8
        call get_global_int("nonlocal_flag", nonlocal_flag)
        call get_global_int("num_roots", num_roots)
        call get_global_bool("restart_from_old_CD", restart_from_old_CD)
        call get_global_bool("acpf_root_follow", acpf_root_follow)
        call get_global_bool("ACPF_ROOT_FOLLOW_HELP", ACPF_ROOT_FOLLOW_HELP)
        call get_global_double("custom_g_val", custom_g_val)
#ifdef _OPENMP
        call get_global_int("numThreads", numberOfThreads)
#endif
        call get_global_int("blockSizeSigma", blockSizeSigma)
        call get_global_int("blockSizeCI", blockSizeCI)
        call get_global_mem("for_buf_blocksizeInts", for_buf_blocksizeInts)
        for_buf_blocksizeInts = for_buf_blocksizeInts / 1024
        call get_global_mem("for_buf_maxmemInts", for_buf_maxmemInts)
        call get_global_bool("for_buf_storeIntegrals", for_buf_storeIntegrals)
        call get_global_mem("for_buf_maxmemCD", for_buf_maxmemCD)
        call get_global_mem("for_buf_maxmemAIJK", for_buf_maxmemAIJK)
        call get_global_bool("for_buf_storeCDVecs", for_buf_storeCDVecs)
        call get_global_mem("for_buf_maxmemSegs", for_buf_maxmemSegs)
        call get_global_int("for_buf_blocksizeSegs", for_buf_blocksizeSegs)
        call get_global_double("ASSUMED_DAV_DISK_SPACE", ASSUMED_DAV_DISK_SPACE)
        call get_global_int_vecs("references", references, num_ref, num_orbitals)
        call get_global_bool("DENSITY FITTING", density_fitting)
        call get_global_int("num_internal", num_internal)
        call get_global_bool("CPP_DECOMPOSED_INTS",CPP_DECOMPOSED_INTS)
        call get_global_bool("CPP_TRANSFORMED_INTS",CPP_TRANSFORMED_INTS)

        write(*,*) "Finished importing global variables"
        
        
        !// Miminum size for 'integral_buffer' (just hold diagonal integrals)
        n2 = num_orbitals-num_internal
        n2 = n2*num_orbitals
        integral_buffer_size = n2
        
        !// DETERMINE MAXIMUM NUMBER OF OPEN SHELLS
          open_shells = 0
          do i = 1, num_ref
              
            ref_singles = 0
            do j = num_inactive+1, num_inactive+num_active
                  if (references(j,i)==1) ref_singles = ref_singles+1
              enddo
              
            open_shells = max(ref_singles+4,open_shells)
            
        enddo
        
        !// ADJUST OPEN SHELLS IF WE DON'T HAVE ENOUGH ELECTRONS
        if (num_elec < open_shells) open_shells = num_elec
        
        !// SET A FEW VARIABLES WHICH ARE FULLY DETERMINED AT THIS POINT
        num_external = num_orbitals - num_internal
        num_extC2 = num_external*(num_external-1)/2

        ! Here we echo all the possible input parameters back to the user
        ! this will most likely include some parameters which were never
        ! specified in the user input ... this is deliberate !
        write(ioOutput,130)
        130 format(/,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/&
                    ,1x,"!//                              ",/&
                    ,1x,"!// TIGER CI IS USING THE FOLLOWING INPUT PARAMETERS         ",/&
                    ,1x,"!//                              ",/&
                    ,1x,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",/)

        ! ALL of the possible input parameters
        write(ioOutput,*) "NUMBER OF ORBITALS               " , num_orbitals
        write(ioOutput,*) "ENERGY TOLERANCE                 " , energy_tol
        write(ioOutput,*) "RESIDUAL TOLERANCE               " , norm_tol
        write(ioOutput,*) "VALENCE CI FLAG                  " , int_to_bool(valence_ci_flag)
        write(ioOutput,*) "RESTART FLAG                     " , restart_flag
        write(ioOutput,*) "NUMBER OF BASIS FUNCTIONS        " , number_bas
        write(ioOutput,*) "NUMBER OF FROZEN ORBITALS        " , num_frozen
        write(ioOutput,*) "NUMBER OF INACTIVE ORBITALS      " , num_inactive
        write(ioOutput,*) "NUMBER OF ACTIVE ORBITALS        " , num_active
        write(ioOutput,*) "SPIN MULTIPLICITY                " , spinM
        write(ioOutput,*) "NUMBER OF ELECTRONS              " , num_elec
        write(ioOutput,*) "NUMBER OF REFERENCES             " , num_ref
        write(ioOutput,*) "NATURAL ORBITAL FLAG             " , int_to_bool(nat_orb_flag)
        write(ioOutput,*) "INTEGRAL THRESHOLD               " , integral_threshold
        write(ioOutput,*) "AO INTEGRAL THRESHOLD            " , ao_integral_threshold
        write(ioOutput,*) "WP DEFAULT RADIUS                " , wp_default_radius
        write(ioOutput,*) "TOV OCCUPIED DEFAULT RADIUS      " , tov_occupied_default_radius
        write(ioOutput,*) "TOV OCCUPIED RADIUS MULTIPLIER   " , tov_occupied_multiplier
        write(ioOutput,*) "WP RADIUS MULTIPLIER             " , wp_multiplier
        write(ioOutput,*) "VIRTUAL OCCUPATION THRESHOLD     " , virtual_threshold
        write(ioOutput,*) "TOV VIRTUAL DEFAULT RADIUS       " , tov_virtual_default_radius
        write(ioOutput,*) "TOV VIRTUAL RADIUS MULTIPLIER    " , tov_virtual_multiplier
        write(ioOutput,*) "TOV CYLINDER RADIUS              " , tov_cylinder_radius
        write(ioOutput,*) "REFERENCE CI FLAG                " , int_to_bool(reference_ci_flag)
        write(ioOutput,*) "CYLINDER RADIUS                  " , wp_cylinder_radius
        write(ioOutput,*) "ACPF FLAG                        " , acpf_flag
        write(ioOutput,*) "CD THRESHOLD                     " , cd_thresh
        write(ioOutput,*) "CD BUFFER                        " , max_mem
        write(ioOutput,*) "MAX MEM SEGS                     " , for_buf_maxmemSegs
        write(ioOutput,*) "INT MEMORY                       " , max_mem_ints
        write(ioOutput,*) "NONLOCAL                         " , int_to_bool(nonlocal_flag)
        write(ioOutput,*) "NUM_ROOTS                        " , num_roots
        !write(ioOutput,*) "LargeCholeskyVecDistanceCutOff   " , LargeCholeskyVecDistanceCutOff
        write(ioOutput,*) "Sphere Print                     " , sphereprint
        write(ioOutput,*) "CD RESTART                       " , restart_from_old_CD
        write(ioOutput,*) "USE ACPF ROOT FOLLOW             " , ACPF_ROOT_FOLLOW
        write(ioOutput,*) "USE ACPF REF ENERGY              " ,  use_input_ref_energy
        if(use_input_ref_energy) write(ioOutput,*) "REF ENERGY              " ,  user_ref_energy
        write(ioOutput,*) "SCRATCH DIRECTORY                " , scratch_directory
        write(ioOutput,*) "INTEGRAL DIRECT MODE             " , integralDirect
        write(ioOutput,*) "FULLY INTEGRAL DIRECT MODE       " , fullyIntegralDirect
        write(ioOutput,*) "LOW MEMORY DIRECT                " , directLowMem
        write(ioOutput,*) "CD VECTORS IN MEMORY             " , cdVecsInMemory
        write(ioOutput,*) "CPP_DECOMPOSED_INTS              " , CPP_DECOMPOSED_INTS
        write(ioOutput,*) "CPP_TRANSFORMED_INTS             " , CPP_TRANSFORMED_INTS

        ! output info about the references
        write(ioOutput,150) num_ref
        150 format(1x,"Number of references...........",i5)

        do i = 1, num_ref
            write(ioOutput,160) i,(references(j,i), j = 1, num_orbitals)
        enddo
        write(ioOutput,*)
        160 format(6x,"Reference #",i5,/,&
                   6x,700i2,/)

        write(*,*)
        write(*,*) "*******************************"
        write(*,*) "Calling Tiger CI..."
        write(*,*) "*******************************"
        write(*,*)
        
        call tiger_ci
        
    end subroutine grab_globals
        
    ! wrapper to call c function for extraction of int variable name to var
    subroutine get_global_int(name, var)
        use iso_c_binding
        implicit none
        
        character(len=*), intent(in)::name
        integer(kind=8), intent(out)::var
        
        call c_get_global_int(name//C_NULL_CHAR, var)
    end subroutine get_global_int
    
    ! wrapper to call c function for extraction of int variable name to var
    subroutine get_global_mem(name, var)
        use iso_c_binding
        implicit none
        
        character(len=*), intent(in)::name
        integer(kind=8), intent(out)::var
        
        call c_get_global_int(name//C_NULL_CHAR, var)
        
        var = var * 1024 * 1024
    end subroutine get_global_mem
    
    ! wrapper to call c function for extraction of double variable name to var
    subroutine get_global_double(name, var)
        use iso_c_binding
        implicit none
        
        character(len=*), intent(in)::name
        real(kind=8), intent(out)::var
        
        call c_get_global_double(name//C_NULL_CHAR, var)
    end subroutine get_global_double
    
    ! wrapper to call c function for extraction of bool variable name to var
    subroutine get_global_bool(name, var)
        use iso_c_binding
        implicit none
        
        character(len=*), intent(in)::name
        logical, intent(out)::var
        integer(kind=8)::c_var
        
        call c_get_global_bool(name//C_NULL_CHAR, c_var)
        
        if (c_var == 0) then
            var = .false.
        else
            var = .true.
        endif
    end subroutine get_global_bool
    
    ! wrapper to call c function for extraction of string variable name to var
    subroutine get_global_string(name, string)
        use iso_c_binding
        implicit none
        
        character(:), allocatable::string
        character(len=*), intent(in)::name
        
        integer(kind=8)::length
        character(len=100)::pass_string
        
        call c_get_global_string(name//C_NULL_CHAR, pass_string, length)
        string = pass_string(1:length)
        
    end subroutine get_global_string
    
    ! wrapper to call c function for extraction of vector of integer vectors variable name to var, given dimensions num_vecs and num_ints
    subroutine get_global_int_vecs(name, vecs, num_vecs, num_ints)
        use iso_c_binding
        implicit none
        
        character(len=*), intent(in)::name
        integer(kind=8), allocatable, intent(out)::vecs(:,:)
        integer(kind=8), intent(in)::num_vecs, num_ints

        integer(kind=8), allocatable :: vec(:)
        integer(kind=8)::i
        
        allocate(vecs(num_ints, num_vecs))
        
        ! call extraction function once for each of the vectors in vector of vectors
        do i = 1, num_vecs
            allocate(vec(num_ints))
            call c_get_global_int_vec(name//C_NULL_CHAR, i-1, vec)
            vecs(:,i) = vec(:)
            deallocate(vec)            
        end do
        
    end subroutine get_global_int_vecs

    ! new subroutine to check I/O return (much like the allocatecheck !)
    !> \brief checks iostat for errors during a read statement
    subroutine read_check(ierr)
        implicit none
        integer, intent(in) :: ierr
        if (ierr == -1) then
            write(*,*) "Fatal Error when reading input file: End Of File"
            call flush(6)
            stop
        else if ( ierr == -2 ) then
            write(*,*) "Fatal Error when reading input file: End Of Record"
            call flush(6)
            stop
        end if
    end subroutine read_check

    !> \brief converts an integer to boolean , 0 false everything else true, for a cleaner output
    logical function int_to_bool (i)
        implicit none
        integer, intent(in) :: i
        if(i == 0) then
            int_to_bool = .false.
            return
        end if
        int_to_bool = .true.
        return
    end function int_to_bool
end module fetch_globals
