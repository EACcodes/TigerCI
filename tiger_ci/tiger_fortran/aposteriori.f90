! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> \brief Handles the calculation of the a posteriori corrections for size extensivity                                 
!>        after an (MR)SDCI calculation.                                                                               
!>                                                                                                                 
!> The following references were used   
!>  - Duch, W. & Diercksen, G. H. F. Size-extensivity corrections in configuration interaction methods. J. Chem. Phys. 101, 3018 (1994)
!> - Szalay, P. G., Müller, T., Gidofalvi, G., Lischka, H. & Shepard, R. Multiconfiguration self-consistent field and multireference configuration interaction methods and applications. Chem. Rev. 112, 108–81 (2012)
!> - Szalay, P. G. in Encycl. Comput. Chem. 1–47 (Wiley Online Library, 2005). at <http://onlinelibrary.wiley.com/doi/10.1002/0470845015.cn0066/full>
!> (Take care with the last one ... the preprint we recieved did have a typo or two in the formulas!) 
module aposteriori_mod
  use global_var_mod
  use molecule_var_mod
  use new_tree_search_mod  
  use initial_davidson_struct
  implicit none
  private
  public :: calculate_a_posteriori_corrections
contains

  !> \brief Calculates a variety of different a posteriori size extensivity corrections (for SDCI wavefuncs)
  !> \param wavefunc  The converged wavefunction
  !> \param guess     The initial guess for the wavefunction (containing the zeroth order wavefunction)
  !> \param total_energy  The total converged energy
  !> \param reference_energy  The energy of the zeroth order wavefunction (not including dynamic correlation)
  subroutine calculate_a_posteriori_corrections(wavefunc, guess, total_energy, reference_energy)
    ! inputs
    real(kind=real8), dimension(:) :: wavefunc
    real(kind=real8), intent(in) :: total_energy, reference_energy
    type(initial_dav_guess) :: guess

    ! local variables
    real(kind=real8) :: correlation_energy, C0squared, C0squared_alt
    integer :: num_corr_elec
    character(len=50) :: title, def1, def2, line


    ! error check 
    if ( num_elec < 4 ) then 
       write(*,*) "You have a really small number of electrons in your system"
       write(*,*) "To the point where (most) of the a posteriori corrections"
       write(*,*) "make absolutely no sense."
       write(*,*) 
       return 
    end if

    ! get our basic constants
    correlation_energy = total_energy - reference_energy 
    num_corr_elec =  num_elec - (2 * num_frozen)

    ! figure out the value of C0squared
    ! just c_0^2 for SDCI
    ! its the more complicated sum c_{ref}^2 for MRSDCI
    C0squared = calculate_ref_coeff_squared(wavefunc)

    ! an alternative formula uses the coefficient of each reference in the converged (MR)SDCI wavefunction
    ! multiplied by the coefficient of that configuration in the zeroth order wavefunction
    C0squared_alt = calculate_ref_coeff_squared_alt(wavefunc, guess)

    ! Now just print out the corrections


    write(title,*)  "(/,9x,'A Posteriori Corrections',9x,/,)"
    write(def1,*)   "(3x,'Define c_0^2 = \sum c_{ref}^2 = ', E15.8)"
    write(def2,*)   "(3x,'Define c_0^2 =  <ref|ci>^2    = ', E15.8)"
    write(line,*)   "(3x,A22,3x,E15.8)"

    write(*,title)
    write(*,def1)  C0squared
    write(*,line) "Davidson              ", davidson(C0squared, correlation_energy)
    write(*,line) "Renormalized Davidson ", renorm_davidson(C0squared, correlation_energy)
    write(*,line) "Davidson-Silver       ", davidson_silver(C0squared, correlation_energy)
    write(*,line) "Meissner              ", meissner(C0squared, correlation_energy, num_corr_elec)
    write(*,line) "Pople                 ",  pople(C0squared, correlation_energy, num_corr_elec)
    write(*,*) 

    write(*,def2) C0squared_alt
    write(*,line) "Davidson              ", davidson(C0squared_alt, correlation_energy)
    write(*,line) "Renormalized Davidson ", renorm_davidson(C0squared_alt, correlation_energy)
    write(*,line) "Davidson-Silver       ", davidson_silver(C0squared_alt, correlation_energy)
    write(*,line) "Meissner              ", meissner(C0squared_alt, correlation_energy, num_corr_elec)
    write(*,line) "Pople                 ",  pople(C0squared_alt, correlation_energy, num_corr_elec)
    write(*,*) 

  end subroutine calculate_a_posteriori_corrections

  !> \brief Calculates \f$ \langle \psi_0 | \psi_{mrsdci} \rangle \f$ where \f$\psi_0\f$ is the zeroth order wavefunction
  real(kind=real8) function calculate_ref_coeff_squared_alt(wavefunc, guess)
    ! inputs
    real(kind=real8), dimension(:) :: wavefunc
    type(initial_dav_guess) :: guess

    ! local variables
    integer :: weight, address, nSpinFunc
    integer :: spin, i 
    integer :: ierr
    real(kind=real8) :: total 
    type(orbital_path) :: path 
    type(graph_search_state) :: graph
    real(kind=real8), dimension(:), allocatable :: zeroth_order

    ! build the zeroth order wavefunction
    allocate(zeroth_order(total_csfs), stat=ierr)
    call allocatecheck(ierr, "zeroth order wavefunction in aposteriori")
    zeroth_order = 0.0
    call pull_initial_guess(guess, zeroth_order, 1)

    ! loop over all the configurations in the wavefunction with no electrons in the virtual space
    call allocate_orbital_path(path)
    total = 0.0 
    call init_tree_search(graph, 0, 0, num_internal, num_elec)
    do while ( get_internal_path(path, graph))
       weight = sum(path%arc_weights(1:num_internal)) + 1
       address = internal_index_vector0(weight)
       if ( address < 0 ) cycle
       if (.not. is_reference(path%occupations(1:num_internal))) cycle
       nSpinFunc = fsn(path%num_singles)
       do spin = 1, nSpinFunc
          i = address + (spin-1)
          total = total + wavefunc(i)*zeroth_order(i)          
       end do
    end do
    calculate_ref_coeff_squared_alt = total * total            ! we square the sum at the end
    deallocate(zeroth_order)
  end function calculate_ref_coeff_squared_alt




  !> \brief Calculates \f$ \sum c_{ref}^2\f$ a necessary input for the a posteriori corrections
  !> \param wavefunc  The converged wavefunction
  real(kind=real8) function calculate_ref_coeff_squared(wavefunc)
    ! inputs
    real(kind=real8), dimension(:) :: wavefunc

    ! local variables
    integer :: weight, address, nSpinFunc
    integer :: spin, i 
    real(kind=real8) :: total 
    type(orbital_path) :: path 
    type(graph_search_state) :: graph

    ! loop over all the configurations in the wavefunction with no electrons in the virtual space
    call allocate_orbital_path(path)
    total = 0.0 
    call init_tree_search(graph, 0, 0, num_internal, num_elec)
    do while ( get_internal_path(path, graph))
       weight = sum(path%arc_weights(1:num_internal)) + 1
       address = internal_index_vector0(weight)
       if ( address < 0 ) cycle
       if (.not. is_reference(path%occupations(1:num_internal))) cycle
       nSpinFunc = fsn(path%num_singles)
       do spin = 1, nSpinFunc
          i = address + (spin-1)
          total = total + wavefunc(i)*wavefunc(i)
       end do
    end do
    calculate_ref_coeff_squared = total 
  end function calculate_ref_coeff_squared

  !> \brief Davidson's original correction
  !> \param C0squared  \f$c_0^2\f$
  !> \param corrE      The converged correlation energy
  real(kind=real8) function davidson( C0squared, corrE)
    real(kind=real8), intent(in) :: C0squared, corrE
    davidson = corrE * ( 1 - C0squared ) 
  end function davidson

  !> \brief Davidson's renormalized correction
  !> \param C0squared  \f$c_0^2\f$
  !> \param corrE      The converged correlation energy
  real(kind=real8) function renorm_davidson( C0squared, corrE)
    real(kind=real8), intent(in) :: C0squared, corrE
    renorm_davidson = corrE * ( (1 - C0squared) / C0squared)
  end function renorm_davidson

  !> \brief The Davidson-Silver correction
  !> \param C0squared  \f$c_0^2\f$
  !> \param corrE      The converged correlation energy
  real(kind=real8) function davidson_silver( C0squared, corrE)
    real(kind=real8), intent(in) :: C0squared, corrE
    davidson_silver = corrE * ( (1 - C0squared ) / ( 2 * C0squared - 1 )) 
  end function davidson_silver

  !> \brief Meissner's correction
  !> \param C0squared  \f$c_0^2\f$
  !> \param corrE      The converged correlation energy
  !> \param N          Number of correlating electrons
  real(kind=real8) function meissner( C0squared, corrE, N)
    integer, intent(in) :: N
    real(kind=real8), intent(in) :: C0squared, corrE
    real(kind=real8) :: term1, term2 
    term1 = (1 - C0squared) / C0squared
    term2 = ( (N-2)*(N-3) ) / real(N*(N-1),real8)
    meissner = corrE * term1 * term2
  end function meissner

  !> \brief Pople's correction
  !> \param C0squared  \f$c_0^2\f$
  !> \param corrE      The converged correlation energy
  !> \param N          Number of correlating electrons
  real(kind=real8) function pople( C0squared, corrE, N)
    real(kind=real8), intent(in) :: C0squared, corrE
    integer, intent(in) :: N
    real(kind=real8) :: term1, term2 ,term3
!    theta = ACOS(C0squared)
!    term1 = sqrt(N*N + 2*N*(TAN(2*theta))**2) - N
!    term2 = 2 * ( 1/COS(2*theta) -1 )
!    pople = corrE * (term1/term2 - 1)
    
    ! using the formulation from:
    ! Duch, W., & Diercksen, G. H. F. (1994).  The Journal of Chemical Physics, 101(4), 3018. doi:10.1063/1.467615
    ! as my implementation with trig functions seems a bit buggy
    term1 = 2*C0squared 
    term2 = (2*C0squared -1)
    term3 = 1 + 8*C0squared * (1-C0squared) / (N*term2*term2)
    pople = term1 / (term2 * (1 + SQRT(term3))) - 1
    pople = corrE * pople  
  end function pople

end module aposteriori_mod
