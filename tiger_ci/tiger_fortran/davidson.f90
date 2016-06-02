! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

!> The new consolidated davidson solvers under the "Great Unification"
!>
!> Some of these routines get rather confusing. We try to reuse memory
!> as much as possible as these vectors can get rather large. To try and 
!> keep things clear, when any subroutine uses the civec or sigmavec as 
!> scratch space (instead of for their actual values) I pass them by 
!> keyword argument ( for instance: call something(scratch=civec))
!> Hopefully this notation will make it clearer when we are reusing
!> memory.


module davidson_mod
  use utilities_mod
  use global_var_mod
  use post_iteration_mod
  use math_utils
  use initial_davidson_struct
  use davidson_io_mod
  use SGGA
  use gtimesc_mod
  use blocked_locks_mod
  use sgga_struct_mod
  use graph_var_mod
  use fortran_timing
  use davidson_restart
  implicit none

  private
  public :: davidson_liu, davidson_liu_gen_eig

contains

  !> \brief Uses the Davidson-Liu (block Davidson) algorithm to solve for the lowest few eigenvalues/vectors of Hc=Ec
  !> \param civec       The c vector (a blocked lock vector type)
  !> \param sigmavec    The Hc vector (a blocked lock vector type)
  !> \param roots      The number of eigenvalues/eigenvectors to solve for
  !> \param energy On exit the energies for the desired eigenvalues
  !> \param guess The initial guess for the solver     
  !> \param sgga_data   Scratch spaces for the SGGA (matrix vector multiplication routine)
  !> \param mult      A pointer to the function which performs the H*c multiplication
  !>
  !> This routines assumes the diagonal elements of the Hamiltonian are already computed and stored in ioDiagStore. 
  subroutine davidson_liu(civec, sigmavec, roots, energy, guess, sgga_data, mult)

    ! inputs
    type(initial_dav_guess) :: guess
    real(real8) :: energy(:)                               
    type(sgga_struct) :: sgga_data
    integer :: roots
    type(blockedLockVectorType) :: civec, sigmavec
    procedure(htimesc), pointer :: mult

    ! local variables
    integer :: ierr
    integer :: i,iRoot
    integer :: maxit
    integer :: iteration
    integer :: num_vec, this_vec
    integer :: index, subspace_size

    logical :: halt

    real(kind=real8) :: norm, time

    real(kind=real8), allocatable, dimension(:) :: eig_vals, energy_change, scratch
    real(kind=real8), allocatable, dimension(:,:) :: h_subspace, h_subspace_copy, eig_vecs, previous_eig_vecs

    integer, dimension(:), allocatable :: root2index
    logical, dimension(:), allocatable :: converged, restarts

    type(clock) :: iteration_time

    write(*,*) " :: Solving a standard eigenvalue problem for one or more roots with the Davidson-Liu algorithm :: "

    ! memory allocations
    maxit = 50 * roots
    energy = 0
    allocate(h_subspace(maxit,maxit), h_subspace_copy(maxit,maxit), eig_vecs(maxit, roots), previous_eig_vecs(maxit, roots), & 
         eig_vals(maxit), energy_change(maxit), scratch(total_csfs), stat=ierr)
    call allocatecheck(ierr, "subspaces in davidson_liu")
    h_subspace = 0
    h_subspace_copy = 0
    eig_vecs = 0
    eig_vals = 0
    energy_change = 0
    scratch = 0

    allocate(root2index(roots), converged(roots), restarts(maxit), stat =ierr)
    call allocatecheck(ierr, "extra stuff for davidson-liu")
    do iRoot = 1, roots
       root2index(iRoot) = iRoot
    end do
    restarts = .false.
    converged = .false.

    ! initial guesses
    ! if you are restarting the I/O is a bit excessive .. oh well its not a serious bottleneck
    do iRoot = 1, roots
       call pull_initial_guess(guess, civec%v, iRoot)
       call davidson_io_write(civec%v, cfile, iRoot)
    end do
    num_vec = roots

    ! Iterative Eigenvalue Solver
    ! ***************************
    do i = 1, maxit
       iteration = i 
       call start_clock(iteration_time)

       do iRoot = 1, roots
          ! skip if we already converged this root
          if (converged(iRoot)) cycle

          !where on disk I need to look for the vector of this iteration(i) associated with this root (iRoot)
          ! this is also the (absolute) number for this vector
          index = root2index(iRoot)

          ! read in the current vector
          call davidson_io_read(civec%v, cfile, index)

          ! Calculate the new H*c product
          sigmavec%v = 0.0
          call mult(civec, sigmavec, sgga_data)
          call davidson_io_write(sigmavec%v, Hcfile, index)

          ! Expand the subspace with the new vector
          call expand_krylov_subspace(h_subspace, sigmavec%v, index, scratch=civec%v)
       end do ! loop over roots

       ! Now we solve the small CI problem on the subspace
       ! { h_subspace } eig_vec = eig_val * eig_vec
       !
       ! Note that the lapack routine will overwrite the matrix with the eigenvectors
       ! so we make a copy of the matrix
       h_subspace_copy(:,:) = h_subspace
       subspace_size = num_vec
       previous_eig_vecs = eig_vecs
       call dsyev_wrapper(.false., .true., subspace_size, h_subspace_copy, maxit, eig_vals) 

       energy_change(1:roots) = abs( eig_vals(1:roots) - energy(1:roots))
       energy(1:roots) =  eig_vals(1:roots)
       eig_vecs(1:num_vec,1:roots) = h_subspace_copy(1:num_vec,1:roots)
       
       if (need_to_restart_subspace(num_vec, roots, total_csfs))then
           write(*,*) " "
           write(*,*) "  RESTARTING THE DAVIDSON SUBSPACE "
           write(*,*) " " 
           restarts(i) = .true.
           call restart_davidson_subspace(roots, converged, num_vec, eig_vecs, previous_eig_vecs, civec%v, sigmavec%v, h_subspace)
       end if

       write(*,9030) i 
       time = get_clock_wall_time(iteration_time)

       do iRoot = 1, roots
          ! skip if we already converged this root
          if (converged(iRoot)) cycle

          ! where on disk I need to look for the vector of this iteration(i) associated with this root (iRoot)
          ! this is also the (absolute) number for this vector
          this_vec = root2index(iRoot)

          ! build the new guess for the eigenvector by writing the eigenvector alpha 
          ! in the space of all configurations
          call build_wavefunction(eig_vecs(:,iRoot), civec%v, subspace_size, scratch=sigmavec%v)
          
          ! run post-iteration checks
          call post_iteration_driver(civec%v, iRoot)

          ! build the residual
          call build_residual_std_eig(civec%v, energy(iRoot), eig_vecs(:,iRoot), subspace_size, &
               Hc_scratch=sigmavec%v, c_scratch=scratch)
          norm = sqrt(dot_product(civec%v,civec%v))

          ! produce the output for this iteration/root
          write(iooutput,9010)
          write(iooutput,9020) i,iRoot,energy(iRoot),norm,energy_change(iRoot),time
          call write_post_iteration_analysis 
          call flush(iooutput)

          ! check convergence
          halt = is_converged(energy_change(iRoot), norm, norm_tol, energy_tol, i , restarts, maxit)
          if (halt) then
             converged(iRoot) = .true. 
             write(*,*) " "
             write(*,*) "CONVERGED ROOT = " , iRoot
             write(*,*) " "
             cycle
          end if

          ! if we get here we aren't converged ... build the next vector by preconditioning the residual
          call build_next_vector_std_eig(civec%v, energy(iRoot), diag_scratch=sigmavec%v)

          ! Orthogonalize against previous vectors then normalize it
          ! Note that I need to orthogonalize against all the vectors, including the ones built for the other
          ! roots this iteration
          call orthonormalize_vector(civec%v, num_vec, scratch=sigmavec%v)

          ! write out the new expansion vector to disk
          call davidson_io_write(civec%v, cfile, num_vec+1)
          num_vec = num_vec+1
          root2index(iRoot) = num_vec
       end do ! loop over roots

       ! check convergence over all roots
       halt = .true.
       do iRoot = 1, roots
          if ( .not. converged(iRoot) ) then
             halt =.false.
          end if
       end do
       if (halt) exit

    end do ! davidson iterations

    ! produce final output for each root
    do iRoot = 1, roots
       call build_wavefunction(eig_vecs(:,iRoot),  civec%v, num_vec, scratch=sigmavec%v)
       call post_iteration_driver(civec%v, iRoot)
       write(*,*) "   Final wavefunction output for root = ", iRoot
       call write_post_iteration_analysis(.true.)
    end do

    ! memory clean up
    deallocate(h_subspace, h_subspace_copy, eig_vals, stat=ierr)
    call deallocatecheck(ierr, "subspaces in expt_acpf_comp_kernel")
    deallocate(root2index, converged, energy_change, stat=ierr)
    call deallocatecheck(ierr, "extras in gen_eig_fake_multi_state_kernel")

9010 format(/,"   -------  --------   ------------------    ---------------    -------------   ----------",/,&
         "    Iter      Root           Energy             |Hc-Ec|            DeltaE        wall time",/,&
         "   -------  --------   ------------------    ---------------    -------------   ----------")
9020 format('  ',i5,5x,i5,5x,d19.12,3x,d15.8,3x,d15.8,1x,f10.2,/)
9030 format(/,'  -----------------------------',/'   Iteration ', i3, ' has finished ',/,'  -----------------------------')
  end subroutine davidson_liu

  !> \brief Solves for multiple roots of the generalized eigenvalue problem [H-E_ref]c=EGc where H and G are matrices and E_ref a scalar
  !> \param c    Scratch space for c
  !> \param Hc Scratch space for storing Hc products
  !> \param Gc     Scratch space for storing Gc products
  !> \param work     An additional scratch space
  !> \param roots                   The number of eigenvalues/eigenvectors to solve for
  !> \param energy               On exit the energies for the desired eigenvalues
  !> \param ref_energy   Energy of the reference or 0th order wavefunction
  !> \param guess                   The initial guess for the solver 
  !> \param sgga_data               SGGA scratch data
  !> \param hmultc   Pointer to the routine that performs the H*c operation
  !> \param gmultc   Pointer to the routine that performs the G*c operation
  subroutine davidson_liu_gen_eig(c, Hc, Gc, work, roots, energy, ref_energy, guess, sgga_data, hmultc, gmultc)
    ! inputs
    type(initial_dav_guess) :: guess
    real(real8) :: energy(:)                               
    type(sgga_struct) :: sgga_data
    integer :: roots
    type(blockedLockVectorType) :: Hc, c
    real(kind=real8), dimension(:) :: Gc, work
    procedure(htimesc), pointer :: hmultc
    procedure(gtimesc), pointer :: gmultc
    real(kind=real8) :: ref_energy

    ! local variables 
    integer :: ierr
    integer :: i,iRoot,r
    integer :: maxit
    integer :: iteration
    integer :: num_vec, this_vec
    integer :: subspace_size

    logical :: halt

    real(kind=real8) :: norm, time

    real(kind=real8), dimension(:), allocatable :: eig_vals, energy_change, alpha
    real(kind=real8), dimension(:,:), allocatable :: h_sub, g_sub, h_sub_copy, g_sub_copy, eig_vecs, previous_eig_vecs

    integer, dimension(:), allocatable :: root2index
    logical, dimension(:), allocatable :: converged, restarts

    type(clock) :: iteration_time


    write(*,*) " :: Solving a generalized eigenvalue equation via Davidson-Liu :: "    

    ! memory allocations
    maxit = 50 * roots
    allocate(h_sub(maxit,maxit), h_sub_copy(maxit,maxit), g_sub(maxit,maxit), g_sub_copy(maxit,maxit), &
         eig_vals(maxit), eig_vecs(maxit, roots), previous_eig_vecs(maxit, roots), alpha(maxit), stat=ierr)
    call allocatecheck(ierr, "subspaces in expt_acpf_comp_kernel")
    h_sub = 0.0
    h_sub_copy = 0.0
    g_sub_copy = 0.0
    g_sub = 0.0 
    eig_vals = 0.0
    eig_vecs = 0.0 

    allocate(root2index(roots), converged(roots), energy_change(roots), restarts(maxit), stat=ierr)
    call allocatecheck(ierr, "extra stuff for excited state execution")
    root2index = 0
    converged = .false.
    restarts = .false.
    energy_change = 0 

    do i = 1, roots 
       root2index(i) = i 
    end do
    converged = .false.
    energy(1:roots) = ref_energy

    ! write each of my initial guesses to disk
    do iRoot = 1, roots
       call pull_initial_guess(guess, c%v, iRoot)
       call davidson_io_write(c%v, cfile, iRoot)
    end do
    num_vec = roots

    ! Iterative Eigenvalue Solver
    ! ***************************
    do i = 1, maxit
       iteration = i 
       call start_clock(iteration_time)

       do iRoot = 1, roots
          ! skip if we already converged this root
          if (converged(iRoot)) cycle

          !where on disk I need to look for the vector of this iteration(i) associated with this root (iRoot)
          ! this is also the (absolute) number for this vector
          this_vec = root2index(iRoot)

          ! read in the current vector
          call davidson_io_read(c%v, cfile, this_vec)

          ! Calculate the new (H-E_ref)*c product
          Hc%v = 0.0
          call HMultC(c, Hc, sgga_data)
          Hc%v(:) = Hc%v(:) - ref_energy * c%v(:)
          call davidson_io_write(Hc%v, Hcfile, this_vec)

          ! Calculate G*c
          Gc(:)=c%v(:)
          call GMultC(Gc)

          ! Expand the subspace with the new vector
          call expand_krylov_subspace(h_sub, Hc%v, this_vec, scratch=c%v)
          call expand_krylov_subspace(g_sub, Gc, this_vec, scratch=c%v)

       end do ! loop over roots

       ! solve the small generalized eigenvalue problem
       ! h_sub * eig_vecs = lambda * g_sub * eig_vecs
       ! 
       ! h_sub and g_sub are my operators on the space spanned by the davidson vectors
       ! eig_vecs is the eigenvector and lambda is the eigenvalue
       !
       ! note that because LAPACK will overwrite the matricies and they are
       ! expensive to recompute I make copies first
       h_sub_copy = h_sub
       g_sub_copy = g_sub
       subspace_size = num_vec
       previous_eig_vecs = eig_vecs
       call dsygv_wrapper(.false., .true., subspace_size, h_sub_copy, maxit, g_sub_copy, maxit, eig_vals)
       
       if (acpf_root_follow) call root_following(eig_vals, h_sub_copy, g_sub, subspace_size, roots)
       do r = 1, roots 
          energy_change(r) = abs( energy(r) - eig_vals(r) )
          energy(r) = eig_vals(r)
          eig_vecs(1:num_vec,r) = h_sub_copy(1:num_vec,r)
       end do
       
       if (need_to_restart_subspace(num_vec, roots, total_csfs))then
           write(*,*) " "
           write(*,*) "  RESTARTING THE DAVIDSON SUBSPACE "
           write(*,*) " " 
           restarts(i) = .true.
           call restart_davidson_subspace(roots, converged, num_vec, eig_vecs, previous_eig_vecs, c%v, Hc%v, h_sub, g_sub)
       end if

       write(*,9030) i
       time = get_clock_wall_time(iteration_time)

       do iRoot = 1, roots

          ! skip if we already converged this root
          if (converged(iRoot)) cycle

          ! where on disk I need to look for the vector of this iteration(i) associated with this root (iRoot)
          ! this is also the (absolute) number for this vector
          this_vec = root2index(iRoot)

          ! build the new guess for the eigenvector by writing the eigenvector eig_vecs 
          ! in the space of all configurations
          call build_wavefunction(eig_vecs(:,iRoot), c%v, subspace_size, scratch=Hc%v)

          ! run post-iteration checks
          call post_iteration_driver(c%v, iRoot) 

          ! build the residual         
          call build_residual_gen_eig(work, energy(iRoot), eig_vecs(:,iRoot), subspace_size, gmultc, &
               Hc_scratch=Hc%v, Gc_scratch=Gc)
          norm = sqrt(dot_product(work,work))
          ! produce some nice output 
          write(iooutput,9010)
          write(iooutput,9020) i,iRoot,energy(iRoot),norm,energy_change(iRoot),time
          write(*,*) " "
          call write_post_iteration_analysis
          call flush(iooutput)

          ! check convergence
          halt = is_converged(energy_change(iRoot), norm, norm_tol, energy_tol, i, restarts, maxit)
          if (halt) then
             converged(iRoot) = .true.
             write(*,*) " "
             write(*,*) "CONVERGED ROOT = " , iRoot
             write(*,*) " "
             cycle
          end if

          ! if we get here we aren't converged ... build the next vector by preconditioning the residual 
          call build_next_vector_gen_eig(work, energy(iRoot), ref_energy, gmultc, diag_scratch=Hc%v, &
               gc_scratch=Gc)

          ! Orthogonalize against previous vectors then normalize it
          ! Note that I need to orthogonalize against all the vectors, including the ones built for the other
          ! roots this iteration
          call orthonormalize_vector(work, num_vec, scratch=Hc%v)

          ! write out the new expansion vector to disk
          call davidson_io_write(work, cfile, num_vec+1)
          num_vec = num_vec+1
          root2index(iRoot) = num_vec
       end do

       ! check convergence over all roots
       halt = .true.
       do iRoot = 1, roots
          if ( .not. converged(iRoot) ) then
             halt =.false.
          end if
       end do
       if (halt) exit

    end do ! Davidson iterations

    ! produce final output for each root
    energy(1:roots) = energy(1:roots) + ref_energy
    do iRoot = 1, roots
       call build_wavefunction(eig_vecs(:,iRoot),  c%v, num_vec, scratch=Hc%v)
       call post_iteration_driver(c%v, iRoot)
       write(*,*) "   Final wavefunction output for root = ", iRoot
       call write_post_iteration_analysis(.true.)
    end do


    ! memory clean up
    deallocate(h_sub, h_sub_copy, g_sub, g_sub_copy, eig_vals, eig_vecs, stat=ierr)
    call deallocatecheck(ierr, "subspaces in expt_acpf_comp_kernel")
    deallocate(root2index, converged, energy_change, stat=ierr)
    call deallocatecheck(ierr, "extras in gen_eig_fake_multi_state_kernel")

9010 format(/,"   -------  --------   ------------------    ---------------    -------------   ----------",/,&
         "    Iter      Root           Energy             |Hc-EGc|           DeltaE        wall time",/,&
         "   -------  --------   ------------------    ---------------    -------------   ----------")
9020 format('  ',i5,5x,i5,5x,d19.12,3x,d15.8,3x,d15.8,1x,f10.2)
9030 format(/,'  -----------------------------',/'   Iteration ', i3, ' has finished ',/,'  -----------------------------')     
  end subroutine davidson_liu_gen_eig



  !> \brief Expands a Krylov subspace of operator A by a new vector x
  !> \param subspace  The Krylov subspace
  !> \param Ax        The produce of the operator A and the new vector x
  !> \param index     The index of vector x (indicies MUST be contiguous)
  !> \param scratch   The new vector x being added to the Krylov subspace. It is overwritten in this routine
  subroutine expand_krylov_subspace(subspace, Ax, index, scratch)
    ! inputs
    integer, intent(in) :: index
    real(kind=real8), dimension(:) :: Ax(:), scratch(:)
    real(kind=real8), dimension(:,:) :: subspace

    ! local variables 
    integer :: iVec 
    real(kind=real8) :: elem

    do iVec = 1 , index
       call davidson_io_read(scratch, cfile, iVec)
       elem = dot_product( scratch, Ax ) 
       subspace(iVec,index) = elem
       subspace(index,iVec) = elem
    end do
  end subroutine expand_krylov_subspace

  !> \brief Builds a residual \f$ r = \left( H - E * I  \right) \psi \f$ for a standard eigenvalue problem
  !> \param residual On exit contains the residual
  !> \param energy   The eigenvalue associated with this eigenvector (coeff)
  !> \param coeff    The eigenvector \psi written in the krylov subspace
  !> \param n        The number of vectors in the subspace  
  !> \param Hc_scratch       Scratch space for storing H * c products
  !> \param c_scratch        Scratch space for storing the expansion vectors
  subroutine build_residual_std_eig(residual, energy, coeff, n, Hc_scratch, c_scratch)
    ! inputs
    integer, intent(in) :: n
    real(kind=real8), dimension(:) :: Hc_scratch, c_scratch,residual, coeff
    real(kind=real8) :: energy

    ! local variables
    integer :: i 
    !real(kind=real8) :: bound

    ! residual = H*c - E*c
    residual = 0.0
    do i = 1, n 
       call davidson_io_read(Hc_scratch, Hcfile, i)
       call davidson_io_read(c_scratch, cfile, i)
       residual = residual + coeff(i) * ( Hc_scratch - energy * c_scratch)
    end do

    ! Calculate a bound on the eigenvalue using the Residual Norm Bound
    ! (see Zhou, Y.; Shepard, R.; Minkoff, M. Computer Physics Communications 2005, 167, 90–102. )
    ! the final wavefunction is still in civec and we'll use sigmavec as scratch space
    !bound = sqrt( dot_product(residual,residual))
    !write(*,*) "** Residual Norm Eignevalue bound = +/- " , bound
  end subroutine build_residual_std_eig

  !> \brief Builds a residual \f$ r = \left( H - E * I  \right) \psi \f$ for a standard eigenvalue problem
  !> \param residual On exit contains the residual
  !> \param energy   The eigenvalue associated with this eigenvector (coeff)
  !> \param coeff    The eigenvector \psi written in the krylov subspace
  !> \param n        The number of vectors in the subspace  
  !> \param Hc_scratch       Scratch space for storing [H-E_ref] * c products
  !> \param c_scratch        Scratch space for storing the G *c products
  subroutine build_residual_gen_eig(residual, energy, coeff, n, gmultc, Hc_scratch, Gc_scratch)
    ! inputs
    integer, intent(in) :: n
    real(kind=real8), dimension(:) :: Hc_scratch, residual, coeff, gc_scratch
    real(kind=real8) :: energy
    procedure(gtimesc), pointer :: gmultc

    ! local variables 
    integer :: i 

    ! residual = (H-E_ref)*c - E*G*c
    residual = 0.0
    do i = 1, n
       call davidson_io_read(Hc_scratch, Hcfile, i)
       call davidson_io_read(Gc_scratch, cfile, i)
       call gmultc(Gc_scratch)
       residual = residual + coeff(i) * ( Hc_scratch - energy * Gc_scratch)
    end do
  end subroutine build_residual_gen_eig


  !> \brief Checks the convergence criteria for all my Davidson-type solvers, produces a nice output if we are finished
  !> \param energy_change  The change in energy during the last iteration
  !> \param norm           The norm of the residual in the last iteration
  !> \param norm_tol       The convergence threshold for the norm of the residual
  !> \param energy_tol     The convergence threshold for the energy change
  !> \param i              The current iteration
  !> \param restarts       An array recording when we have restarted
  !> \param maxit          The maximum allowed number of iterations
  logical function is_converged(energy_change, norm, norm_tol, energy_tol, i , restarts, maxit)
    implicit none
    integer, intent(in) :: i, maxit
    logical, intent(in) :: restarts(:)
    real(kind=real8), intent(in) :: energy_change, norm, norm_tol, energy_tol

    is_converged = .false.
    
    ! Never converge just after a restart (energy change if often quite low
    ! the iteration after a restart but the residual is still large)
    if ( restarts(i)) then 
        return 
    else if ( i>1 ) then 
        if (restarts(i-1))then
          return 
        end if
    end if
    
    if ( norm < norm_tol ) then 
       write(*,*) "*************************************"
       write(*,*) "  norm of the residual has converged "
       write(*,*) "*************************************"
       is_converged = .true.
    else if ( abs(energy_change) < energy_tol ) then
       ! edge case .. the first iteration always has a 0 energy change
       if ( i > 1 ) then
          write(*,*) "******************************"
          write(*,*) "  energy change has converged "
          write(*,*) "******************************"
          is_converged = .true.
       end if
    else if ( i == maxit ) then
       write(*,*) "*************************************************"
       write(*,*) "diagonalization routine has failed to converge ! "
       write(*,*) "        maxit = " ,maxit, " exceeded    "
       write(*,*) "*************************************************"
       is_converged = .true.
    end if
  end function is_converged


  !> \brief Builds a new vector to expand the krylov subspace using Davidson's formula for a standard eigenvalue problem
  !> \param residual  The residual for the current state
  !> \param energy    The energy of this state
  !> \param diag_scratch      Scratch space for  diagonal elements of the Hamiltonian
  subroutine build_next_vector_std_eig(residual, energy, diag_scratch)
    ! inputs
    real(kind=real8), dimension(:) :: diag_scratch
    real(kind=real8), dimension(:), intent(inout) :: residual
    real(kind=real8) :: energy

    ! local variables
    real(kind=real8) :: precond
    integer :: i

    ! Get the diagonal integrals
    rewind(iodiagstore)
    read(iodiagstore) diag_scratch

    do i = 1, total_csfs
       precond = diag_scratch(i) - energy 
       if ( precond == 0.0 ) then
          residual(i) = 0.0
       else
          residual(i) = -1.0 * residual(i) / precond
       end if
    end do
  end subroutine build_next_vector_std_eig

  !> \brief Builds a new vector to expand the krylov subspace using Davidson's formula for a generalized eigenvalue problem
  !> \param residual  The residual for the current state
  !> \param energy    The energy of this state
  !> \param ref_energy The reference energy
  !> \param gmultc      A pointer to the routine that performs G*c
  !> \param diag_scratch      Scratch space for diagonal elements of the Hamiltonian
  !> \param gc_scratch        Scratch space for a G*c product
  subroutine build_next_vector_gen_eig(residual, energy, ref_energy, gmultc, diag_scratch, gc_scratch)
    ! inputs 
    real(kind=real8), dimension(:) :: diag_scratch, gc_scratch
    real(kind=real8), dimension(:), intent(inout) :: residual
    real(kind=real8) :: energy, ref_energy
    procedure(gtimesc), pointer :: gmultc

    ! local variables 
    real(kind=real8) :: precond
    integer :: i 

    ! Get the diagonal integrals
    rewind(iodiagstore)
    read(iodiagstore) diag_scratch

    ! Get the g values
    gc_scratch = 1.0
    call gmultc(gc_scratch)

    do i = 1, total_csfs
       precond = diag_scratch(i) - ref_energy - energy * gc_scratch(i)
       if ( precond == 0.0 ) then
          residual(i) = 0.0
       else
          residual(i) = -1.0 * residual(i) / precond
       end if
    end do
  end subroutine build_next_vector_gen_eig


  !> \brief Builds the wavefunction by a linear combination of the Krylov vectors
  !> \param coeff   The wavefunction expressed as coefficients of the Krylov vectors
  !> \param wavefunc On exit contains the wavefunction
  !> \param scratch  A scratch vector
  !> \param N        The number of vectors in the Krylov subspace
  subroutine build_wavefunction(coeff, wavefunc, N, scratch)
    ! inputs
    real(kind=real8), dimension(:) :: coeff, wavefunc, scratch 
    integer, intent(in) :: N

    ! local variables
    integer :: i 
    real(kind=real8) :: norm

    wavefunc = 0.0 
    do i =1, N
       call davidson_io_read(scratch, cfile, i)
       wavefunc = wavefunc + coeff(i) * scratch
    end do
    ! formally we shouldn't have to renormalize. Still doesn't hurt to be careful
    ! when using floating points
    norm = dot_product(wavefunc, wavefunc)
    norm = sqrt(norm)
    wavefunc = (1/norm) * wavefunc 
  end subroutine build_wavefunction

  !> \brief Orthogonalizes a new vector in our krylov subpsace against the previous vectors on disk
  !> \param vector  The new vector, on exit it has orthonormalized
  !> \param scratch Scratch space
  !> \param n       The number of vectors in the subspace  
  !>
  !> We perform the orthnormalization with Modified Gram Schmidt (MGS)
  subroutine orthonormalize_vector(vector, n, scratch)
    implicit none
    ! inputs
    integer, intent(in) :: n 
    real(kind=real8), dimension(:) :: vector, scratch

    ! local variables
    real(kind=real8) :: overlap, norm
    integer :: i 

    do i = 1, n 
       call davidson_io_read(scratch, cfile, i)
       overlap = dot_product( vector, scratch)
       vector = vector - overlap * scratch
       norm = sqrt( dot_product(vector,vector) )
       if (abs(norm) < 1E-14) then
          vector = 0.0
       else
          vector = vector / norm
       end if
    end do
  end subroutine orthonormalize_vector

  !> \brief Root following, aka "The Columbus Trick". Reorders the eigenvalues/eigenvectors on 
  !> the krylov subspace to better match the state(s) we are solving for
  !> \param evals The eigenvalues
  !> \param evecs The eigenvectors
  !> \param metric   The metric the eigenvectors are normalized under ( <i|metric|j> = \delta_ij )
  !> \param n     The current dimension of the krylov subspace
  !> \param nstates The number of states we are solving for
  subroutine root_following(evals, evecs, metric, n, nstates)
    implicit none

    ! inputs
    real(kind=real8), intent(in) :: metric(:,:)
    real(kind=real8), intent(inout) :: evals(:), evecs(:,:)
    integer, intent(in) :: n, nstates

    ! local variables
    integer :: i, istate, ierr, j, lda
    real(kind=real8) :: overlap, max_overlap
    logical, allocatable :: used(:)
    real(kind=real8), allocatable :: copy_evals(:), tmp(:), best_guess(:)
    real(kind=real8), allocatable :: copy_evecs(:,:)
    character(len=4) :: number_string
    character(len=20) :: vec_format

    ! external functions
    real(kind=real8), external :: ddot 

    ! header output
200 format(/," Best fit for state ", I5, " found in krylov subspace eigenvector " , I5)
201 format(" Root Following " /, " -------------- ")
    write(*,201)
100 format ("ACPF root following help enabled")
101 format ("  (vectors written in the basis of the krylov subspace)")
102 format ("  eigvector ", I4, "     =")
    if (ACPF_ROOT_FOLLOW_HELP) then
       write(number_string, '(I0)') n
       vec_format = '(' // number_string // 'F7.3)'
       write(*,100)      
       write(*,101)
       do i = 1, n 
          write(*,102, advance='no') i 
          write(*, vec_format) evecs(1:n,i)
       end do
       write(*,*)
    end if

    allocate(tmp(n), copy_evals(n), copy_evecs(n,n), used(n), best_guess(n), stat=ierr)
    call allocatecheck(ierr, "temporaries in root_following")
    used = .false.

    ! For each state, find the eigenvector which best matches
    do istate = 1, nstates
       max_overlap = 0.0
       j = 0
       do i = 1, n
          ! our best guess for what the state should look like comes for our 
          ! initial guess for this state, which is the ith kryov vector for the ith state
          best_guess(:) = 0.0
          best_guess(istate) = 1.0

          ! calculate the overlap as <ref | metric | eigen_state>
          lda = size(metric,1)
          call dgemv( 'n', n, n, 1.0, metric, lda, evecs(1,i), 1, 0.0, tmp, 1)  
          overlap =  abs(ddot(n, tmp, 1, best_guess, 1))  

          if (overlap> max_overlap) then 
             j = i 
             max_overlap = overlap
          end if
       end do

       ! at this point the jth eigenvector overlaps the best with the ith state
       if (used(j)) then
          write(*,*) " "
          write(*,*) " Weird root following error"
          write(*,*) " One eigenvector has the highest overlap with multiple target states"
          write(*,*) " Be careful ...."
          write(*,*) " "
       end if
       write(*,200) istate, j
       copy_evecs(1:n,istate) = evecs(1:n,j)
       copy_evals(istate) = evals(j)
       used(j) = .true.
    end do

    ! Copy the data back into the evecs, evals
    do istate = 1, nstates
       evecs(1:n, istate) = copy_evecs(1:n, istate)
       evals(istate) = copy_evals(istate)
    end do

    deallocate(tmp, best_guess, used, copy_evecs, copy_evals)
  end subroutine root_following

end module davidson_mod
