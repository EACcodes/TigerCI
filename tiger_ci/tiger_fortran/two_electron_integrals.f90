! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
! Here we handle the two electron integrals using their analytic formulas.
! We make no use of the cholesky decomposition. Note that this code is 
! STRICTLY for debugging. It is horribly inefficient. 

#ifndef TIGER_GAMESS
module two_electron_integrals
  use cholesky_structs
  use io_unit_numbers
  use global_var_mod
  use molecule_var_mod
  use utilities_mod
  use get_integrals 
  use get_basis_data
  use molecular_orbital_mod
  use IOBuffer
#ifdef TIGER_USE_OMP
  use omp_lib
#endif

  implicit none

  real(kind=real8), dimension(:,:,:,:), allocatable :: ao_ints, mo_ints

contains

  !> \brief Checks the accuracy of our MO cholesky vectors
  subroutine check_MO_CD(cho_data)
    implicit none

    integer :: i,j,k,l
    integer :: ij, kl, file 
    integer :: thread, numthreads, ierr
    real(kind=real8) :: cholesky_value
    real(kind=real8), allocatable :: vec1(:,:), vec2(:,:)
    type(cholesky_data), intent(in) :: cho_data

    if (.not. allocated(mo_ints)) call bootstrap(.true.)
#ifdef TIGER_USE_OMP
    numthreads = numberOfThreads
#else
    numthreads = 1
#endif
    allocate(vec1(numcho, numthreads), vec2(numcho, numthreads), stat=ierr)
    call allocatecheck(ierr, "cholesky vectors in check_MO_CD")

    file = mo_int_no        !gfortran OMP implementation doesn't like module variables in the omp pragma?
    !$omp parallel &
    !$omp default(none) &
    !$omp shared(vec1, vec2, cho_data, file, num_orbitals) &
    !$omp private(i,j,k,l,ij,kl,cholesky_value, thread)
    thread = 1
#ifdef TIGER_USE_OMP
    thread = OMP_get_thread_num() + 1
#endif
    !$omp do schedule(static)
    do i = 1, num_orbitals
       do j = 1, i
          ij = i*(i-1)/2 + j
          ij = cho_data%mo_ind_inv(ij)
          if (ij == 0) continue
          call for_double_buf_readblock(file, ij, vec1(:,thread), thread)

          do k = 1, num_orbitals
             do l = 1, k 
                kl = k*(k-1)/2 + j 
                kl = cho_data%mo_ind_inv(kl)
                if (kl == 0) continue
                call for_double_buf_readblock(file, kl, vec2(:,thread), thread)

                cholesky_value = dot_product(vec1(:,thread), vec2(:,thread))
                call check_integral(cholesky_value, i, j, k, l)  
             end do
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(vec1, vec2)
  end subroutine check_MO_CD

  !> \brief Checks the accuracy of our AO cholesky vectors
  subroutine check_AO_CD(AOCho, P, pairMap, activePairs)
    implicit none

    real(kind=real8), intent(in) :: AOCho(:,:)
    integer, intent(in) :: P(:), activePairs
    type(pairIndex), intent(in) :: pairMap

    integer :: i,j,k,l
    integer :: ij, kl, unpivoted_ij, unpivoted_kl
    integer :: thread
    real(kind=real8) :: cholesky_value


    if (.not. allocated(ao_ints)) call bootstrap(.false.)

    write(*,*) "Checking validity of the AO cholesky decomposition"
    write(*,*) "**************************************************"

    !$omp parallel &
    !$omp default(none) &
    !$omp shared(activePairs, P, pairMap, AOCho) &
    !$omp private(cholesky_value, i, j, k, l, ij, kl, unpivoted_ij, unpivoted_kl, thread)
#ifdef TIGER_USE_OMP
    thread = OMP_get_thread_num() + 1  
#else
    thread = 1
#endif
    !$omp do schedule(static)
    do ij = 1, activePairs
       unpivoted_ij = P(ij)
       i = pairMap%pair2ij(unpivoted_ij)%i
       j = pairMap%pair2ij(unpivoted_ij)%j

       do kl = 1, ij
          unpivoted_kl = P(kl)
          k = pairMap%pair2ij(unpivoted_kl)%i
          l = pairMap%pair2ij(unpivoted_kl)%j

          cholesky_value = dot_product(AOCho(ij,:), AOCho(kl,:))
          call check_integral_ao(cholesky_value,i,j,k,l)
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine check_AO_CD


  !> \brief This runs if we try to check an integral value, but we haven't computed the integral tables ao_ints or mo_ints yet
  subroutine bootstrap(mo_flag)
    implicit none
    logical, intent(in) :: mo_flag

    ! AO integrals
    if (.not. allocated(ao_ints)) call get_two_electron_integrals()

    ! Transform from AO to MO 
    if (mo_flag .and. (.not. allocated(mo_ints))) call ao2mo_quarter(molecular_orbitals)
  end subroutine bootstrap

  subroutine get_two_electron_integrals
    implicit none

    ! Variables
    !**********
    integer :: ierr
    integer :: iShell, jShell, kShell, lShell
    integer :: ijShell, klShell
    integer :: iFunc, jFunc, kFunc, lFunc
    integer :: nShell, nTInt
    integer :: biggest
    integer :: ip, jp, kp, lp
    integer :: i,j,k,l,ijkl
    
    integer :: p_i, p_ij, p_ijk

    integer, dimension(number_bas) :: numbfshl

    real(kind=real8), dimension(:), allocatable :: TInt

    ! callback routine
    external :: integral2_ucla

    ! Memory allocations
    !********************

    ! Build the ugly 2e ints matrix
    allocate(ao_ints(number_bas,number_bas,number_bas,number_bas),stat=ierr)
    call allocatecheck(ierr,"ao_ints  ")
    ao_ints = 0.0

    numbfshl = 0
    nShell = 0 
    ! get the basis set info
    call basis_func_data(numbfshl, number_bas, nShell)

    ! allocate the TInt array to hold the integrals for a given set of shells based on the
    ! maximum number of functions in all shells
    biggest = maxval(numbfshl)
    nTInt =  biggest*biggest*biggest*biggest
    allocate(TInt(nTInt), stat = ierr)
    call allocatecheck(ierr,"TInt      ")
    TInt= 0.0

    ! Integral evaluations 
    !***********************      
    ! loop over shells for the eval_ijkl call at each iteration 
    ! we evaluate all the integrals for ishell,jshell,kshell,lshell
    do iShell = 1, nShell
       do jShell= 1, iShell
          do kShell = 1, nShell
             do lShell= 1, kShell
                ijShell = iShell*(iShell+1)/2 + jShell
                klShell = kShell*(kShell+1)/2 + lShell

                ! The way I'm doing loops I occasionally hit something unphysical ... whoops
                if ( jShell > nShell .or. lShell > nShell) cycle

                ! By permutational symmetry I only need ijShell <= klShell (or the opposite .. just not both)
                if ( ijShell <= klShell ) then

                   ! evaluate integrals
                   ! lets try the default Tiger callback routine 
                   TInt = 0.0
                   call eval_ijkl(iShell,jShell,kShell,lShell,TInt,nTInt)

                   ! the following are the number of basis functions before this shell
                   ! you can consider them conceptually to be "pointers" to the current position
                   ! in the basis set
                   ip = sum( numbfshl(1:iShell-1))
                   jp = sum( numbfshl(1:jShell-1))
                   kp = sum( numbfshl(1:kShell-1))
                   lp = sum( numbfshl(1:lShell-1))

                   ! pointer to the TInt array for the (ij|kl) integral
                   ijkl = 1

                   ! loop over all basis functions in the shells we evaluated
                   do iFunc = 1, numbfshl(iShell)
                      p_i = (iFunc - 1) * numbfshl(jShell)
                      do jFunc = 1, numbfshl(jShell)
                         p_ij = (p_i + jFunc - 1) * numbfshl(kShell)
                         do kFunc = 1, numbfshl(kShell)
                            p_ijk = (p_ij + kFunc - 1) * numbfshl(lShell)
                            do lFunc = 1, numbfshl(lShell)

                               ijkl = p_ijk + lFunc - 1
                               ! These are the absolute basis function numbers for this integral
                               i = ip + iFunc
                               j = jp + jFunc
                               k = kp + kFunc
                               l = lp + lfunc

                               ! In case I ever want to look into the ordering again ...
                               !100 format("(",i5, 1X,i5,"|",i5,1X,i5,") = ",f10.8)
                               !write(*,100) i, j, k, l, TInt(ijkl) 

                               ! Record integral (ij|kl)
                               ao_ints(i,j,k,l) = TInt(ijkl + 1)

                               ! Normally you would just want to record that one integral
                               ! and ignore all the equivalent ones via permutational symmetry
                               ! but for debugging purposes (making the code as readable as possible)
                               ! we are going to record all the equivalent ones explicitly
                               !(kl|ij)
                               ao_ints(k,l,i,j) = TInt(ijkl + 1)
                               !(ji|kl)
                               ao_ints(j,i,k,l) = TInt(ijkl + 1)
                               !(ij|lk)
                               ao_ints(i,j,l,k) = TInt(ijkl + 1)
                               !(ji|lk)
                               ao_ints(j,i,l,k) = TInt(ijkl + 1)
                               !(kl|ji)
                               ao_ints(k,l,j,i) = TInt(ijkl + 1)
                               !(lk|ij)
                               ao_ints(l,k,i,j) = TInt(ijkl + 1)
                               !(lk|ji)
                               ao_ints(l,k,j,i) = TInt(ijkl + 1)
                            end do !lFunc
                         end do !kFunc
                      end do ! jFunc 
                   end do ! iFunc 
                end if ! ij <= kl

             end do
          end do
       end do
    end do

    ! Final debugging check
    !        123 format("(",i5, 1X,i5,"|",i5,1X,i5,") = ",f10.8)
    !        do iFunc = 1, number_bas
    !            do jFunc = 1, number_bas
    !                do kFunc = 1, number_bas
    !                    do lFunc = 1, number_bas
    !                        write(*,123) iFunc,jFunc,kFunc,lFunc,ao_ints(iFunc,jFunc,kFunc,lFunc)
    !                    end do
    !                end do
    !            end do
    !        end do
    !        call flush(6)


  end subroutine get_two_electron_integrals

  subroutine ao2mo(ao_ints, mo_coeffs )

    !**************************************************************
    !//  AO2MO Transforms integrals in AO basis to MO basis
    !//  Does a direct transformation with an inefficient ON^8 scaling
    !// 
    !//  WRITTEN BY VICTOR OYEYEMI, 2012
    !**************************************************************


    implicit none

    !MO indices
    integer::i,j,k,l

    !AO indices
    integer::a,b,c,d
    real::C_ai,C_bj,C_ck,C_dl
    real (kind=8) :: ao_ints(:,:,:,:)
    real (kind=8) :: mo_coeffs(:,:)

    i = 0
    j = 0
    k = 0
    l = 0
    a = 0
    b = 0
    c = 0
    d = 0


    allocate(mo_ints(num_orbitals,num_orbitals,num_orbitals,num_orbitals))
    mo_ints=0

    do i = 1,num_orbitals
       do j = 1,num_orbitals
          do k = 1,num_orbitals
             do l = 1,num_orbitals

                !Comput sum
                do a = 1,num_orbitals
                   C_ai=mo_coeffs(a,i)
                   do b = 1,num_orbitals
                      C_bj=mo_coeffs(b,j)
                      do c = 1,num_orbitals
                         C_ck=mo_coeffs(c,k)
                         do d = 1,num_orbitals
                            C_dl=mo_coeffs(d,l)
                            mo_ints(l,k,j,i) = mo_ints(l,k,j,i) + C_ai*C_bj*C_ck*C_dl*ao_ints(d,c,b,a)
                         enddo
                      enddo
                   enddo
                enddo

             enddo
          enddo
       enddo
    enddo

  end subroutine ao2mo

  subroutine ao2mo_quarter( mo_coeffs )

    !**************************************************************
    !//  AO2MO Transforms integrals in AO basis to MO basis
    !//  Does quarter ransformations with ON^5 scaling
    !// 
    !//  WRITTEN BY VICTOR OYEYEMI, 2012
    !**************************************************************


    implicit none

    !MO indices
    integer::i,j,k,l

    !AO indices
    integer::a,b,c,d
    real::C_ai,C_bj,C_ck,C_dl
    real (kind=8) :: mo_coeffs(:,:)
    real (kind=8), allocatable :: T1(:,:,:,:), T2(:,:,:,:),T3(:,:,:,:) !temporary quarter transforms

    i = 0
    j = 0
    k = 0
    l = 0
    a = 0
    b = 0
    c = 0
    d = 0


    allocate(mo_ints(num_orbitals,num_orbitals,num_orbitals,num_orbitals))
    allocate(T1(num_orbitals,num_orbitals,num_orbitals,num_orbitals))
    allocate(T2(num_orbitals,num_orbitals,num_orbitals,num_orbitals))
    allocate(T3(num_orbitals,num_orbitals,num_orbitals,num_orbitals))
    mo_ints=0
    T1=0
    T2=0
    T3=0

    write(*,*) "Beginning first quarter transform"
    call flush(6)
    do i = 1,num_orbitals
       do b = 1,num_orbitals
          do c = 1,num_orbitals
             do d = 1,num_orbitals

                !Compute sum
                do a = 1,num_orbitals
                   C_ai=mo_coeffs(a,i)
                   T1(d,c,b,i) = T1(d,c,b,i) + C_ai*ao_ints(d,c,b,a)
                enddo

             enddo
          enddo
       enddo
    enddo
    write(*,*) "Finished first quarter transform"
    write(*,*) "Beginning second quarter transform"
    call flush(6)
    do i = 1,num_orbitals
       do j = 1,num_orbitals
          do c = 1,num_orbitals
             do d = 1,num_orbitals

                !Compute sum
                do b = 1,num_orbitals
                   C_bj=mo_coeffs(b,j)
                   T2(d,c,j,i) = T2(d,c,j,i) + C_bj*T1(d,c,b,i)
                enddo

             enddo
          enddo
       enddo
    enddo
    write(*,*) "Finished second quarter transform"
    write(*,*) "Beginning third quarter transform"
    call flush(6)
    do i = 1,num_orbitals
       do j = 1,num_orbitals
          do k = 1,num_orbitals
             do d = 1,num_orbitals

                !Compute sum
                do c = 1,num_orbitals
                   C_ck=mo_coeffs(c,k)
                   T3(d,k,j,i) = T3(d,k,j,i) + C_ck*T2(d,c,j,i)
                enddo

             enddo
          enddo
       enddo
    enddo
    write(*,*) "Finished third quarter transform"
    write(*,*) "Beginning fourth quarter transform"
    call flush(6)
    do i = 1,num_orbitals
       do j = 1,num_orbitals
          do k = 1,num_orbitals
             do l = 1,num_orbitals

                !Compute sum
                do d = 1,num_orbitals
                   C_dl=mo_coeffs(d,l)
                   mo_ints(l,k,j,i) = mo_ints(l,k,j,i) + C_dl*T3(d,k,j,i)
                enddo

             enddo
          enddo
       enddo
    enddo
    write(*,*) "Finished fourth quarter transform"
    call flush(6)
    deallocate(T1,T2,T3)

! Final debugging check
 !           123 format("(",i5, 1X,i5,"|",i5,1X,i5,") = ",f10.8)
            do i = 1, number_bas
                do j = 1, number_bas
                    do k = 1, number_bas
                        do l = 1, number_bas
!                            write(*,123) i,j,k,l,mo_ints(i,j,k,l)
                        end do
                    end do
                end do
            end do
            call flush(6)


  end subroutine ao2mo_quarter


  subroutine check_integral ( value, i, j, k, l )
    implicit none

    ! inputs
    real(kind=real8) :: value        ! Value of the integral to be checked
    integer :: i,j,k,l               ! MO integral (ij|kl)

    ! variables
    real(kind=real8) :: diff

    if ( .not. allocated(mo_ints)) then
       call bootstrap(.true.)
    end if

    diff = abs( value - mo_ints(i,j,k,l))
201 format ('ERROR ON INTEGRAL (',I5,1X,I5,'|',I5,1X,I5,'), VALUE = ', E12.5, ' ANALYTIC = ', E12.5, ' DIFFERENCE = ', E12.5)
203 format ('INTEGRAL OK (',I5,1X,I5,'|',I5,1X,I5,')')
    if ( diff > 1.0D-7 ) then 
       write(*,201) i,j,k,l,value, mo_ints(i,j,k,l),diff
    else 
       write(*,203) i,j,k,l
    end if
  end subroutine check_integral

  subroutine check_integral_ao ( value, i, j, k, l )
    implicit none

    ! inputs
    real(kind=real8) :: value        ! Value of the integral to be checked
    integer :: i,j,k,l               ! AO integral (ij|kl)

    ! variables
    real(kind=real8) :: diff

    if ( .not. allocated(ao_ints)) then
       call bootstrap(.false.)
    end if

    diff = abs( value - ao_ints(i,j,k,l))
201 format ('ERROR ON INTEGRAL (',I5,1X,I5,'|',I5,1X,I5,'), VALUE = ', E12.5, ' ANALYTIC = ', E12.5, ' DIFFERENCE = ', E12.5)
203 format ('INTEGRAL OK (',I5,1X,I5,'|',I5,1X,I5,')')
    if ( diff > 1.0D-7 ) then 
       write(*,201) i,j,k,l,value, ao_ints(i,j,k,l), diff
    else
       write(*,203) i,j,k,l
    end if
  end subroutine check_integral_ao


end module two_electron_integrals
#endif
