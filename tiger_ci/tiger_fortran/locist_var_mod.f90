! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!> Just a bunch of definitions of data and scratch arrays to be used in the rest of the code.
!> Johannes M. Dieterich
module locist_var_mod

  use global_var_mod
  use utilities_mod
  
  type locist_scratch
  
      ! integers
      integer::num_virt
      integer::number_paths
  
      ! vectors
      integer,dimension(:),allocatable::virt_lam_allow
      integer,dimension(:),allocatable::virt_mu_allow
      real(real8),dimension(:),allocatable::sigma_nm1_mu
  
      ! matrices
      real(real8),dimension(:,:),allocatable::integrals_abcd
      real(real8),dimension(:,:),allocatable::sigma_nm2_lambda
      real(real8),dimension(:,:),allocatable::sigma_nm2_mu
  end type locist_scratch
  
  contains
  
  subroutine allocLocScratch(loc_scr,num_external)
  
  implicit none
  type(locist_scratch)::loc_scr
  integer::num_external,allocatestatus
  
  allocate(loc_scr%integrals_abcd(num_external,num_external), &
         stat = allocatestatus)
  if(allocatestatus.ne.0) then
     write(*,*) "ERROR: Failure to allocate int_abcd. ",allocatestatus
     flush(6)
     stop
  endif
  
  end subroutine allocLocScratch
  
  subroutine deallocLocScratch(loc_scr)
  
  implicit none
  type(locist_scratch)::loc_scr
  
  call try_deallocate_real_2D(loc_scr%integrals_abcd, "integrals_abcd")
  call try_deallocate_real_2D(loc_scr%sigma_nm2_lambda, "sigma_nm2_lambda")
  call try_deallocate_real_2D(loc_scr%sigma_nm2_mu, "sigma_nm2_mu")
  call try_deallocate_real_1D(loc_scr%sigma_nm1_mu, "sigma_nm1_mu")
  call try_deallocate_int_1D(loc_scr%virt_lam_allow, "virt_lam_allow")
  call try_deallocate_int_1D(loc_scr%virt_mu_allow, "virt_mu_allow")
  
  loc_scr%num_virt = 0 
  loc_scr%number_paths = 0
  end subroutine deallocLocScratch
  
  subroutine copyLocScratch(this, that)
  
  type(locist_scratch)::this,that
  
  ! simple first: the integers
  that%number_paths = this%number_paths
  that%num_virt = this%num_virt
  
  if(allocated(this%virt_mu_allow)) then
     call copyArrayInt(this%virt_mu_allow,that%virt_mu_allow)
  endif
  if(allocated(this%virt_lam_allow)) then
     call copyArrayInt(this%virt_lam_allow,that%virt_lam_allow)
  endif
  if(allocated(this%sigma_nm1_mu)) then
     call copyArray(this%sigma_nm1_mu,that%sigma_nm1_mu)
  endif
  if(allocated(this%integrals_abcd)) then
     call copyMatrix(this%integrals_abcd,that%integrals_abcd)
  endif
  if(allocated(this%sigma_nm2_lambda)) then
     call copyMatrix(this%sigma_nm2_lambda,that%sigma_nm2_lambda)
  endif
  if(allocated(this%sigma_nm2_mu)) then
     call copyMatrix(this%sigma_nm2_mu,that%sigma_nm2_mu)
  endif
  
  end subroutine copyLocScratch

end module locist_var_mod
