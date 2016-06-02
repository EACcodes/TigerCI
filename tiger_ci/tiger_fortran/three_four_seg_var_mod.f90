! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module three_four_seg_var_mod

  use global_var_mod
  
  implicit none
  
#ifdef TIGER_USE_OMP
  !// THESE ARRAYS ARE *ONLY* USED TO COMMUNICATE INFORMATION (OFFSETS ETC)
  !// IN THE OPENMP CASE WITHIN THIS MODULE. DO NOT USE FOR ANYTHING ELSE!
  !// ALSO, DO NOT PUT ANY OTHER GLOBAL MODULE VARIABLES HERE WITHOUT 
  !// PROPER JUSTIFICATION!
  integer,dimension(:),allocatable,public:: omp_offsets_abcd_icount,omp_offsets_iabc_icount
#endif
  
  type threefourmod1vars
    integer::iabc_c
    integer::a_ic
    integer::b_ic
    integer::i_ind
    integer::abcd_rec
    integer::iabc_rec
    integer::iter_count
    
    integer,dimension(:),allocatable::the_virtuals
  
    real(real8),dimension(:,:),allocatable::scep_ci_1
    real(real8),dimension(:,:),allocatable::scep_sigma_1
    real(real8),dimension(:,:),allocatable::Px
    real(real8),dimension(:,:),allocatable::integral_buf  !FOR (IA|IB) AND (IA|IA) INTEGRALS
  end type threefourmod1vars

  type threefourmod2vars
    integer::iter_count
    integer::ijab_rec,ijaj_rec
    integer::ij_ind

    logical:: ijka_flag
  
    !// SCEPPER MATRICES AND COUPLING COEFFICIENTS
    real(real8),dimension(:,:),allocatable::scep_integral_C
    real(real8),dimension(:,:),allocatable::scep_integral_EX
    real(real8),dimension(:,:),allocatable::Px
    real(real8),dimension(:,:),allocatable::Py            
    real(real8),dimension(:,:),allocatable::scep_ci_1
    real(real8),dimension(:,:),allocatable::scep_ci_2
    real(real8),dimension(:,:),allocatable::scep_sigma_1
    real(real8),dimension(:,:),allocatable::scep_sigma_2 
                                            
    real(real8),dimension(:,:),allocatable::aijk_buf  
                                            
    real(real8),dimension(:,:),allocatable::intgrl_sin
  end type threefourmod2vars
  
  contains
   
  subroutine allocMod2Vars(mod2vars,num_external,fsnopen)
  implicit none
  
  type(threefourmod2vars)::mod2vars
  integer::num_external,allocatestatus,fsnopen
  
  allocatestatus = 0
  
  allocate(mod2vars%scep_ci_1(num_external,num_external),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate scep_ci_1. ",allocatestatus
    stop
  endif
  allocate(mod2vars%scep_ci_2(num_external,num_external),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate scep_ci_2. ",allocatestatus
    stop
  endif
  allocate(mod2vars%scep_sigma_1(num_external,num_external),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate scep_sigma_1. ",allocatestatus
    stop
  endif
  allocate(mod2vars%scep_sigma_2(num_external,num_external),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate scep_sigma_2. ",allocatestatus
    stop
  endif
  allocate(mod2vars%scep_integral_EX(num_external,num_external),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate scep_integral_EX. ",allocatestatus
    stop
  endif
  allocate(mod2vars%scep_integral_C(num_external,num_external),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate scep_integral_C. ",allocatestatus
    stop
  endif
  allocate(mod2vars%intgrl_sin(num_external,num_external),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate intgrl_sin. ",allocatestatus
    stop
  endif
  allocate(mod2vars%Px(fsnopen,fsnopen),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate Px. ",allocatestatus
    stop
  endif
  allocate(mod2vars%Py(fsnopen,fsnopen),stat=allocatestatus)
  if(allocatestatus.ne.0) then
    write(*,*) "ERROR: Failure to allocate Py. ",allocatestatus
    stop
  endif
  
  end subroutine allocMod2vars
  
  subroutine zeroMod2Vars(mod2vars)
  implicit none
  
  type(threefourmod2vars)::mod2vars
  
  mod2vars%scep_ci_1 = 0
  mod2vars%scep_ci_2 = 0
  mod2vars%scep_sigma_1 = 0
  mod2vars%scep_sigma_2 = 0
  mod2vars%scep_integral_ex = 0
  mod2vars%scep_integral_c= 0
  mod2vars%intgrl_sin = 0  
  mod2vars%Px = 0
  mod2vars%Py = 0
  
  end subroutine zeroMod2Vars
  
  subroutine deallocMod1Vars(mod1vars)
  
  implicit none
  
  type(threefourmod1vars)::mod1vars
  integer::dealloccheck
  
  if(allocated(mod1vars%the_virtuals)) then
    deallocate(mod1vars%the_virtuals,stat=dealloccheck)
    if(dealloccheck.ne.0) then
      write(*,*) "ERROR: Deallocating the virtuals in threefourmod1vars. ",dealloccheck
      stop
    endif
  endif
  if(allocated(mod1vars%scep_ci_1)) then
    deallocate(mod1vars%scep_ci_1,stat=dealloccheck)
    if(dealloccheck.ne.0) then
      write(*,*) "ERROR: Deallocating scep_ci_1 in threefourmod1vars. ",dealloccheck
      stop
    endif
  endif
  if(allocated(mod1vars%scep_sigma_1)) then
    deallocate(mod1vars%scep_sigma_1,stat=dealloccheck)
    if(dealloccheck.ne.0) then
      write(*,*) "ERROR: Deallocating scep_sigma_1 in threefourmod1vars. ",dealloccheck
      stop
    endif
  endif
  if(allocated(mod1vars%Px)) then
    deallocate(mod1vars%Px,stat=dealloccheck)
    if(dealloccheck.ne.0) then
      write(*,*) "ERROR: Deallocating Px in threefourmod1vars. ",dealloccheck
      stop
    endif
  endif
  if(allocated(mod1vars%integral_buf)) then
    deallocate(mod1vars%integral_buf,stat=dealloccheck)
    if(dealloccheck.ne.0) then
      write(*,*) "ERROR: Deallocating integral_buf in threefourmod1vars. ",dealloccheck
      stop
    endif
  endif
  end subroutine deallocMod1Vars
  
  subroutine deallocMod2Vars(mod2vars)
  
  implicit none
  
  type(threefourmod2vars)::mod2vars
  integer::dealloccheck
  
  deallocate(mod2vars%scep_ci_1,mod2vars%scep_ci_2,mod2vars%scep_sigma_1,mod2vars%scep_sigma_2,&
             mod2vars%scep_integral_EX,mod2vars%scep_integral_C,mod2vars%intgrl_sin,stat=dealloccheck)
             
  if(dealloccheck.ne.0) then
     write(*,*) "ERROR: Deallocating scepper matrices in threefourmod2vars. ",dealloccheck
     stop
  endif
  
  deallocate(mod2vars%Px,mod2vars%Py,stat=dealloccheck)
  
  if(dealloccheck.ne.0) then
     write(*,*) "ERROR: Deallocating Px/Py matrices in threefourmod2vars. ",dealloccheck
     stop
  endif
  
  end subroutine deallocMod2Vars

end module three_four_seg_var_mod
