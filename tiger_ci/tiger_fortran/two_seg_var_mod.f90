! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module two_seg_var_mod

  use global_var_mod
  use utilities_mod

  type twoModVars
    integer::iter_count
    integer::i_ind
    integer::twostwo_rec
    integer::twosingles_rec

    real(real8), dimension(:), allocatable::Gij_matrix  
    real(real8), dimension(:),allocatable::aibi_debug
    real(real8), dimension(:,:), allocatable::Hab

    real(real8), dimension(:,:), allocatable::scep_ci_1,scep_ci_2

  end type twoModVars
  
  contains
  
  subroutine copyTwoModVars(this,that)
  implicit none
  type(twoModVars)::this,that
  
  ! scalar values
  that%iter_count = this%iter_count
  that%i_ind = this%i_ind
  that%twostwo_rec = this%twostwo_rec
  that%twosingles_rec = this%twosingles_rec
  
  ! arrays
  if(allocated(this%Gij_matrix)) then
    call copyArray(this%Gij_matrix,that%Gij_matrix)
  endif
  if(allocated(this%aibi_debug)) then
    call copyArray(this%aibi_debug,that%aibi_debug)
  endif
  
  ! matrices
  if(allocated(this%scep_ci_1)) then
    call copyMatrix(this%scep_ci_1,that%scep_ci_1)
  endif
  if(allocated(this%scep_ci_2)) then
    call copyMatrix(this%scep_ci_2,that%scep_ci_2)
  endif
  
  ! tensor
  if(allocated(this%Hab)) then
    call copyMatrix(this%Hab,that%Hab)
  endif
  
  end subroutine copyTwoModVars

end module two_seg_var_mod
