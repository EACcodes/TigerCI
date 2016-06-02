! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module blocked_locks_mod
  use global_var_mod
  use utilities_mod
#ifdef TIGER_USE_OMP
  use omp_lib
#endif

  type blockedLockType
#ifdef TIGER_USE_OMP
      integer(kind=OMP_lock_kind),dimension(:),allocatable::blocked_locks
     integer::blockSize
     integer::blocks
#endif 
  end type blockedLockType
       

  type blockedLockVectorType
     real(kind=real8), dimension(:), allocatable :: v
     type(blockedLockType) :: l
  end type blockedLockVectorType
  
  
  contains
  
  subroutine allocateBlockedLockVector(foo, length, blockSize)
      implicit none
      ! inputs
      type(blockedLockVectorType),intent(inout) :: foo
      integer,intent(in) :: length, blockSize
      ! local variables
      integer :: ierr
      allocate(foo%v(length), stat=ierr)
      call allocatecheck(ierr, "foo%vector")
#ifdef TIGER_USE_OMP
      call setupLocks(foo%l, length, blockSize)
#else
      ! nothing to be done w/ blockSize
      if(blockSize < 0) then
        write(*,*) "Unused blocksize is negative! ",blockSize
      endif
#endif
  end subroutine allocateBlockedLockVector
  
  subroutine deallocateBlockedLockVector(foo)
      type(blockedLockVectorType),intent(inout) :: foo
      integer :: ierr
      deallocate(foo%v,stat=ierr)
      call deallocatecheck(ierr, "foo%v")
#ifdef TIGER_USE_OMP
      call deallocateLocks(foo%l)
#endif
  end subroutine
      

#ifdef TIGER_USE_OMP
  subroutine setupLocks(locks, length, blockSize)
      implicit none
      ! input
      type(blockedLockType),intent(inout) :: locks
      integer, intent(in) :: length, blockSize
      ! local variables
      integer :: ierr, i
      
      locks%blockSize = blockSize
      locks%blocks    = (length-1)/blockSize + 1
      allocate(locks%blocked_locks(locks%blocks), stat=ierr)
      call allocatecheck(ierr, "locks%blocked_locks")
      do i = 1, locks%blocks
          call OMP_init_lock(locks%blocked_locks(i))
      end do
  end subroutine setupLocks
  
  subroutine deallocateLocks(locks)
      type(blockedLockType),intent(inout):: locks
      integer :: i 
      do i = 1, locks%blocks
          call OMP_destroy_lock(locks%blocked_locks(i))
      end do
  end subroutine
  
  function setLock(locks, point)
      implicit none
      ! input
      type(blockedLockType),intent(inout) :: locks
      integer,intent(in) :: point
      ! output
      integer :: setLock
      setLock = (point-1)/locks%blockSize+1
      call OMP_set_lock(locks%blocked_locks(setLock))
  end function
  
  subroutine unsetLock(locks, point)
      implicit none
      ! input
      type(blockedLockType),intent(inout) :: locks
      integer,intent(in) :: point
      call OMP_unset_lock(locks%blocked_locks(point))
  end subroutine unsetLock
  
  subroutine setLocks(locks, startp, endp)
      implicit none
      ! input
      type(blockedLockType),intent(inout) :: locks
      integer,intent(in) :: startp, endp
      ! local variables
      integer :: startbl, endbl, i
      
      startbl=(startp-1)/locks%blockSize+1
      endbl=(endp-1)/locks%blockSize+1
  
  do i=startbl,endbl
    call OMP_set_lock(locks%blocked_locks(i))
  enddo
  end subroutine setLocks
  
  subroutine unsetLocks(locks, startp, endp)
      implicit none
      ! input
      type(blockedLockType),intent(inout) :: locks
      integer,intent(in) :: startp, endp
      ! local variables
      integer :: startbl, endbl, i
      
        startbl=(startp-1)/locks%blockSize+1
  endbl=(endp-1)/locks%blockSize+1
  
  do i=startbl,endbl
    call OMP_unset_lock(locks%blocked_locks(i))
  enddo
  end subroutine unsetLocks
#endif
  
end module
