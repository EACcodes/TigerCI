! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module stack_mod
    use utilities_mod
    implicit none
    
    ! default starting size for a stack
    integer, parameter, private :: START = 1
    
    ! an integer stack
    type IntegerStackType
        integer :: head
        integer, dimension(:), allocatable :: s
    end type
    
    ! a stack of integer stacks
    type IntegerStackStackType      
        integer :: head
        type(IntegerStackType), dimension(:), allocatable :: s
    end type
    
    ! interfaces allowing client code to use stacks of multiple types
    ! using the name for each basic operation
    interface push_stack
        module procedure push_stack_integer_stack
        module procedure push_stack_integer_stack_stack
    end interface
        
    interface pop_stack
        module procedure pop_stack_integer_stack
        module procedure pop_stack_integer_stack_stack
    end interface
    
    interface is_stack_empty
        module procedure is_stack_empty_integer_stack
        module procedure is_stack_empty_integer_stack_stack
    end interface

    interface copy_stack
        module procedure copy_stack_integer_stack
        !module procedure copy_stack_integer_stack_stack
    end interface

    interface print_stack
        module procedure print_stack_integer_stack
        module procedure print_stack_integer_stack_stack
    end interface

    interface empty_stack
       module procedure empty_stack_integer_stack
    end interface empty_stack

    interface get_stack_size
       module procedure get_stack_size_integer_stack
    end interface

    contains
    
      integer function get_stack_size_integer_stack(stack)
        type(IntegerStackType) ,intent(in) :: stack 
        get_stack_size_integer_stack = stack%head
      end function  get_stack_size_integer_stack
        

    subroutine push_stack_integer_stack(stack, val)
        ! inputs
        type(IntegerStackType) :: stack 
        integer :: val 
        
        ! local variables
        integer, dimension(:), allocatable :: tmp
        integer :: old_size, new_size
        integer :: ierr

        ! handle the stack memory
        if (.not. allocated(stack%s) ) then 
            allocate(stack%s(START), stat=ierr)
            call allocatecheck(ierr, "start of integer stack")
            stack%head = 0
        else if ( size(stack%s) < stack%head+1 ) then
            ! we need to allocate a larger stack 
            old_size = size(stack%s)
            new_size = old_size * 2
            if (new_size == 0) new_size = START           ! ugly, I think I have cases where the memory gets move_alloc'd but I still use the stack variable
            allocate(tmp(new_size), stat=ierr)
            call allocatecheck(ierr, "allocate for integer stack expansion")
            tmp(1:old_size) = stack%s(1:old_size)
            call move_alloc(tmp, stack%s)
        end if
        
        ! add the new value to the stack
        if ( size(stack%s) == 0 ) call cause_floating_point_exception
        stack%head = stack%head + 1
        stack%s(stack%head) = val 
    end subroutine push_stack_integer_stack
    
    integer function pop_stack_integer_stack(stack)
        ! inputs
        type(IntegerStackType) :: stack 
        if ( (.not. allocated(stack%s)) .or. (stack%head == 0) ) then
            pop_stack_integer_stack = 0 
        else
            pop_stack_integer_stack = stack%s(stack%head)
            stack%head = stack%head - 1
        end if
    end function pop_stack_integer_stack
    
    logical function is_stack_empty_integer_stack(stack)
        ! inputs
        type(IntegerStackType) :: stack
        if ( (.not. allocated(stack%s)) .or. (stack%head == 0) ) then
            is_stack_empty_integer_stack = .true.
        else
            is_stack_empty_integer_stack =.false.
        end if
    end function is_stack_empty_integer_stack
    
    subroutine copy_stack_integer_stack(src, dest)
        ! inputs
        type(IntegerStackType) :: src, dest
        
        ! local variables
        integer, dimension(:), allocatable :: tmp
        integer :: ierr

        ! error checking
        if (.not. allocated(src%s) ) then 
           write(*,*) "FATAL ERROR copy_stack_integer_stack"
           write(*,*) "Tried to copy an unallocated stack"
           stop
        end if
        
        ! set the memory first
        if ( (.not. allocated(dest%s)) .or. (dest%head < src%head)) then
            allocate(tmp(src%head), stat=ierr)
            call allocatecheck(ierr, "allocating a stack during a stack copy")

            tmp(1:src%head) = src%s(1:src%head)
            if ( allocated(dest%s) ) deallocate(dest%s)
            call move_alloc(tmp, dest%s)
        end if
        
        ! copy values
        dest%head = src%head
        dest%s(1:src%head) = src%s(1:src%head)
    end subroutine copy_stack_integer_stack
    
    subroutine print_stack_integer_stack(stack)
        ! inputs
        type(IntegerStackType) :: stack
        
        ! local variables
        integer :: i
        if (.not. allocated(stack%s) ) return 
        do i = stack%head, 1, -1
            write(*,*) stack%s(i)
        end do
    end subroutine print_stack_integer_stack

    subroutine empty_stack_integer_stack(stack)
      ! inputs
      type(IntegerStackType) :: stack 
      stack%head = 0 
    end subroutine empty_stack_integer_stack
    
    subroutine test_integer_stack()
        type(IntegerStackType) :: stack, stack2
        integer :: i
        integer, parameter :: n = 100
        do i = 1, n
           call push_stack(stack,i)
        end do
        call copy_stack(stack, stack2)
        write(*,*) "Testing the stack"
        do i = 1, n
           write(*,*) pop_stack(stack) 
        end do
        write(*,*) "Testing the stack copy"
        do i = 1, n 
            write(*,*) pop_stack(stack2)
        end do
  end subroutine test_integer_stack
    
  subroutine push_stack_integer_stack_stack(stack, item)
      ! inputs
      type(IntegerStackStackType) :: stack 
      type(IntegerStackType) :: item 
      
      ! local variables
      type(IntegerStackType), dimension(:), allocatable :: tmp
      integer :: old_size, new_size
      integer :: ierr, i
      
      ! handle the memory for the stack 
      if (.not. allocated(stack%s) ) then 
          allocate(stack%s(START), stat=ierr)
          call allocatecheck(ierr, "stack of a stack of integer stacks")
          stack%head = 0 
      else if ( size(stack%s) == stack%head ) then 
          ! we need to allocate a larger stack
          ! note that this is going to require moving a fairly large amount of memory
          ! this might need to be optimized later ...
          old_size = size(stack%s) 
          new_size = old_size * 2
          allocate(tmp(new_size), stat=ierr)
          call allocatecheck(ierr, "allocate for a stack of integer stacks expansion")
          do i = 1, old_size
              call move_alloc(stack%s(i)%s, tmp(i)%s)
              tmp(i)%head = stack%s(i)%head 
          end do
          call move_alloc(tmp, stack%s)
      end if 
          
      ! add the new value to the stack
      stack%head = stack%head + 1
      call copy_stack(item, stack%s(stack%head))
  end subroutine push_stack_integer_stack_stack
  
  type(IntegerStackType) function pop_stack_integer_stack_stack(stack)
    ! inputs
    type(IntegerStackStackType) :: stack 
    if ( (.not. allocated(stack%s)) .or. (stack%head == 0)) then 
        write(*,*) "Error you tried to pop from an empty stack"
        write(*,*) "and I haven't implemented a good solution for this yet"
        write(*,*) "for a stack of integer stacks ..."
        stop
    end if
    call copy_stack(stack%s(stack%head), pop_stack_integer_stack_stack)
    stack%head = stack%head - 1
  end function
  
  logical function is_stack_empty_integer_stack_stack(stack)
    ! inputs
    type(IntegerStackStackType) :: stack 
    if ( (.not. allocated(stack%s)) .or. (stack%head==0))then 
        is_stack_empty_integer_stack_stack = .true.
    else
        is_stack_empty_integer_stack_stack = .false.
    end if
  end function 
  
!  subroutine copy_stack_integer_stack_stack(stack)
!      ! inputs
!      type(IntegerStackStackType) :: stack 
!      write(*,*) "The copy_stack_integer_stack_stack is not implemented"
!      stop
!  end subroutine
  
  subroutine print_stack_integer_stack_stack(stack)
      ! inputs
      type(IntegerStackStackType) :: stack 
      
      ! local variables
      integer :: i 
      
      write(*,*) "Printing out a stack of integer stacks"
      if ( .not. allocated(stack%s) ) return 
      do i = stack%head, 1, -1 
          write(*,*) "Stack #", i
          call print_stack(stack%s(i))
      end do
  end subroutine print_stack_integer_stack_stack
  
  subroutine test_integer_stack_stack()
      type(IntegerStackStackType) :: stack_of_stacks
      type(IntegerStackType) :: stack1, stack2, stack3
      integer, parameter :: n = 10
      integer :: i 
      
      ! set up the integer stacks
      do i = 1, n 
          call push_stack(stack1, i)
          call push_stack(stack2, i*2)
      end do
      
      ! push the integer stacks onto their stack 
      call push_stack(stack_of_stacks, stack1)
      call push_stack(stack_of_stacks, stack2)
      
      ! test the stack print 
      write(*,*) "Testing the stack print"
      call print_stack(stack_of_stacks)
      
      ! test a pop
      write(*,*) "Testing the pop"
      do i = 1, 2 
          stack3 = pop_stack(stack_of_stacks)
          call print_stack(stack3)
      end do
  end subroutine

end module stack_mod
