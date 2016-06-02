! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module arraylist_integer

! this is a VERY BASIC ArrayList implementation for Fortran.
! resizing strategy by doubling size mimicks the JRE's behavior.
! other than that: A TEMPORA! O MORES! O FORTRAN!
! author: Johannes Dieterich (who gets grey over missing features in this "language")

private

public :: construct,append,addElement,changeElement,resetIterator,getElement,getNext,destroy,arrayContent_integer

type arrayContent_integer
   integer,dimension(:),allocatable,private::storage
   integer,private::highestUsed
   integer,private::iterCounter
end type arrayContent_integer

contains

subroutine construct(content,startSize)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  integer,intent(in)::startSize
  
  integer::status
  
  allocate(content%storage(startSize),stat=status)
  if(status /= 0) then
    write(*,*) "Couldn't allocate storage provider for integer ArrayList. Stopping. ",status
    flush(6)
    stop
  endif
  
  content%highestUsed = 0
  
end subroutine construct

function append(content,element)
  implicit none
  type(arrayContent_integer),intent(inout)::content
  integer,intent(in)::element
  integer::append
  
  integer::currSize,position
  
  position = content%highestUsed+1
  currSize = size(content%storage,1)
  if(currSize < position) then
    call resize(content,position)
  endif
  
  content%storage(position) = element
  content%highestUsed = position
  append = position
  
end function append

subroutine resetIterator(content)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  
  content%iterCounter = 1
  
end subroutine resetIterator

function getNext(content)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  integer::getNext
  
#ifdef DEBUG_TIGER
  if(content%iterCounter > size(content%storage,1)) then
    write(*,*) "DEBUG: Trying to get element beyond end of list. STOPPING!"
    flush(6)
    stop
  endif
#endif
  
  getNext = content%storage(content%iterCounter)
  content%iterCounter = content%iterCounter + 1
  
end function getNext

subroutine addElement(content,element,position)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  integer,intent(in)::element
  integer,intent(in)::position
  
  integer::currSize
  
#ifdef DEBUG_TIGER
  if(content%highestUsed+1 < position) then
     write(*,*) "DEBUG: Adding element in non-sequential fashion at position ",position
     write(*,*) "DEBUG: Previous was ",content%highestUsed
     write(*,*) "DEBUG: Are you REALLY sure you wanted to do this?"
  endif
#endif
  
  currSize = size(content%storage,1)
  if(currSize < position) then
    call resize(content,position)
  endif
  
  content%storage(position) = element
  content%highestUsed = max(content%highestUsed,position)
  
end subroutine addElement

function changeElement(content,newElement,position)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  integer,intent(in)::newElement,position
  integer::changeElement
  
  ! we assume that the person calling this knows what she is doing and the list is long enough
  changeElement = content%storage(position)
  content%storage(position) = newElement
  
end function changeElement

function getElement(content,position)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  integer,intent(in)::position
  integer::getElement
  
  ! again we assume the user to know what she is doing
  getElement = content%storage(position)
  
end function getElement

subroutine destroy(content)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  
  integer::status
  
  deallocate(content%storage,stat=status)
  if(status /= 0) then
    write(*,*) "Deallocation of storage provider of integer ArrayList failed. ",status
  endif
  
end subroutine destroy

recursive subroutine resize(content,minNewSize)
  implicit none
  
  type(arrayContent_integer),intent(inout)::content
  integer,intent(in)::minNewSize
  
  integer,dimension(:),allocatable::backup
  integer::status,currSize,newSize
  
  ! always double the size
  currSize = size(content%storage,1)
  newSize = 2*currSize
  
  ! allocate backup
  allocate(backup(currSize),stat=status)
  if(status /= 0) then
    write(*,*) "Couldn't allocate backup storage provider for integer ArrayList. Stopping. ",status
    flush(6)
    stop
  endif
  
  ! copy
  backup(:)=content%storage(:)
  
  ! deallocate and reallocate
  deallocate(content%storage,stat=status)
  if(status /= 0) then
    write(*,*) "Couldn't deallocate storage provider for integer ArrayList. Stopping. ",status
    flush(6)
    stop
  endif
  allocate(content%storage(newSize),stat=status)
  if(status /= 0) then
    write(*,*) "Couldn't reallocate storage provider for integer ArrayList. Stopping. ",status
    flush(6)
    stop
  endif
  
  ! copy back over
  content%storage(1:currSize) = backup(:)
  
  ! deallocate backup
  deallocate(backup,stat=status)
  if(status /= 0) then
    write(*,*) "Couldn't deallocate backup storage provider for integer ArrayList. Stopping. ",status
    flush(6)
    stop
  endif
  
  ! check if we are big enough yet. If not: recurse
  if(newSize < minNewSize) then
     call resize(content,minNewSize)
  endif
  
end subroutine resize

end module arraylist_integer
