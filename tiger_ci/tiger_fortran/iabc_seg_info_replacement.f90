! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module iabc_seg_info_replacement
 use IOBuffer
 use global_var_mod

type pseudo_iabc_seg_info
  integer::blockSize
  integer::fileUnitInts
  integer::maxRec
end type pseudo_iabc_seg_info

contains

subroutine setup_iabc_seg(file,maxmem,unitnoInts,blockSize,numberOfThreads,poolIDInts)
  implicit none
  type(pseudo_iabc_seg_info),intent(inout)::file
  integer,intent(in)::maxmem,unitnoInts,blockSize,numberOfThreads
  integer,intent(out)::poolIDInts

#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: About to create one buffer for the iabc seg data."
  flush(6)
#endif
  
  ! setup buffer
  call for_int_buf_construct(maxmem/(blockSize*8),blockSize,0,numberOfThreads,poolIDInts)
  call for_int_buf_openfile(poolIDInts,unitnoInts,scratch_directory // 'iabc_seg_info.dat',len(scratch_directory)  + 17)
  
  file%blockSize = blockSize
  file%fileUnitInts = unitnoInts
  file%maxRec = 0
  
#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: Created one buffer for the iabc seg data."
  flush(6)
#endif

end subroutine setup_iabc_seg

subroutine iabc_seg_writeData(pos,file,threadID,i_ind,lambda_weight,mu_weight,internal_singles,ibar,loop_type)
  implicit none

  type(pseudo_iabc_seg_info)::file
  integer,intent(in)::pos,threadID
  integer,intent(in)::i_ind,lambda_weight,mu_weight,internal_singles,ibar,loop_type
  integer::tmp
  
  tmp = 6*(pos-1)+1
    
  call for_int_buf_writeElement(file%fileUnitInts,tmp+1 ,i_ind           ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+2 ,lambda_weight   ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+3 ,mu_weight       ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+4 ,internal_singles,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+5 ,ibar            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+6 ,loop_type       ,threadID)
  
  file%maxRec = max(file%maxRec,tmp)

end subroutine iabc_seg_writeData

subroutine iabc_seg_readData(pos,file,threadID,i_ind,lambda_weight,mu_weight,internal_singles,ibar,loop_type,flag)
  implicit none

  type(pseudo_iabc_seg_info)::file
  integer,intent(in)::pos,threadID
  integer,intent(out)::i_ind,lambda_weight,mu_weight,internal_singles,ibar,loop_type,flag
  integer::tmp
  
  tmp = 6*(pos-1)+1
  
  if(file%maxRec < tmp) then
    flag = 1
    return
  else
    flag = 0
  endif
  
  call for_int_buf_readElement(file%fileUnitInts,tmp+1 ,i_ind           ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+2 ,lambda_weight   ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+3 ,mu_weight       ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+4 ,internal_singles,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+5 ,ibar            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+6 ,loop_type       ,threadID)
  
end subroutine iabc_seg_readData

end module iabc_seg_info_replacement
