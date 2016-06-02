! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module abcd_seg_info_replacement
 use IOBuffer
 use global_var_mod


type pseudo_abcd_seg_info
  integer::blockSize
  integer::fileUnitInts
  integer::maxRec
end type pseudo_abcd_seg_info

contains

subroutine setup_abcd_seg(file,maxmem,unitnoInts,blockSize,numberOfThreads,poolIDInts)
  implicit none
  type(pseudo_abcd_seg_info),intent(inout)::file
  integer,intent(in)::maxmem,unitnoInts,blockSize,numberOfThreads
  integer,intent(out)::poolIDInts

#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: About to create one buffer for the abcd seg data."
  flush(6)
#endif
  
  ! setup buffer
  call for_int_buf_construct(maxmem/(blockSize*8),blockSize,0,numberOfThreads,poolIDInts)
  call for_int_buf_openfile(poolIDInts,unitnoInts,scratch_directory // 'abcd_seg_info.dat',len(scratch_directory)+17)
  
  file%blockSize = blockSize
  file%fileUnitInts = unitnoInts
  file%maxRec = 0
#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: Created one buffer for the abcd seg data."
  flush(6)
#endif

end subroutine setup_abcd_seg

subroutine abcd_seg_writeData(pos,file,threadID,ab_2,head_weight,single_spin,single_spinm2,total_spin,a_check2,b_check2)
  implicit none

  type(pseudo_abcd_seg_info)::file
  integer,intent(in)::pos,threadID
  integer,intent(in)::ab_2,head_weight,single_spin,single_spinm2,total_spin,a_check2,b_check2
  integer::tmp
  
  tmp = 7*(pos-1)+1
    
  call for_int_buf_writeElement(file%fileUnitInts,tmp+1 ,ab_2         ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+2 ,head_weight  ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+3 ,single_spin  ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+4 ,single_spinm2,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+5 ,total_spin   ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+6 ,a_check2     ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+7 ,b_check2     ,threadID)
  
  file%maxRec = max(file%maxRec,tmp)

end subroutine abcd_seg_writeData

subroutine abcd_seg_readData(pos,file,threadID,ab_2,head_weight,single_spin,single_spinm2,total_spin,a_check2,b_check2,flag)
  implicit none

  type(pseudo_abcd_seg_info)::file
  integer,intent(in)::pos,threadID
  integer,intent(out)::ab_2,head_weight,single_spin,single_spinm2,total_spin,a_check2,b_check2,flag
  integer::tmp
  
  tmp = 7*(pos-1)+1
  
  if(file%maxRec < tmp) then
    flag = 1
    return
  else
    flag = 0
  endif
  
  call for_int_buf_readElement(file%fileUnitInts,tmp+1 ,ab_2         ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+2 ,head_weight  ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+3 ,single_spin  ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+4 ,single_spinm2,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+5 ,total_spin   ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+6 ,a_check2     ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+7 ,b_check2     ,threadID)
  
end subroutine abcd_seg_readData

end module abcd_seg_info_replacement
