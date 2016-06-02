! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module two_seg_info_replacement
 use IOBuffer
 use global_var_mod

type pseudo_two_seg_info
  integer::blockSize
  integer::fileUnitInts
  integer::maxRec
end type pseudo_two_seg_info

contains

subroutine setup_two_seg(file,maxmem,unitnoInts,blockSize,numberOfThreads,poolIDInts)
  implicit none
  type(pseudo_two_seg_info),intent(inout)::file
  integer,intent(in)::maxmem,unitnoInts,blockSize,numberOfThreads
  integer,intent(out)::poolIDInts

#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: About to create one buffer for the two seg data."
  flush(6)
#endif
  
  ! setup buffer
  call for_int_buf_construct(maxmem/(blockSize*8),blockSize,0,numberOfThreads,poolIDInts)
  call for_int_buf_openfile(poolIDInts,unitnoInts,scratch_directory // 'two_seg_info.dat',len(scratch_directory) + 16)
  
  file%blockSize = blockSize
  file%fileUnitInts = unitnoInts
  file%maxRec = 0
  
#ifdef DEBUG_TIGER
  write(*,*) "DEBUG: Created one buffer for the two seg data."
  flush(6)
#endif

end subroutine setup_two_seg

subroutine two_seg_writeData(pos,file,threadID,ij_2,s1,s2,s3,s4,weight,rt_loop_weight,num_singles,loop_type,start_elec)
  implicit none

  type(pseudo_two_seg_info)::file
  integer,intent(in)::pos,threadID
  integer,intent(in)::ij_2,s1,s2,s3,s4,weight,rt_loop_weight,num_singles,loop_type,start_elec
  integer::tmp

  tmp = 10*(pos-1)+1
    
  call for_int_buf_writeElement(file%fileUnitInts,tmp+1 ,ij_2          ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+2 ,s1            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+3 ,s2            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+4 ,s3            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+5 ,s4            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+6 ,weight        ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+7 ,rt_loop_weight,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+8 ,num_singles   ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+9 ,loop_type     ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+10,start_elec    ,threadID)
  
  file%maxRec = max(file%maxRec,tmp)
  
end subroutine two_seg_writeData

subroutine two_seg_readData(pos,file,threadID,ij_2,s1,s2,s3,s4,weight,rt_loop_weight,num_singles,loop_type,start_elec,flag)
  implicit none

  type(pseudo_two_seg_info)::file
  integer,intent(in)::pos,threadID
  integer,intent(out)::ij_2,s1,s2,s3,s4,weight,rt_loop_weight,num_singles,loop_type,start_elec,flag
  integer::tmp

  tmp = 10*(pos-1)+1
  
  if(file%maxRec < tmp) then
    flag = 1
    return
  else
    flag = 0
  endif
  
  call for_int_buf_readElement(file%fileUnitInts,tmp+1 ,ij_2          ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+2 ,s1            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+3 ,s2            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+4 ,s3            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+5 ,s4            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+6 ,weight        ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+7 ,rt_loop_weight,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+8 ,num_singles   ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+9 ,loop_type     ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+10,start_elec    ,threadID)

end subroutine two_seg_readData

end module two_seg_info_replacement
