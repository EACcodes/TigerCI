! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module twostwo_seg_info_replacement
use IOBuffer
use global_var_mod

type pseudo_twostwo_seg_info
  integer::blockSize
  integer::fileUnitInts
  integer::fileUnitReals
  integer::maxRec
end type pseudo_twostwo_seg_info

contains

subroutine setup_twostwo_seg(file,maxmem,unitnoInts,unitnoReals,blockSize,numberOfThreads,poolIDInts,poolIDReals)
  implicit none
  type(pseudo_twostwo_seg_info),intent(inout)::file
  integer,intent(in)::maxmem,unitnoInts,unitnoReals,blockSize,numberOfThreads
  integer,intent(out)::poolIDInts,poolIDReals

#ifdef DEBUG_TWOTWO_TIGER
  write(*,*) "DEBUG: About to create two buffers for the twostwo seg data."
  flush(6)
#endif
  
  ! setup buffer
  call for_int_buf_construct(maxmem/(blockSize*8),blockSize,0,numberOfThreads,poolIDInts)
  call for_int_buf_openfile(poolIDInts,unitnoInts,scratch_directory//'2s2_seg_info.dat',len(scratch_directory)+16)
  ! here, we just want one tenth of the above buffer, seems fair
  call for_double_buf_construct(maxmem/(blockSize*real8*10),blockSize,0,numberOfThreads,poolIDReals)
  call for_double_buf_openfile(poolIDReals,unitnoReals,scratch_directory//'2s2_seg_info_real.dat',len(scratch_directory)+21)
  
  file%blockSize = blockSize
  file%fileUnitInts = unitnoInts
  file%fileUnitReals = unitnoReals
  file%maxRec = 0
#ifdef DEBUG_TWOTWO_TIGER
  write(*,*) "DEBUG: Created two buffers for the twostwo seg data."
  flush(6)
#endif

end subroutine setup_twostwo_seg

subroutine twostwo_seg_writeData(pos,file,threadID,ij_2,lambda_weight,mu_weight,ibar,jbar,internal_singles,&
         Gij_internal,loop_type,start_elec,i,j)
  implicit none

  type(pseudo_twostwo_seg_info),intent(inout)::file
  integer,intent(in)::pos,threadID
  integer,intent(in)::ij_2,lambda_weight,mu_weight,ibar,jbar,internal_singles,&
         loop_type,start_elec,i,j
  real(real8),intent(in)::Gij_internal
  integer::tmp
  
  tmp = 10*(pos-1)+1
  
#ifdef DEBUG_TWOTWO_TIGER
  write(*,*) "DEBUG: TWOSTWO writing at position ",pos,tmp
  write(*,*) "DEBUG: TWOSTWO Data to be written ",ij_2,lambda_weight,mu_weight,ibar,jbar,internal_singles,&
         Gij_internal,loop_type,start_elec,i,j
 flush(6)
#endif
  
  call for_int_buf_writeElement(file%fileUnitInts,tmp+1 ,ij_2            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+2 ,lambda_weight   ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+3 ,mu_weight       ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+4 ,ibar            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+5 ,jbar            ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+6 ,internal_singles,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+7 ,loop_type       ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+8 ,start_elec      ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+9 ,i               ,threadID)
  call for_int_buf_writeElement(file%fileUnitInts,tmp+10,j               ,threadID)
  
  call for_double_buf_writeElement(file%fileUnitReals,pos,Gij_internal   ,threadID)
  
  file%maxRec = max(file%maxRec,tmp)
  
end subroutine twostwo_seg_writeData

subroutine twostwo_seg_readData(pos,file,threadID,ij_2,lambda_weight,mu_weight,ibar,jbar,internal_singles,&
         Gij_internal,loop_type,start_elec,i,j,flag)
  implicit none

  type(pseudo_twostwo_seg_info),intent(in)::file
  integer,intent(in)::pos,threadID
  integer,intent(out)::ij_2,lambda_weight,mu_weight,ibar,jbar,internal_singles,&
         loop_type,start_elec,i,j,flag
  real(real8)::Gij_internal
  integer::tmp
  
  tmp = 10*(pos-1)+1
  
  if(file%maxRec < tmp) then
    flag = 1
    return
  else
    flag = 0
  endif
  
  call for_int_buf_readElement(file%fileUnitInts,tmp+1 ,ij_2            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+2 ,lambda_weight   ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+3 ,mu_weight       ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+4 ,ibar            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+5 ,jbar            ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+6 ,internal_singles,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+7 ,loop_type       ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+8 ,start_elec      ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+9 ,i               ,threadID)
  call for_int_buf_readElement(file%fileUnitInts,tmp+10,j               ,threadID)
  
  call for_double_buf_readElement(file%fileUnitReals,pos,Gij_internal   ,threadID)

#ifdef DEBUG_TWOTWO_TIGER
  write(*,*) "DEBUG: TWOSTWO reading at position ",pos,tmp
  write(*,*) "DEBUG: TWOSTWO Data read ",ij_2,lambda_weight,mu_weight,ibar,jbar,internal_singles,&
         Gij_internal,loop_type,start_elec,i,j
 flush(6)
#endif
  
end subroutine twostwo_seg_readData

end module twostwo_seg_info_replacement
