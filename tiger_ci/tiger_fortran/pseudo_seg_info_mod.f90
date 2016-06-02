! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module pseudo_seg_info_mod

  use twostwo_seg_info_replacement
  use two_seg_info_replacement
  use abcd_seg_info_replacement
  use iabc_seg_info_replacement
  use global_var_mod
  use io_unit_numbers
  use IOBuffer  

  implicit none

  type(pseudo_twostwo_seg_info)::pseudo_twostwo_info
  type(pseudo_two_seg_info)::pseudo_two_info
  type(pseudo_abcd_seg_info)::pseudo_abcd_info
  type(pseudo_iabc_seg_info)::pseudo_iabc_info

  contains

  subroutine setupPseudoInfoFiles()

    implicit none
    
    ! initialize the two_seg_info_replacement pseudofile
    call setup_two_seg(pseudo_two_info,for_buf_maxmemSegs,two_seg_no,for_buf_blocksizeSegs,numberOfThreads,for_buf_twoSeg_poolID)
#ifdef DEBUG_TIGER
    write(*,*) "DEBUG: Initialized pseudo_two_info file with a maxSize of ",for_buf_maxmemSegs
#endif

    ! initialize the abcd_seg_info_replacement pseudofile
    call setup_abcd_seg(pseudo_abcd_info,for_buf_maxmemSegs,abcd_seg_no,for_buf_blocksizeSegs,numberOfThreads,for_buf_abcdSeg_poolID)
#ifdef DEBUG_TIGER
    write(*,*) "DEBUG: Initialized pseudo_abcd_info file with a maxSize of ",for_buf_maxmemSegs
#endif

    ! initialize the iabc_seg_info_replacement pseudofile
    call setup_iabc_seg(pseudo_iabc_info,for_buf_maxmemSegs,iabc_seg_no,for_buf_blocksizeSegs,numberOfThreads,for_buf_iabcSeg_poolID)
#ifdef DEBUG_TIGER
    write(*,*) "DEBUG: Initialized pseudo_iabc_info file with a maxSize of ",for_buf_maxmemSegs
#endif

    ! setup the twostwo_seg_info_replacement pseudofile
    call setup_twostwo_seg(pseudo_twostwo_info,for_buf_maxmemSegs,twostwoseg_ints,twostwoseg_reals, &
         for_buf_blocksizeSegs,numberOfThreads,for_buf_twostwoSeg_poolID,for_buf_twostwoSegReal_poolID)
#ifdef DEBUG_TIGER
    write(*,*) "DEBUG: Initialized pseudo_twostwo_info file with a maxSize of ",for_buf_maxmemSegs
#endif

  end subroutine setupPseudoInfoFiles
  
  subroutine deletePseudoInfoFile()
  
    implicit none
    
    ! close all our small buffer pools
    call for_int_buf_removepool(for_buf_twostwoSeg_poolID)
    write(ioOutput,*) "DEBUG: Took care of buffer I."
    flush(ioOutput)
    call for_double_buf_removepool(for_buf_twostwoSegReal_poolID)
    write(ioOutput,*) "DEBUG: Took care of buffer II."
    flush(ioOutput)
    call for_int_buf_removepool(for_buf_abcdSeg_poolID)
    write(ioOutput,*) "DEBUG: Took care of buffer III."
    flush(ioOutput)
    call for_int_buf_removepool(for_buf_iabcSeg_poolID)
    write(ioOutput,*) "DEBUG: Took care of buffer IV."
    flush(ioOutput)
    call for_int_buf_removepool(for_buf_twoSeg_poolID)
    write(ioOutput,*) "DEBUG: Took care of buffer V."
    flush(ioOutput)
  
  end subroutine deletePseudoInfoFile

end module pseudo_seg_info_mod
