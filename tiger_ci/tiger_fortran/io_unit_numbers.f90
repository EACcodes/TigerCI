! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!*******************************************************************
  !> \brief This is a collection of unit numbers used within the
  !> TigerCI code. Please use the identifiers here and extend when new
  !> files are needed (and remove old ones).
  !>
  !> \date 2012
  !> \author Johannes M. Dieterich
!*******************************************************************
module io_unit_numbers

   implicit none
   
  integer, parameter::ioOutput    = 6
  integer, parameter::ioError     = 6 
  integer, parameter::ioInput     = 5
  
  integer, parameter::ioii        = 53
  integer, parameter::ioiijj      = 55
  integer, parameter::ioIntegral  = 50
  integer, parameter::ioijij      = 72

  integer, parameter ::io_output_buff = 73
  
  integer, parameter::ioDiag      = 74
  integer, parameter::ioB         = 75
  integer, parameter::ioAB        = 76
  integer, parameter::ioCIFinal   = 77
  integer, parameter::ioAltCsf    = 78
  integer, parameter::ioDiagStore = 80

  integer, parameter::iodiag1     = 96
  integer, parameter::iojdtemp    = 124
  integer, parameter::iojdtemp2   = 126
  integer, parameter::iojdtemp3   = 127
  integer, parameter::ioRef       = 128
   
   ! files for the 2e-4i-integrals, all of them currently buffered
   ! real numbered w/ semi-uniform access
   integer,parameter::cd_iiab_no = 320
   integer,parameter::cd_iabc_no = 340
   integer,parameter::cd_abcd_no = 345
   integer,parameter::cd_iajb_no = 351
   integer,parameter::cd_iajb_unsorted_no = 251
   integer,parameter::cd_ijab_no = 352
   integer,parameter::cd_iajk_no = 366
   integer,parameter::cd_ikaj_no = 369
   integer,parameter::cd_ijka_no = 371
   
   ! files w/ integer records (random access)
   integer,parameter::iajb_ind_no = 353
   integer,parameter::iajk_bl_no = 368
   integer,parameter::ikaj_bl_no = 370
   integer,parameter::ijka_bl_no = 372
   
   ! files with no specified record size (sequential access)
   integer,parameter::cd_iaib_no = 322
   integer,parameter::cd_ijpp_no = 330
   integer,parameter::cd_ipjp_no = 331
   integer,parameter::cd_iapp_no = 335
   integer,parameter::cd_ipap_no = 336
   integer,parameter::rediaia_no = 338
   integer,parameter::rediaib_no = 339
   
   integer,parameter::iajb_bl_no = 354
   integer,parameter::cd_ijja_no = 360
   integer,parameter::cd_ijia_no = 361
   integer,parameter::cd_ijkl_no = 365

   integer,parameter::cd_test_no = 307
   integer,parameter::ijka_test_no = 308

   integer,parameter::cd_restart_no = 362
   integer,parameter::mo_int_no = 305
   integer,parameter::cvec_no = 302
   integer,parameter::cvec_index_no = 304
   
   integer,parameter::twosingles_no=404
   integer,parameter::twostwoseg_ints=603
   integer,parameter::twostwoseg_reals=604
   integer,parameter::abcd_seg_no=605
   integer,parameter::iabc_seg_no=606
   integer,parameter::two_seg_no=607
   
   integer,parameter::ioci_info   = 96
   integer,parameter::ioints =  244
   
   integer,parameter::ioH = 405
   
#ifdef TIGER_USE_OMP
   integer,parameter::threeiseg_no= 609
#endif
   
end module io_unit_numbers
