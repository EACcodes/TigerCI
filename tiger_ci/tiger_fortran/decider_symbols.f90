! Copyright (c) 1998-2004, University of California, Los Angeles, Emily A. Carter
!                       2004-2016, Princeton University, Emily A. Carter
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of
!     conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!     of conditions and the following disclaimer in the documentation and/or other
!     materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!     used to endorse or promote products derived from this software without specific
!     prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
! CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
! OF THE POSSIBILITY OF SUCH DAMAGE.
!> \brief Contains all the symbols used in the decider
!> 
!> Placed here to make it simple to use/add new symbols

module decider_symbols
    implicit none
    
    integer, parameter :: DEBUG                 = 0
    integer, parameter :: BLD_HEAD_ELEM         = 1
    integer, parameter :: UNSUPPORTED           = 3
    integer, parameter :: ADD_CNSTR_LVL         = 4
    integer, parameter :: PATH_WRITER           = 5
    integer, parameter :: FOUR_INT_COMPL_VEC_MO = 6
    integer, parameter :: HOLES_INTERNAL_2      = 7
    integer, parameter :: HOLES_INTERNAL_1      = 8
    integer, parameter :: GTIMESC_NO_EXT        = 9
    integer, parameter :: RECORD_REF            = 10
    integer, parameter :: PURELY_EXT_COMP       = 11
    integer, parameter :: PURELY_INT_COMP       = 12
    integer, parameter :: INT_EXT_COMP          = 13
    integer, parameter :: DIAG_COMP             = 14
    integer, parameter :: LMO_EXT_2_COMP        = 15
    integer, parameter :: LMO_INT_2_COMP        = 16
    integer, parameter :: LMO_1_INT_COMP        = 17
    integer, parameter :: LMO_1_EXT_COMP        = 18
    integer, parameter :: LMO_3_INT_COMP        = 19
    integer, parameter :: LMO_2_INT_COMP        = 20
    integer, parameter :: LMO_34_EXT_LOOPS      = 21
    integer, parameter :: SEG_4_S               = 22
    integer, parameter :: SEG_3_S               = 23
    integer, parameter :: INT_3_SEG_ALL         = 24
    integer, parameter :: CHO_IABC              = 25
  
end module decider_symbols
