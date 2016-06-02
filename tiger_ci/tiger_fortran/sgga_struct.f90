! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module sgga_struct_mod
  use cholesky_structs
  use locist_var_mod
  use three_four_seg_var_mod
  use two_seg_var_mod
  implicit none

  type sgga_struct
     type(cholesky_data)::cho_data
     type(locist_scratch)::loc_scr
     type(threefourmod1vars)::mod1vars
     type(threefourmod2vars)::mod2vars
     type(twoModVars)::twoVars
  end type sgga_struct
  

end module sgga_struct_mod
