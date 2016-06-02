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
!*****************************************************************
!>  \brief Performs \f$ H_{diag} * c \f$
!>
!>  THIS ROUTINE COMPUTES THE CONTRIBUTIONS FROM THE
!>  PURELY DIAGONAL ELEMENTS AND THE ELEMENTS ASSOCIATED WITH DIFFERENT
!>  SPIN FUNCTIONS FROM THE SAME SPATIAL FUNCTION.
!>
!> \date 2013
!> \author Derek Walter
!> \author Arun Venkathnathan
!> \author Johannes M Dieterich (added BLAS statement where needed, removed multiplication by identity, cleaned up code)
!*****************************************************************
module diagtimesc_mod
    contains
subroutine diagtimesc(civec, sigmavec, loc_scr)

  use global_var_mod
  use graph_var_mod
  use spin_var_mod
  use ci_utilities_mod
  use scepmaker_mod
  use blocked_locks_mod
  use locist_mod
  use locist_var_mod,only:locist_scratch
#ifdef TIGER_USE_OMP
  use omp_lib
#endif

  implicit none
  
  type(locist_scratch)::loc_scr

  integer::sf1,sf2
  integer::ab_address,aa_address,a_address,v_address
  integer::interaction_index
  integer::v1_start,v2_start
  integer::a_start_1,a_start_2
  integer::ab_start_1,ab_start_2
  integer::aa_start_1,aa_start_2
  integer::rec_end 
  integer::internal_weight
  integer::length2C2

  integer::mode_1
  integer::singles
  integer::a,b
  integer::length
 
  type(blockedLockVectorType) :: civec, sigmavec
  real(kind=real8)::full_val
  real(kind=real8)::full_mat_sing
  real(kind=real8),dimension(num_external,num_external)::scep_ci_1,scep_ci_2,scep_sigma_1,scep_sigma_2
  real(kind=real8),dimension(num_external,num_external)::Hab_aa
  real(kind=real8)::Hab,Haa

  !// NOW WE ARE GOING TO LOOP OVER PATHS IN THE ORBITAL GRAPH AND 
  !// COMPUTE THE CONTRIBUTION FOR ALL THE FUNCTIONS ASSOCIATED
  !// WITH THE OCCUPATION PATTERN.

  rec_end = 0 
  call flush(ioOutput)


  !// REWIND THE DATA FILE NOW 
  rewind iodiag1 
  interaction_index = -1

  do 
     read(iodiag1,iostat=rec_end) interaction_index
     if (rec_end < 0) exit    
     if (rec_end > 0) then
        write(ioOutput,*) "error code: ", rec_end
        stop
     endif

     if (interaction_index == 2) then 
        read(iodiag1,iostat=rec_end) ab_address,aa_address,internal_weight,sf1,sf2,Hab,Haa,singles

        if (rec_end < 0) exit    
        if (rec_end > 0) then
           write(ioOutput,*) "error code: ", rec_end
           stop
        endif

        if (ab_address > 0 .and. aa_address > 0) then
        
           if (internal_index_vector2(internal_weight) < 0) cycle

           length = num_allowed_virtuals(internal_weight,"D")

           do a=1,length   
              do b=1,a-1
                 Hab_aa(a,b) = Hab
                 Hab_aa(b,a) = Hab
              enddo
              Hab_aa(a,a) = Haa
           enddo

           if (sf1 > fsn(singles)) then
              do a=1,length
                 Hab_aa(a,a) = 0.0d0
              enddo
           endif

           length2C2 = length*(length-1)/2
                                 
           mode_1 = -1
           if (sf1 <= fsn(singles)) mode_1 = 1

           ab_start_1 = ab_address + (sf1-1)*length2C2
           aa_start_1 = aa_address + (sf1-1)*length

           ab_start_2 = ab_address + (sf2-1)*length2C2
           aa_start_2 = aa_address + (sf2-1)*length

           call scepper_diag(scep_ci_2,civec%v,ab_start_2, &
                length, mode_1, civec%v,aa_start_2)

           scep_sigma_1(1:length,1:length) = 0.5d0*Hab_aa(1:length,1:length)*&
                scep_ci_2(1:length,1:length)

           call unscepper_diag(scep_sigma_1,sigmavec%v, ab_start_1, &
                length, mode_1, sigmavec%v,aa_start_1)


           if (sf1 /= sf2) then
           
              call scepper_diag(scep_ci_1,civec%v,ab_start_1, &
                   length, mode_1,civec%v,aa_start_1)
           
              !////////////////////////////////////////////////////////////////////////////////!
              !// NOW FOR THE SIGMA_{MU} !SEE THE DIMENSION OF SIGMA_NM2_MU
              scep_sigma_2(1:length,1:length) = &
                   0.5d0*Hab_aa(1:length,1:length)*scep_ci_1(1:length,1:length)

              !// NOW WE HAVE SIGMAVECTOR IN PAO BASIS 

              !////////////////////////////////////////////////////////////////////////////////!

              call unscepper_diag(scep_sigma_2,sigmavec%v, ab_start_2, &
                   length, mode_1, sigmavec%v,aa_start_2)

           endif
        elseif (ab_address > 0 .and. aa_address < 0) then
        
           if (internal_index_vector2(internal_weight) < 0) cycle

           length = num_allowed_virtuals(internal_weight,"D")

           do a=1,length
              do b=1,a-1 
                 Hab_aa(a,b) = Hab 
                 Hab_aa(b,a) = Hab
              enddo
              Hab_aa(a,a) = 0.0d0 
           enddo

           mode_1 = 1

           ab_start_1 = ab_address + (sf1-1)*length
           ab_start_2 = ab_address + (sf2-1)*length

           call scepper(scep_ci_2,civec%v, ab_start_2, length, mode_1)

           !////////////////////////////////////////////////////////////////////////////////!
           !// NOW FOR THE SIGMA_{LAMBDA} !SEE THE DIMENSION OF SIGMA_NM2_LAMBDA
           loc_scr%sigma_nm2_lambda(1:length,1:length) = &
                0.5d0*Hab_aa(1:length,1:length)*scep_ci_2(1:length,1:length)

           !// NOW WE HAVE SIGMAVECTOR IN PAO BASIS 

           call unscepper(loc_scr%sigma_nm2_lambda,sigmavec%v,ab_start_1,length, mode_1)



           if (sf1 /= sf2) then
           
              call scepper(scep_ci_1,civec%v, ab_start_1, length, mode_1)

              !////////////////////////////////////////////////////////////////////////////////!
              !// NOW FOR THE SIGMA_{MU} !SEE THE DIMENSION OF SIGMA_NM2_MU
              loc_scr%sigma_nm2_mu(1:length,1:length) = &
                   0.5d0*Hab_aa(1:length,1:length)*scep_ci_1(1:length,1:length)

              !// WE HAVE TO INCORPORATE THE PAO_OVERLAP MATRIX TWICE.
              !// FOR THE SIGMA_NM2_MU WE HAVE TO TAKE THE TRANSPOSE OF PAO OVERLAP
              !// MATRIX ...

              !// NOW WE HAVE SIGMAVECTOR IN PAO BASIS 
              call unscepper(loc_scr%sigma_nm2_mu,sigmavec%v, ab_start_2, length, mode_1)

           endif

        endif

     elseif (interaction_index  == 1) then 

        read(iodiag1,iostat=rec_end) a_address,internal_weight,sf1,sf2,full_mat_sing,singles
        if (rec_end < 0) exit    
        if (rec_end > 0) then
           write(ioOutput,*) "error code: ", rec_end
           stop
        endif

        if (internal_index_vector1(internal_weight) < 0) cycle

        length  = num_allowed_virtuals(internal_weight,"S")

        a_start_1 = a_address + (sf1-1)*length
        a_start_2 = a_address + (sf2-1)*length

        !////////////////////////////////////////////////////////////////////////////////!
        !// EVALUATE SIGMA_LAMBDA IN PAO BASIS !SEE THE DIMENSION OF SIGMA_NM2_LAMBDA
        call daxpy(length,full_mat_sing,civec%v(a_start_2),1,sigmavec%v(a_start_1),1)

        !////////////////////////////////////////////////////////////////////////////////!
        !// EVALUATE SIGMA_MU IN PAO BASIS !SEE THE DIMENSION OF SIGMA_NM2_MU
        if (sf1 /=sf2) then
           !// USE THE TRANSPOSE TO COMPUTE SIGMA_MU
           call daxpy(length,full_mat_sing,civec%v(a_start_1),1,sigmavec%v(a_start_2),1)
        endif


     elseif (interaction_index  == 0) then 

        read(iodiag1,iostat=rec_end) v_address,internal_weight,sf1,sf2,full_val,singles

        v1_start = v_address + sf1 - 1
        v2_start = v_address + sf2 - 1 


        sigmavec%v(v1_start)   =  sigmavec%v(v1_start) + full_val*civec%v(v2_start)


        if (sf1 /= sf2) then
           sigmavec%v(v2_start) = sigmavec%v(v2_start) + full_val*civec%v(v1_start)
        endif


     endif
     
  enddo
  
end subroutine diagtimesc
end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

