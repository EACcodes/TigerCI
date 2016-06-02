      subroutine spxfrm(u,s,h,w,n,m)
c     symmetric-packed matrix transformation		mar 83
c     written by ron shepard, argonne national lab
c     u = input real n by m matrix = u( n,m )
c     s = input real packed symmetric matrix of dimension n
c     h = outout real packed symmetric matrix of dimension m
c     w = scratch vector of length n
c     n = input scalar dimension of s and w
c     m = input scalar dimension of h
c     
c     computes h = ut.s.u
c     
      implicit double precision (a-h,o-z)
      double precision u(*),s(*),h(*),w(n)
c     
      if( n.le.0 .or. m.le.0 ) return
      ij = 1
      j1 = 1
      do j=1,m
        w(1) = s(1)*u(j1)
        if( n.le.1 ) then
c         
          do i=1,j
            h(ij) = u(i)*w(1)
            ij = ij+1
          end do
c           
        else
c           
          ix = 2
          do k=2,n
            w(k) = ddot(k,s(ix),1,u(j1),1)
            ix = ix+k
          end do
c             
              ix = 2
              do k=2,n
                km1 = k-1
                call daxpy(km1,u(j1+km1),s(ix),1,w,1)
                ix = ix+k
              end do
c               
              ix = 1
              do k=1,j
                h(ij) = ddot(n,u(ix),1,w,1)
                ij = ij + 1
                ix = ix+n
              end do
c                 
         endif
c               
        j1 = j1+n
      end do
      
      return
      end
