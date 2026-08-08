C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
       subroutine cwork5521(
     *                      y,na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,
     *                      ne1,ne2,
     *                      a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,
     *                      inda,indb,indc,indd,inde, 
     *
     *                      x1,nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,
     *                      nt1,nt2,
     *                      m1,m2,n1,n2,r1,r2,s1,s2,t1,t2,
     *                      indx1,
     * 
     *                      x2,ni1,ni2,nj1,nj2,i1,i2,j1,j2,
     *                      indx2,
     * 
     *                      cind, flopcount, scr1, scr2, scr3)
C
c      5 index output array  : y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,ne1:ne2)
c      5 index operand array : x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2,nt1,nt2)
c      2 index operand array : x2(ni1:ni2,nj1:nj2)
c      1 index contraction.
c
c   I. e., contract one of the indices of two four and two dimensional
C   arrays to produce a five dimensional array.
c--------------------------------------------------------------------------
      implicit none
      include 'trace.h'
      include 'saved_data.h'

      integer na1,na2,nb1,nb2,nc1,nc2,nd1,nd2,ne1,ne2,
     *        nm1,nm2,nn1,nn2,nr1,nr2,ns1,ns2,nt1,nt2,
     *        ni1,ni2,nj1,nj2
      integer a1,a2,b1,b2,c1,c2,d1,d2,
     *        e1,e2,
     *        n1,n2,m1,m2,r1,r2,s1,s2,t1,t2,
     *        i1,i2,j1,j2
      integer inda, indb, indc, indd, inde, indf, indx1(5), indx2(2)
      integer nrow,ncol,nsum
      integer cind(1), flopcount

      double precision y(na1:na2,nb1:nb2,nc1:nc2,nd1:nd2,ne1:ne2)
      double precision x1(nm1:nm2,nn1:nn2,nr1:nr2,ns1:ns2,nt1:nt2)
      double precision x2(ni1:ni2,nj1:nj2)
      double precision scr1(*), scr2(*), scr3(*)
      double precision etemp, thresh, max1, max2   
c
c Temporary arrays used in array permutations 
c 
c
c END Temporary arrays used in array permutations 

c     integer p1, p2, q1, q2
      integer na, nb, nc, nd, ne, nf, ng, nh, ni, nj, nk, nl, np, nq
      integer p0(5), q1(5), q2(2)

      integer i1a, i1b, i1c, i1d, i1p, i1q
      integer i2a, i2b, i2c, i2d, i2p, i2q
      integer a, b, c, d, e, f, g, h, i, j, l, p, q
      integer aa, bb, cc, dd, ee, ff, ii, jj, ll  
      integer yy, xx, y1, y2, y3, y4, yTOT, nTOT, inc, xinc    
      integer ystart, yend 

#ifdef __crayx1
      integer*8 lda, ldb, ldc, m, nn, k
#else
      integer lda, ldb, ldc, m, nn, k
#endif
      integer n
      character*1 op_a, op_b
      logical b2pq
      logical debug

      debug = cwork_debug
      flopcount = 0

      if (.false.) then

c   o assert x1 is contiguous
      if (m1.ne.nm1.or.n1.ne.nn1.or.r1.ne.nr1.or.s1.ne.ns1.or.
     &    m2.ne.nm2.or.n2.ne.nn2.or.r2.ne.nr2.or.s2.ne.ns2) then
         print *, '@CWORK5521: Assertion failed.'
         print *, '            The x1 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert x2 is contiguous
      if (i1.ne.ni1.or.j1.ne.nj1.or.i2.ne.ni2.or.j2.ne.nj2) then
         print *, '@CWORK5521: Assertion failed.'
         print *, '            The x2 tensor is not contiguous.'
         call abort_job()
      end if

c   o assert y is contiguous
      if (a1.ne.na1.or.b1.ne.nb1.or.c1.ne.nc1.or.d1.ne.nd1.or.
     &    a2.ne.na2.or.b2.ne.nb2.or.c2.ne.nc2.or.d2.ne.nd2.or. 
     &    e1.ne.ne1.or.e2.ne.ne2) then
         print *, '@CWORK5521: Assertion failed.'
         print *, '            The product tensor is not contiguous.'
         call abort_job()
      end if

      end if

c     go to 99 
c Try optimized version !!! VFL 

      do i = 1, 5
         p0(i) = i 
         q1(i) = 0 
      enddo
      do i = 1, 2
         q2(i) = 0 
      enddo 

      do i = 1, 5
         if (indx1(i) .eq. inda) q1(i) = p0(1) 
         if (indx1(i) .eq. indb) q1(i) = p0(2) 
         if (indx1(i) .eq. indc) q1(i) = p0(3) 
         if (indx1(i) .eq. indd) q1(i) = p0(4) 
         if (indx1(i) .eq. inde) q1(i) = p0(5) 
      enddo 

      do i = 1, 2
         if (indx2(i) .eq. inda) q1(i) = p0(1) 
         if (indx2(i) .eq. indb) q2(i) = p0(2) 
         if (indx2(i) .eq. indc) q2(i) = p0(3) 
         if (indx2(i) .eq. indd) q2(i) = p0(4) 
         if (indx2(i) .eq. inde) q2(i) = p0(5) 
      enddo 

      do i = 1, 5
      do j = 1, 2 
         if (indx1(i) .eq. indx2(j)) then 
            q1(i) = 7 ! indx1(i) 
            q2(j) = 7 ! indx1(i) 
         endif 
      enddo 
      enddo 

      If (q1(5) .ne. 7 .and. q2(1) .ne. 7) then
         Write(6,"(4a)") " Only one arrangement,", 
     *                    " x1(a,b,c.d,e,f)* x(f,g)",
     *                    " is available. Permute the indices before", 
     *                    " calling this contraction"
         call abort_job()
      Endif 

      nTOT = (e2-e1+1) 
     *     * (d2-d1+1) 
     *     * (c2-c1+1) 
     *     * (b2-b1+1) 
     *     * (a2-a1+1) 

#ifdef _DEBUG_LVL0
      Write(6,"(a,5i3)") "indas", inda,indb,indc,indd,inde
      Write(6,"(5i3)") (indx1(i),i=1,5)
      Write(6,"(5i3)") (q1(i),i=1,5)
      Write(6,"(6i3)") (indx2(i),i=1,2)
      Write(6,"(6i3)") (q2(i),i=1,2)
      Write(6,"(a,6i3)") "i1 and j2", i1,j2,ni1,ni2,nj1,nj2
      Write(6,"(6(1x,F10.5))") ((x2(i1+i,j1+j),i=0,ni2),j=0,nj2)
#endif 
C T(a,b,s1,c,s3) = X1(a,b,s1,c,s2) * x2(s2,s3)

      if ((q1(1).eq.1).and.(q1(2).eq.2).and.(q1(3).eq.3).and.
     +    (q1(4).eq.4).and.(q1(5).eq.7) .and.
     +    (q2(1).eq.7).and.(q2(2).eq.5)) then

          nrow = (m2-m1+1) * (n2-n1+1) * (r2-r1+1) * (s2-s1+1) 
          ncol = (j2-j1+1)
          nsum = (i2-i1+1)

          lda  = (nm2-nm1+1) * (nn2-nn1+1) * (nr2-nr1+1) * (ns2-ns1+1)
          ldb  = (ni2-ni1+1)
          ldc  = (na2-na1+1) * (nb2-nb1+1) * (nc2-nc1+1) * (nd2-nd1+1)

          Call dgemm('N','N',nrow,ncol,nsum,1.0D0,x1(m1,n1,r1,s1,t1),
     +                lda,x2(i1,j1),ldb,0.0D0,y(a1,b1,c1,d1,e1),ldc)

      Endif 

      Return 
      End

     
