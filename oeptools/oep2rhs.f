










      subroutine oep2rhs (
     &     nbas, nocc, naux,
     &     isdff, isdfw, isdww, isfdf, iswdw, isdffdiag, hybrid,
     &     ener, eval, oep1au, 
     &     auxpp, auxph, auxhh, f1pp, f1ph, f1hh,
     &     ppph, pphh, phph, phhh, intpp, intph, inthh, intr12, rhs,
     &     dor12)
c
c This routine calculate the right hand side of the second-order OEP
c equation and the second-order energy.
c
c For the detail on the first- and second-order OEP equation see
c document REPORTS/OEP2.DERIVATION.TEX
c
c To make the routine more universial the contribution of the each
c diagram can be neglected by setting the corresponding flag ISxxx to
c FALSE. E.g., by the exchange-only calculation can be emulating by
c setting all the flags to FALSE.
c
c In: nbas, nocc, naux, eval, auxpp, auxph, auxhh, f1pp, f1ph, f1hh,
c     ppph, pphh, phph, phhh
c Scr: intpp, intph, inthh,
c Out: rhs
c
c Igor Schweigert, Jan 2004
c $Id: oep2rhs.FPP,v 1.4 2008/06/06 18:09:17 taube Exp $
c
      implicit none
c
c     Arguments
c
      logical
     &     isdff, isdfw, isdww, isfdf, iswdw, isdffdiag, dor12
c     
      integer
     &     nbas, nocc, naux
c
      double precision
     &     hybrid, ener, eval (nbas), oep1au (naux), 
     &     auxpp (*), auxph (*), auxhh (*),
     &      f1pp (*),  f1ph (*),  f1hh (*),
     &      ppph (*),  pphh (*),  phph (*), phhh (*),
     &     intpp (*), intph (*), inthh (*), rhs (*)
c
c     Local variables
c
      integer
     &     n, n1, n2, n3, n4, a, b, c, i1, i2, i3, i4
c
      double precision
     &     e1, e2, x
c     R12 variables
c     Need 4 byte integers to connect with Noga's code
      integer*4
     &     nh4,np4
c
      integer
     &     unit_mos 
c     
      double precision           
     &     er12,r12fact,intr12(*)
c
c     Define the one- and two-electron array offset operators
c
c
c This file defines two-electron integral index statement functions of
c various kinds.  It must be included right after the variables
c declarations and before the first executable statement, since it
c contains the variable declarations, statement function definitions and
c variable definitions (that are executable statements)
c
c For the operators to work, two macros N_BAS and N_OCC have to be
c defined to the names of the variables that give the number of basis
c functions and the number of occupied orbitals correspondingly. So a
c routine that includes these definitions must have these variables
c declared and defined. This piece of code also affects the caller name
c space since it needs the statement function dummy arguments to be
c declared. The implementation could be much clear if Fortran allowed
c the external functions to be specified as inline. Most of the modern
c compiler will do that automatically during the optimization step, but
c it is not guaranteed. 
c
c See also the documentation file OEPINTIND.TEX for the description of
c the offset operators for the two-electron integral indices.
c 
c Igor Schweigert, Jan 2004
c $Id: oepintind.com,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
c
c     Include machine-dependent definitions for VMOL-related operators.
c     


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c
c     Declare the operators and the dummy variables. This will affect
c     the namespace of the caller, so be careful.
c     
      integer
     &     i, j, k, l, n_a, n_p, n_h,
     &     i_aa, i_ma, i_ah, i_pa, i_hh, i_pp, i_ph,
     &     n_aa, n_ma, n_ah, n_pa, n_hh, n_pp, n_ph,
     &     i_aaah, i_paah, i_ppah, i_phah, 
     &     n_aaah, n_paah, n_ppah, n_phah, 
     &     i_ppph, i_pphh, i_phph, i_phhh,
     &     n_ppph, n_pphh, n_phph, n_phhh,
     &     i1_aa, i1_pp, i1_hh, i1_ma, i1_ah, i1_pa, i1_ph,
     &     i2_aa, i2_pp, i2_hh, i2_ma, i2_ah, i2_pa, i2_ph,
     &     i1_aaah, i1_paah, i1_ppah, i1_phah, 
     &     i2_aaah, i2_paah, i2_ppah, i2_phah,
     &     i3_aaah, i3_paah, i3_ppah, i3_phah, 
     &     i4_aaah, i4_paah, i4_ppah, i4_phah,
     &     i1_ppph, i1_pphh, i1_phph, i1_phhh,
     &     i2_ppph, i2_pphh, i2_phph, i2_phhh,
     &     i3_ppph, i3_pphh, i3_phph, i3_phhh,
     &     i4_ppph, i4_pphh, i4_phph, i4_phhh,
     &     iupki, iupkj, iupkk, iupkl,
     &     i_aa_cr, i_pp_cr, i_hh_cr,
     &     i_aaah_cr, i_ppah_cr,
     &     i_ppph_cr, i_pphh_cr, i_phph_cr, i_phhh_cr

c
c     Define the operators that pack two orbital indices into the
c     combined index.
c     
      i_aa (i, j) = i + (j*(j-1))/2
      i_hh (i, j) = i + (j*(j-1))/2
      i_pp (i, j) = i + (j*(j-1))/2
c
      i_ma (i, j) = i + n_a * (j-1)
      i_ah (i, j) = i + n_a * (j-1)
      i_pa (i, j) = i + n_p * (j-1)
      i_ph (i, j) = i + n_p * (j-1)
c
c     Define the operators that pack four orbital indices into the
c     combined index.
c 
      i_aaah (i, j, k, l) = i_aa (i, j) + n_aa * (i_ah (k, l) - 1)
      i_paah (i, j, k, l) = i_pa (i, j) + n_pa * (i_ah (k, l) - 1)
      i_ppah (i, j, k, l) = i_pp (i, j) + n_pp * (i_ah (k, l) - 1)
      i_phah (i, j, k, l) = i_ph (i, j) + n_ph * (i_ah (k, l) - 1)
      i_ppph (i, j, k, l) = i_pp (i, j) + n_pp * (i_ph (k, l) - 1)
      i_pphh (i, j, k, l) = i_pp (i, j) + n_pp * (i_hh (k, l) - 1)
      i_phph (i, j, k, l) = i_aa (i_ph (i, j), i_ph (k, l))
      i_phhh (i, j, k, l) = i_ph (i, j) + n_ph * (i_hh (k, l) - 1)
c
c     Define the operators that pack two orbital indices into the
c     combined index. These operators ensure the proper range of
c     indices, and hence a bit slower. Some operators do not impose the
c     inequality condition on their arguments, so they do not have the
c     "_cr" versions.
c
      i_aa_cr (i, j) = min(i,j) + (max(i,j)*(max(i,j)-1))/2
      i_hh_cr (i, j) = min(i,j) + (max(i,j)*(max(i,j)-1))/2
      i_pp_cr (i, j) = min(i,j) + (max(i,j)*(max(i,j)-1))/2
      i_aaah_cr (i, j, k, l) = i_aa_cr (i, j) + n_aa * (i_ah (k, l) - 1)
      i_ppah_cr (i, j, k, l) = i_pp_cr (i, j) + n_pp * (i_ah (k, l) - 1)
      i_ppph_cr (i, j, k, l) = i_pp_cr (i, j) + n_pp * (i_ph (k, l) - 1)
      i_pphh_cr (i, j, k, l) = i_pp_cr (i, j) + n_pp * (i_hh_cr(k,l)- 1)
      i_phph_cr (i, j, k, l) = i_aa_cr (i_ph (i, j), i_ph (k, l))
      i_phhh_cr (i, j, k, l) = i_ph (i, j) + n_ph * (i_hh_cr (k, l) - 1)
c     
c     Define the operators that unpack a combined index into the two
c     orbital indices.
c     
      i2_aa (i) = 1 + (-1 + int (dsqrt (8.d0*i+0.999d0)))/2
      i2_pp (i) = 1 + (-1 + int (dsqrt (8.d0*i+0.999d0)))/2
      i2_hh (i) = 1 + (-1 + int (dsqrt (8.d0*i+0.999d0)))/2
      i2_ma (i) = (i-1) / n_a + 1
      i2_ah (i) = (i-1) / n_a + 1
      i2_pa (i) = (i-1) / n_p + 1
      i2_ph (i) = (i-1) / n_p + 1
      i1_aa (i, j) = i - (j*(j-1))/2
      i1_pp (i, j) = i - (j*(j-1))/2
      i1_hh (i, j) = i - (j*(j-1))/2
      i1_ma (i, j) = i - n_a * (j-1)
      i1_ah (i, j) = i - n_a * (j-1)
      i1_pa (i, j) = i - n_p * (j-1)
      i1_ph (i, j) = i - n_p * (j-1)
c
c     Define the operators that unpack a combined index into the four
c     orbital indices.
c     
      i4_aaah (i) = i2_ah ((i-1)/n_aa + 1)
      i4_paah (i) = i2_ah ((i-1)/n_pa + 1)
      i4_ppah (i) = i2_ah ((i-1)/n_pp + 1)
      i4_phah (i) = i2_ah ((i-1)/n_ph + 1)
      i4_ppph (i) = i2_ph ((i-1)/n_pp + 1)
      i4_pphh (i) = i2_hh ((i-1)/n_pp + 1)
      i4_phph (i) = i2_ph (i2_aa (i) )
      i4_phhh (i) = i2_hh ((i-1)/n_ph + 1)
c 
      i3_aaah (i, j) = i1_ah ((i-1)/n_aa + 1, j)
      i3_paah (i, j) = i1_ah ((i-1)/n_pa + 1, j)
      i3_ppah (i, j) = i1_ah ((i-1)/n_pp + 1, j)
      i3_phah (i, j) = i1_ah ((i-1)/n_ph + 1, j)
      i3_ppph (i, j) = i1_ph ((i-1)/n_pp + 1, j)
      i3_pphh (i, j) = i1_hh ((i-1)/n_pp + 1, j)
      i3_phph (i, j) = i1_ph (i2_aa (i),  j)
      i3_phhh (i, j) = i1_hh ((i-1)/n_ph + 1, j)
c 
      i2_aaah (i, j, k) = i2_aa (i - n_aa * (i_ah (j, k) - 1))
      i2_paah (i, j, k) = i2_pa (i - n_pa * (i_ah (j, k) - 1))
      i2_ppah (i, j, k) = i2_pp (i - n_pp * (i_ah (j, k) - 1))
      i2_phah (i, j, k) = i2_ph (i - n_ph * (i_ah (j, k) - 1))
      i2_ppph (i, j, k) = i2_pp (i - n_pp * (i_ph (j, k) - 1))
      i2_pphh (i, j, k) = i2_pp (i - n_pp * (i_hh (j, k) - 1))
      i2_phph (i, j, k) = i2_ph (i1_aa (i, i_ph (j,k)))
      i2_phhh (i, j, k) = i2_ph (i - n_ph * (i_hh (j, k) - 1))
c 
      i1_aaah (i, j, k, l) = i1_aa (i - n_aa * (i_ah (k, l) - 1), j)
      i1_paah (i, j, k, l) = i1_pa (i - n_pa * (i_ah (k, l) - 1), j)
      i1_ppah (i, j, k, l) = i1_pp (i - n_pp * (i_ah (k, l) - 1), j)
      i1_phah (i, j, k, l) = i1_ph (i - n_ph * (i_ah (k, l) - 1), j)
      i1_ppph (i, j, k, l) = i1_pp (i - n_pp * (i_ph (k, l) - 1), j)
      i1_pphh (i, j, k, l) = i1_pp (i - n_pp * (i_hh (k, l) - 1), j)
      i1_phph (i, j, k, l) = i1_ph (i1_aa (i, i_ph (k, l)), j)
      i1_phhh (i, j, k, l) = i1_ph (i - n_ph * (i_hh (k, l) - 1), j)
c
c     Define the operators that unpack the VMOL index into the four 2e
c     integral indices.
c
      iupki (i) = iand (i,ialone)
      iupkj (i) = iand (ishft(i,-ibitwd),ialone)
      iupkk (i) = iand (ishft(i,-2*ibitwd),ialone)
      iupkl (i) = iand (ishft(i,-3*ibitwd),ialone)

c
c This file defines two-electron integral index statement functions of
c various kinds.  It must be included right after the variables
c declarations and before the first executable statement, since it
c contains the variable declarations, statement function definitions and
c variable definitions (that are executable statements)
c
c For the operators to work, two macros nbas and nocc have to be
c defined to the names of the variables that give the number of basis
c functions and the number of occupied orbitals correspondingly. So a
c routine that includes these definitions must have these variables
c declared and defined. This piece of code also affects the caller name
c space since it needs the statement function dummy arguments to be
c declared. The implementation could be much clear if Fortran allowed
c the external functions to be specified as inline. Most of the modern
c compiler will do that automatically during the optimization step, but
c it is not guaranteed. 
c
c See also the documentation file OEPINTIND.TEX for the description of
c the offset operators for the two-electron integral indices.
c 
c Igor Schweigert, Jan 2004
c $Id: oepintind.com,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
c
c     Define the sizes of arrays. Note these definitions rely on two
c     macros that has to be defined, nbas and nocc
c     
      n_a = nbas
      n_h = nocc
      n_p = n_a - n_h
c     
      n_aa = i_aa (n_a, n_a)
      n_pp = i_pp (n_p, n_p)
      n_hh = i_hh (n_h, n_h)
      n_ma = i_ma (n_a, n_a)
      n_ah = i_ah (n_a, n_h)
      n_pa = i_pa (n_p, n_a)
      n_ph = i_ph (n_p, n_h)
c 
      n_aaah = i_aaah (n_a, n_a, n_a, n_h)
      n_paah = i_paah (n_p, n_a, n_a, n_h)
      n_ppah = i_ppah (n_p, n_p, n_a, n_h)
      n_phah = i_phah (n_p, n_h, n_a, n_h)
      n_ppph = i_ppph (n_p, n_p, n_p, n_h)
      n_pphh = i_pphh (n_p, n_p, n_h, n_h)
      n_phph = i_phph (n_p, n_h, n_p, n_h)
      n_phhh = i_phhh (n_p, n_h, n_h, n_h)
c

c
c     Messages to the user
c     
 5905 format (
     &     t3,'@OEP2RHS-W: Reference eigenvalues are degenerate.',/,
     &     t3,'The HOMO-LUMO gap is ',f20.10)
c     
c     Calculate the MOMEs of the one-particle part of the perturbation
c     by subtracting the first-order OEP MOMEs from non-local exchange
c     MOMEs and scale by the 1 minus the fraction of the non-local
c     exchange.
c
      call dgemm (
     &     'n', 'n', n_pp, 1, naux,
     &     -1.d0+hybrid, auxpp, n_pp,
     &     oep1au, naux,
     &     1.d0-hybrid, f1pp, n_pp)
      call dgemm (
     &     'n', 'n', n_ph, 1, naux,
     &     -1.d0+hybrid, auxph, n_ph,
     &     oep1au, naux,
     &     1.d0-hybrid, f1ph, n_ph)
      call dgemm (
     &     'n', 'n', n_hh, 1, naux,
     &     -1.d0+hybrid, auxhh, n_hh,
     &     oep1au, naux,
     &     1.d0-hybrid, f1hh, n_hh)
c     
c     Calculate the second-order contribution to the total energy: <>
c     Calculate the singles contibutions. <> Calculate the doubles
c     contribution.
c     
      e1 = ener
      do n1 = 1, n_ph
         i2 = i2_ph (n1)
         i1 = i1_ph (n1, i2)
         ener = ener +
     &        2.d0 * f1ph (n1) * f1ph (n1) / (eval (i2) - eval (i1+n_h))
      enddo
c$$$      print *,'E2s = ', ener - e1, ' GAP = ', eval (n_h) - eval (n_h+1)
      if (eval (n_h) - eval (n_h+1).gt.-1.d-3)
     &     write (6, 5905) eval (n_h) - eval (n_h+1)
c      e1r12=ener
c      do n1 = 1, n_hh
c         ener = ener + 2.d0 * v1(n1) * c1(n1)
c      end do
c$$$      print *,'E2r12S = ', ener - e1r12
c     
      e2 = ener
      do n1 = 1, n_phph
         i4 = i4_phph (n1)
         i3 = i3_phph (n1, i4)
         i2 = i2_phph (n1, i3, i4)
         i1 = i1_phph (n1, i2, i3, i4)
         e1 = eval (i2) + eval (i4) - eval (i1+n_h) - eval (i3+n_h)
         n2 = i_phph_cr (i1, i4, i3, i2)
         ener = ener +
     &        phph (n1) * (2.d0 * phph (n1) - phph (n2)) / e1
         if (i_ph (i1, i2) .ne. i_ph (i3, i4)) then
            n2 = i_phph_cr (i3, i2, i1, i4)
            ener = ener +
     &           phph (n1) * (2.d0 * phph (n1) - phph (n2)) /e1
         endif
      enddo
c$$$      print *, 'E2d = ', ener - e2
c      e2r12 = ener
c      do n1 = 1, n_hhhh
c      do n1 = 1, n_phph
c         i4 = i4_phph (n1)
c         i3 = i3_phph (n1, i4)
c         i2 = i2_phph (n1, i3, i4)
c         i1 = i1_phph (n1, i2, i3, i4)
c         e1 = eval (i2) + eval (i4) - eval (i1+n_h) - eval (i3+n_h)
c         n2 = i_phph_cr (i1, i4, i3, i2)
c         ener = ener +
c     &        phph (n1) * (2.d0 * phph (n1) - phph (n2)) / e1
c         if (i_ph (i1, i2) .ne. i_ph (i3, i4)) then
c            n2 = i_phph_cr (i3, i2, i1, i4)
c            ener = ener +
c     &           phph (n1) * (2.d0 * phph (n1) - phph (n2)) /e1
c         endif
c      enddo
c$$$      print *, 'E2d = ', ener - e2
c     
c     Initialize the intermediates.
c
      call zero (intpp, n_pp)
      call zero (intph, n_ph)
      call zero (inthh, n_hh)
c
c     DFF - SINGLES PH. Calculate the contribution of the DFF diagram.
c     
      if (.not.isdff .or. hybrid.eq.1d0) goto 190
      do n = 1, n_ph
         i = i2_ph (n)
         b = i1_ph (n, i)
         do a = 1, n_p
            if (.not.isdffdiag.and.a.eq.b) goto 110
            e1 = eval (i) - eval (b+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_pp_cr (b, a)
            n2 = i_ph    (a, i)
            intph (n) = intph (n) + 4.d0 * 
     &           f1pp (n1) * f1ph (n2) /
     &           (e1 * e2)
 110        continue
         enddo
      enddo
c
      do n = 1, n_ph
         j = i2_ph (n)
         a = i1_ph (n, j)
         do i = 1, n_h
            if (.not.isdffdiag.and.i.eq.j) goto 120
            e1 = eval (j) - eval (a+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_hh_cr (j, i)
            n2 = i_ph    (a, i)
            intph (n) = intph (n) - 4.d0 *
     &           f1hh (n1) * f1ph (n2) /
     &           (e1 * e2)
 120        continue
         enddo
      enddo     
 190  continue
c
c     DFW - SINGLE and DOUBLES PH. Calculate the contribution of the DFW
c     and DWF diagrams
c     
      if (.not.isdfw .or. hybrid.eq.1d0) goto 290
      do n = 1, n_ph
         j = i2_ph (n)
         b = i1_ph (n, j)
         do a = 1, n_p
            do i = 1, n_h
               e1 = eval (j) - eval (b+n_h)
               e2 = eval (i) - eval (a+n_h)
               n1 = i_ph (a, i)
               n2 = i_pphh_cr (b, a, i, j)
               n3 = i_phph_cr (b, i, a, j)
               intph (n) = intph (n) - 4.d0 *
     &              f1ph (n1) * (pphh (n2) + phph (n3)) /
     &              (e1 * e2)
            enddo
         enddo
      enddo
 290  continue
c
c     DWW - DOUBLES PH. Calculate the contribution of the DWW diagram.
c
      if (.not.isdww) goto 390
      do n = 1, n_ph
         i = i2_ph (n)
         c = i1_ph (n, i)
         do a = 1, n_p
            do b = 1, a
               do j = 1, n_h
                  e1 = eval (i) - eval (c+n_h)
                  e2 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_ppph_cr (a, c, b, j)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_ppph_cr (b, c, a, j)
                  n4 = i_phph_cr (a, j, b, i)
                  if (b.eq.a) e1 = 2.d0 * e1
                  intph (n) = intph (n) + 4.d0 * (
     &                 ppph (n1) * phph (n2) * 2.d0 -
     &                 ppph (n1) * phph (n4) -
     &                 ppph (n3) * phph (n2) +
     &                 ppph (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c
      do n = 1, n_ph
         k = i2_ph (n)
         a = i1_ph (n, k)
         do b = 1, n_p
            do i = 1, n_h
               do j = 1, i
                  e1 = eval (k) - eval (a+n_h)
                  e2 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_phhh_cr (b, j, k, i)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_phhh_cr (b, i, k, j)
                  n4 = i_phph_cr (a, j, b, i)
                  if (j.eq.i) e1 = 2.d0 * e1
                  intph (n) = intph (n) - 4.d0 * (
     &                 phhh (n1) * phph (n2) * 2.d0 -
     &                 phhh (n1) * phph (n4) -
     &                 phhh (n3) * phph (n2) +
     &                 phhh (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c     
 390  continue

c
c     FDF - SINGLES PP, HH. Calculate the contribution of the FDF diagram
c     
      if (.not.isfdf .or. hybrid.eq.1d0) goto 490
      do n = 1, n_pp
         a = i2_pp (n)
         b = i1_pp (n, a)
         do i = 1, n_h
            e1 = eval (i) - eval (b+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_ph (b, i)
            n2 = i_ph (a, i)
            intpp (n) = intpp (n) + 2.d0 * 
     &           f1ph (n1) * f1ph (n2) /
     &           (e1 * e2)
         enddo
      enddo
c     
      do n = 1, n_hh
         i = i2_hh (n)
         j = i1_hh (n, i)
         do a = 1, n_p
            e1 = eval (j) - eval (a+n_h)
            e2 = eval (i) - eval (a+n_h)
            n1 = i_ph (a, j)
            n2 = i_ph (a, i)
            inthh (n) = inthh (n) - 2.d0 * 
     &           f1ph (n1) * f1ph (n2) /
     &           (e1 * e2)
         enddo
      enddo
 490  continue

c
c     WDW - DOUBLES PP, HH. Calculate the contribution of the WDW diagram.
c
      if (.not.iswdw) goto 590
      do n = 1, n_pp
         a = i2_pp (n)
         c = i1_pp (n, a)
         do b = 1, n_p
            do i = 1, n_h
               do j = 1, i
                  e1 = eval (i) + eval (j) - eval (c+n_h) - eval (b+n_h)
                  e2 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_phph_cr (c, i, b, j)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_phph_cr (c, j, b, i)
                  n4 = i_phph_cr (a, j, b, i)
                  if (i.eq.j) e1 = 2.d0 * e1
                  intpp (n) = intpp (n) + 2.d0 * (
     &                 phph (n1) * phph (n2) * 2.d0 -
     &                 phph (n1) * phph (n4) -
     &                 phph (n3) * phph (n2) +
     &                 phph (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c
      do n = 1, n_hh
         i = i2_hh (n)
         k = i1_hh (n, i)
         do j = 1, n_h
            do a = 1, n_p
               do b = 1, a
                  e1 = eval (i) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  e2 = eval (k) + eval (j) - eval (a+n_h) - eval (b+n_h)
                  n1 = i_phph_cr (a, k, b, j)
                  n2 = i_phph_cr (a, i, b, j)
                  n3 = i_phph_cr (a, j, b, k)
                  n4 = i_phph_cr (a, j, b, i)
                  if (a.eq.b) e1 = 2.d0 * e1
                  inthh (n) = inthh (n) - 2.d0 * (
     &                 phph (n1) * phph (n2) * 2.d0 -
     &                 phph (n1) * phph (n4) -
     &                 phph (n3) * phph (n2) +
     &                 phph (n3) * phph (n4) * 2.d0) /
     &                 (e1 * e2)
               enddo
            enddo
         enddo
      enddo
c$$$      write (6,*) 'INTHH 1: ', inthh (1)
c$$$      inthh (1) = 0.d0
c$$$      do a = 1, n_p
c$$$         do b = 1, n_p
c$$$            n = i_phph_cr (a, 1, b, 1)
c$$$            e1 = 2.d0 * eval (1) - eval (a+1) - eval (b+1)
c$$$            inthh (1) = inthh (1) - 2.d0 * phph (n)**2 / e1**2
c$$$         enddo
c$$$      enddo
c$$$      write (6,*) 'INTHH 1: ', inthh (1)
 590  continue
      
      if (dor12) then
C       R12 contribution
C       Completely separate read of R12 OV_OO File test      
        unit_mos=8
        open(unit=unit_mos,file='OV_OO_BLOCK_F12',ACCESS='SEQUENTIAL',
     &       STATUS='OLD',FORM='UNFORMATTED')
        read(unit=unit_mos) er12,nh4,np4
        n1=nh4
        n2=np4
        read(unit=unit_mos) (intr12(i),i=1,n1*n2)
        close(unit=unit_mos)
C       R12 array is stored hole-particle
C       OEP expects particle-hole
C       Need to check proper factor... don't know exactly how to check this yet
        r12fact=1.0d0      
        do n = 1, n_ph
C n1 is the index of the transposed matrix               
          i = i2_ph (n)
          a = i1_ph (n, i)
          n1 = i + n_h * (a-1)   
          intph(n) = intph(n)+r12fact*intr12(n1)
        enddo
      endif       
       
c
c     Calculate the OEP r.h.s. in auxiliary basis by contracting the RHS
c     intermediate with the auxiliary MOMEs: <> Scale the diagonal
c     elements of auxiliary PP and HH MOMEs by a factor of one-half. <>
c     Contract the intermediates with the auxiliary MOMEs, type by
c     type. <> Scale the diagonal elements of auxiliary PP and HH MOMEs
c     by a factor of two.
c
c     Note that the fact the full matrix is twice as long as PP, PH, HH
c     triangular matrices is taken into account by factor of two in the
c     DGEMMs. Hence, in order to prevent the double-counting of the
c     diagonal elements (PP and HH), they are scaled by a factor of
c     one-half.
c
      do i = 1, naux
         do j = 1, n_p
            n = i_pp (j, j) + n_pp * (i-1)
            auxpp (n) = .5d0 * auxpp (n)
         enddo
         do j = 1, n_h
            n = i_hh (j, j) + n_hh * (i-1)
            auxhh (n) = .5d0 * auxhh (n)
         enddo
      enddo
c
      call dgemm (
     &     'n', 'n', 1, naux, n_pp,
     &     2.d0, intpp, 1,
     &     auxpp, n_pp,
     &     0.d0, rhs, 1)
c
      call dgemm (
     &     'n', 'n', 1, naux, n_ph,
     &     1.d0, intph, 1,
     &     auxph, n_ph,
     &     1.d0, rhs, 1)
c
      call dgemm (
     &     'n', 'n', 1, naux, n_hh,
     &     2.d0, inthh, 1,
     &     auxhh, n_hh,
     &     1.d0, rhs, 1)
c     
      do i = 1, naux
         do j = 1, n_p
            n = i_pp (j, j) + n_pp * (i-1)
            auxpp (n) = 2.d0 * auxpp (n)
         enddo
         do j = 1, n_h
            n = i_hh (j, j) + n_hh * (i-1)
            auxhh (n) = 2.d0 * auxhh (n)
         enddo
      enddo
c     
      return
      end
