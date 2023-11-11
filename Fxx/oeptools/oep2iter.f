










      subroutine oep2iter (
     &     nbas, nocc,
     &     icore, maxmem, usdmem,
     &     scfh1e, scfevl, scfevc, scfden, scfham, scfenr, dor12)
c
c This routine solves the Optimized Effective Potential equations (in
c the first order and second order) to get the local OEP exchange and
c correlation potentials.
c
c The OEP solution will be obtained as an expansion over auxiliary basis
c functions. The OEP local exchange and/or correlation AOMEs will be
c then calculated and added to core-Hamiltonian AOMEs (supplied in
c SCFH1E) and Coulomb AOMEs to form the SCF Hamiltonian (passed back
c through SCFHAM). The corresponding OEP energy is added to the
c nuclear-nuclear repulsion energy (supplied in SCFENR) and passed back
c through SCFENR.
c
c Since we need a lot of different scratch arrays to solve the OEP
c equation, the routine itself allocates and releases memory, using
c available core passed in through ICORE. The amount of memory available
c is given by MAXMEM. That is all the routine needs to know, however, in
c a case of insufficient memory, we need to let user know how much total
c memory is needed, so the amount of memory already used in VSCF is
c given by USDMEM.
c
c In:  nbas, nocc, naux, maxmem, usdmem,
c      scfh1e, scfevl, scfevc, scfden, scfenr, scfham
c Scr: icore
c Out: scfenr and scfham are modified
c
c Igor Schweigert, Jan 2003 
c $Id: oep2iter.FPP,v 1.3 2008/06/06 18:09:17 taube Exp $
c
      implicit none
c
c     Common block
c
c symm2.com : begin

c This is initialized in vscf/symsiz.

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end


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
c This common block defines the OEP variables of general scope, mostly
c flags and their default values.
c
c Igor Schweigert, Feb 2004
c $Id: oep.com,v 1.8 2008/06/06 18:09:17 taube Exp $
c 
c
      logical
     &     OEP1, OEP2, OEP1ON, OEP2ON, OEP2ALWAYS, OEP2ITR, 
     &     OEPDFF, OEPDFW, OEPDWW, OEPFDF, OEPWDW, BOEPR12,
     &     R12ON
c     
      integer
     &     OEP1ONTHR, OEP1ONTHRDEF, OEP2ONTHR, OEP2ONTHRDEF,
     &     OEP2OFFTHR, OEP2OFFTHRDEF, OEP2UPDATE, OEP2UPDATEDEF,
     &     OEP2ALWAYSTHR, OEP2ALWAYSTHRDEF, NAUX,
     &     OEPSVDTHR, OEPSVDTHRDEF, OEP1HYBRIDINT, OEP1HYBRIDDEF,
     &     R12ONTHR, R12ONTHRDEF
c
      double precision OEP1HYBRID
c     
      character*(50)
     &     OEPH0, OEPH0DEF, OEPAXP, OEPAXPDEF
     
c      character*1
c     &     asv_null  
c     
      common /OEPCONFIG/
     &     OEP1, OEP2, OEP1ON, OEP2ON, OEP2ALWAYS, OEP2ITR,
     &     BOEPR12, R12ON, OEPDFF, OEPDFW, OEPDWW, OEPFDF, OEPWDW,
     &     OEP1ONTHR, OEP2ONTHR, OEP2OFFTHR, OEP2ALWAYSTHR, 
     &     OEPSVDTHR, NAUX, OEP2UPDATE, R12ONTHR,OEP1HYBRID, 
     &     OEPH0, OEPAXP
c     
      parameter (OEP1ONTHRDEF  = 1)
      parameter (OEP2ONTHRDEF  = 2)
      parameter (OEP2OFFTHRDEF = 20)
      parameter (OEP2ALWAYSTHRDEF = 6)
      parameter (OEP2UPDATEDEF = 1)
      parameter (OEPSVDTHRDEF       = 8)
      parameter (OEP1HYBRIDDEF = 0)
      parameter (R12ONTHRDEF  = 2)
c
      data
     &     OEPH0DEF, OEPAXPDEF   /'OEP', 'NONE'/
c     
      
c
c     The input parameters
c
      integer
     &     nbas, nocc,
     &     icore (*), maxmem, usdmem
c
      double precision
     &     scfh1e (itriln (nirrep+1)),
     &     scfevl (nbas), scfevc (nbas, nbas),
     &     scfden (itriln (nirrep+1)),
     &     scfham (itriln (nirrep+1)), scfenr
     
      logical 
     &      dor12
c
c     Local variables
c
      logical
     &     diagonal
c     
      integer
     &     x, n, 
     &     i010, i020, i030, i040, i050, i060, itop, coretop,
     &     auxaa, auxpp, auxph, auxhh,
     &     hfxaa, hfxpp, hfxph, hfxhh,
     &     ippph, ipphh, iphph, iphhh,
     &     resinv, oep1au, oep2au, refevl
c
      integer lnbuf, luint
      parameter (lnbuf=600)
      parameter (luint=10)
c
      double precision homo
c
c     Define the one- and two-electron integrals array sizes.
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
c     Initialize the top of the core at the begging of the ICORE array
c     (it has been already been shifted in VSCF.)
c
      coretop = 1
c
c     Calculate all the MO integrals (auxiliary and two-electron) needed
c     to solve the OEP2 equation: <> Allocate memory for the auxiliary,
c     one- and two-electron MO integrals and the scratch arrays
c     necessary to read and transform the AO integrals. <> Check if
c     there is enough memory available. <> Initialize the SCF
c     Hamiltonian with the core Hamiltonian. <> Read and transform
c     AOMEs.
c
c     Memory pointers legend: i010 - PAAH, i020 - AAAH, i030 - PPAH,
c     i040 - PHAH, i050 - integral buffer, i060 - index buffer.
c
c
      auxaa = coretop
      auxpp = auxaa + iintfp * n_aa * naux
      auxph = auxpp + iintfp * n_pp * naux 
      auxhh = auxph + iintfp * n_ph * naux 
      hfxaa = auxhh + iintfp * n_hh * naux 
      hfxpp = hfxaa + iintfp * n_aa
      hfxph = hfxpp + iintfp * n_pp
      hfxhh = hfxph + iintfp * n_ph
      ippph = hfxhh + iintfp * n_hh
      ipphh = ippph + iintfp * n_ppph
      iphph = ipphh + iintfp * n_pphh
      iphhh = iphph + iintfp * n_phph
      coretop = iphhh + iintfp * n_phhh
      i010  = ippph
      i020  = i010  + iintfp * n_paah
      i030  = i020  
      i040  = i030  + iintfp * n_ppah
      i050  = i020  + iintfp * n_aaah
      i060  = i050  + iintfp * lnbuf
      itop  = i060  + lnbuf
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Transformation of the auxiliary and 
     &     two-electron AO integrals',
     &     itop-1+usdmem, maxmem+usdmem)
c     
      if (OEPH0 (1:4) .eq. 'SEMI') then
         call oepsemcan (
     &        nbas, nocc, luint, lnbuf, icore (i060),
     &        icore (i050), scfevc, scfden, scfh1e,
     &        icore (hfxaa), icore (hfxpp), icore (hfxhh),
     &        icore (i010), icore (i020))
      endif
c     
      call dcopy (n_aa, scfh1e, 1, scfham, 1)
c 
      call oep2ints (
     &     .true., nbas, nocc, naux, luint, lnbuf, OEPAXP, icore (i060), 
     &     icore (i050), scfevc, scfham, 
     &     icore (auxaa), icore (auxpp), icore (auxph), icore (auxhh),
     &     icore (hfxaa), icore (hfxpp), icore (hfxph), icore (hfxhh),
     &     icore (i010 ), icore (i020 ), icore (i030 ), icore (i040),
     &     icore (ippph), icore (ipphh), icore (iphph), icore (iphhh))
c
c$$$      print *, 'OEP2INTS: AUXPH'
c$$$      call output (
c$$$     &     icore (auxph), 1, naux*n_ph, 1, 1, naux*n_ph, 1, 1)
c$$$      print *, 'OEP2INTS: SCFH'
c$$$      call output (
c$$$     &     scfham, 1, n_aa, 1, 1, n_aa, 1, 1)
c$$$      print *, 'OEP2INTS: HFXAA'
c$$$      call output (
c$$$     &     icore (hfxaa), 1, n_aa, 1, 1, n_aa, 1, 1)
c$$$      print *, 'OEP2INTS: HFXPH'
c$$$      call output (
c$$$     &     icore (hfxph), 1, n_ph, 1, 1, n_ph, 1, 1)
c
c     Calculate the reference Hamiltonian: <> Allocate memory for the
c     eigenvalues of the reference Hamiltonian. <> Check if there is
c     enough memory available. <> Calculate the reference Hamiltonian.
c     
      refevl = coretop
      coretop = refevl + iintfp * nbas
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Calculating the reference Hamiltonian.',
     &     itop-1+usdmem, maxmem+usdmem)
c     
      call oeprefham (
     &     nbas, nocc, naux, OEPH0,
     &     scfevl, scfevc, scfham,
     &     icore (hfxaa), icore (refevl))
c     
c     Calculate the inverse of the response function: <> Allocate memory
c     for the inverse and the scratch array. <> Check if there is enough
c     memory available. <> Construct the response function and calculate
c     its inverse.
c     
      resinv = coretop
      coretop = resinv + iintfp * naux * naux
      i010 = coretop
      i020 = i010 + iintfp * naux
      itop = i020 + iintfp * naux * naux
c     
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Calculating the response inverse using SVD.',
     &     Itop-1+usdmem, maxmem+usdmem)
c     
      call oepresinv (
     &     nbas, nocc, naux, OEPSVDTHR,
     &     icore (refevl), icore (auxph), 
     &     icore (i010), icore (i020), icore (resinv))
c
c     Solve the first-order OEP equation: <> Allocate memory for the
c     OEP1 solution, and the scratch arrays. <> Check if there is enough
c     memory. <> Calculate the first-order energy and write it to
c     JOBARC. <> If the local exchange potential is used, calculate the
c     first-order RHS and contract it with the response inverse to get
c     the the first-order OEP solution. If the non-local exchange
c     potential is used, set the OEPX solution to zero and add the
c     non-local exchange AOMEs to the SCF Hamiltonian. <> Store the OEPX
c     in the auxiliary representation to JOBARC.
c     
      oep1au = coretop
      coretop = oep1au + iintfp * naux
      i010 = coretop
      itop = i010 + iintfp * naux
c
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Solving the first-order OEP equation.',
     &     itop-1+usdmem, maxmem+usdmem)
c
      call getrec (20, 'JOBARC', 'NUCREP  ', iintfp, scfenr)
      call oep1enr (nbas, scfden, scfh1e, scfham, icore (hfxaa), scfenr)
      call putrec (20, 'JOBARC', 'OEP1ENER', iintfp, scfenr)
c
      call oep1rhs (
     &     nbas, nocc, naux,
     &     icore (refevl), icore (auxph), icore (hfxph), icore(i010))
      call xgemm(
     &     'n', 'n', naux, 1, naux,
     &     1.d0, icore (resinv), naux,
     &     icore (i010), naux,
     &     0.d0, icore (oep1au), naux)
c
      homo = 0.d0
      call oep1homo (
     &     nbas, nocc, naux, OEP1HYBRID,
     &     scfevc, icore (hfxaa), icore (oep1au),
     &     icore (auxaa), homo)
      call putrec (20, 'JOBARC', 'OEPXHOMO', iintfp, homo)
c
c     Solve the second-order OEP equation: <> Allocate memory. <> Check
c     if there is enough memory available. <> Check if diagonal elements
c     have to be included. <> Calculate the RHS of the second-order OEP
c     equation and the second-order energy. <> Store the second-order
c     OEP energy to JOBARC. <> Contract the response inverse with the
c     RHS to get the second-order solution. <> Store it to JOBARC. <>
c     Add the first-order solution to the second-order solution. <>
c     Store it to JOBARC. <> Calculate the AOMEs of the OEP and add them
c     to the SCF Hamiltonian.
c
c     Note that the non-local exchange MOMEs are overwritten in OEP2RHS!
c
c     Memory pointer legend: i010, i020, i030 - PP, PH, HH intemediates,
c     i040 - the second-order OEP RHS.
c
      oep2au = coretop
      coretop = oep2au + iintfp * naux
      i010 = coretop
      i020 = i010 + iintfp * n_pp
      i030 = i020 + iintfp * n_ph
      i040 = i030 + iintfp * n_hh      
      i050 = i040 + iintfp * naux
      if (dor12) then
        itop = i050 + iintfp * n_ph
      else
        itop = i050
      endif  
c
      if (itop .gt. maxmem) call insmem (
     &     'OEP2ITER-F: Solving the second-order OEP equation.',
     &     itop-1+usdmem, maxmem+usdmem)
c
      diagonal = OEPH0 (1:4) .eq. 'OEP '

      call oep2rhs (
     &     nbas, nocc, naux,
     &     OEPDFF, OEPDFW, OEPDWW,
     &     OEPFDF, OEPWDW, diagonal, OEP1HYBRID,
     &     scfenr, icore (refevl), icore (oep1au),
     &     icore (auxpp), icore (auxph), icore (auxhh),
     &     icore (hfxpp), icore (hfxph), icore (hfxhh),
     &     icore (ippph), icore (ipphh), icore (iphph), icore (iphhh),
     &     icore (i010 ), icore (i020 ), icore (i030 ), icore (i050 ), 
     &     icore (i040 ), dor12)

      call putrec (20, 'JOBARC', 'OEP2ENER', iintfp, scfenr)
c     
      call xgemm(
     &     'n', 'n', naux, 1, naux,
     &     1.d0, icore (resinv), naux,
     &     icore (i040), naux,
     &     0.d0, icore (oep2au), naux)
c     
c$$$      write (6,*) 'OEP2 RHS'
c$$$      call output (
c$$$     &     icore (i040), 1, 9, 1, 1, naux, 1, 1)
c$$$      write (6,*) 'OEP2 RESINV'
c$$$      do n=1, naux
c$$$         write (6,*)
c$$$     &        n,
c$$$     &        deref (1+naux*(n-1), icore (resinv)),
c$$$     &        deref (n, icore (i040))
c$$$      enddo
c$$$      write (6,*) 'OEP2 SOLUTION'
c$$$      call output (
c$$$     &     icore (oep2au), 1, 9, 1, 1, naux, 1, 1)
c
      call dscal (naux, 1d0-OEP1HYBRID, icore (oep1au), 1)
c     
      call putrec (
     &     20, 'JOBARC', 'OEP1AU  ', iintfp*naux, icore (oep1au))
      call putrec (
     &     20, 'JOBARC', 'OEPXAU  ', iintfp*naux, icore (oep1au))
      call putrec (
     &     20, 'JOBARC', 'OEP2AU  ', iintfp*naux, icore (oep2au))
c
      call daxpy (naux, 1.d0, icore (oep1au), 1, icore (oep2au), 1)
c     
      call putrec (
     &     20, 'JOBARC', 'OEPAU   ', iintfp*naux, icore (oep2au))
c     
      call xgemm(
     &     'n', 'n', n_aa, 1, naux,
     &     1.d0, icore (auxaa), n_aa,
     &     icore (oep2au), naux,
     &     1.d0, scfham, n_aa)
c
      call daxpy (
     &     n_aa,
     &     OEP1HYBRID, icore (hfxaa), 1,
     &     scfham, 1)
c
c$$$      call daxpy (
c$$$     &     n_aa,
c$$$     &     1.d0, icore (hfxaa), 1,
c$$$     &     scfham, 1)
c     
      return
      end
      
      function deref (n, a)
      implicit none
      integer n
      double precision a (*), deref
      deref = a (n)
      return
      end

      subroutine ao2mo (n, c, x, y, z)
      implicit none
      integer n
      double precision c (*), x (*), y (*), z (*)
      call expnd2 (x, z, n)
      call xgemm (
     &     't', 'n', n, n, n,
     &     1.d0, c, n,
     &     z, n,
     &     0.d0, y, n)
      call xgemm (
     &     'n', 'n', n, n, n,
     &     1.d0, y, n,
     &     c, n,
     &     0.d0, z, n)
      return
      end

c$$$      if (OEPH0 (1:4) .eq. 'SEMI') then
c$$$         i010 = coretop
c$$$         i020 = i010 + iintfp * n_aa
c$$$         i030 = i020 + iintfp * n_ma
c$$$         i040 = i030 + iintfp * n_pp
c$$$         i050 = i040 + iintfp * n_ph
c$$$         i060 = i050 + iintfp * n_hh
c$$$         call dcopy (n_aa, scfham, 1, icore (i010), 1)
c$$$         call daxpy (n_aa, 1.d0, icore (hfxaa), 1, icore (i010), 1)
c$$$         call oepao2mo (
c$$$     &        n_a, n_h, scfevc,
c$$$     &        icore (i010), icore (i020),
c$$$     &        icore (i030), icore (i040), icore (i050))
c$$$         print *, 'OEP2ITER: SEMI PP'
c$$$         call output (icore (i030), 1, n_pp, 1, 1, n_pp, 1, 1)
c$$$         print *, 'OEP2ITER: FOCK PH'
c$$$         call output (icore (i040), 1, n_p, 1, n_h, n_p, n_h, 1)
c$$$         print *, 'OEP2ITER: FOCK HH'
c$$$         call output (icore (i050), 1, n_hh, 1, 1, n_hh, 1, 1)
c$$$      endif
