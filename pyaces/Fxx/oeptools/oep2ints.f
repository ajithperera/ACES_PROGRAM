










      subroutine oep2ints (
     &     isaux, nbas, nocc, naux, luint, lnbuf, axptype, ibuf, 
     &     buf, evec, scfh, 
     &     auxaa, auxpp, auxph, auxhh,
     &     hfxaa, hfxpp, hfxph, hfxhh,
     &     paah, aaah, ppah, phah,
     &     ppph, pphh, phph, phhh)
c
c This routine takes care of all the OEP2-related integrals. It reads
c the AO two-electron integrals from the VMOL integral file and
c trasforms them to MO basis. As the result (a) AOMEs of Coulomb
c potential and non-local exchange operator (b) MOMEs of the non-local
c exchange operator (c) the MO two-electron integrals (if requested) are
c constructed.
c
c If ISAUX is false than than the auxiliary integrals are not read.
c
c The Coulomb AOMEs are added to the SCF hamiltonian matrix (supplied
c via SCFH.)
c
c For the details on the AO to MO transformation of the two-electron
c integrals, see OEPLIB/OEPINTTRN.F (.TEX).
c
c In:  isaux, nbas, nocc, luint, lnbuf, auxtype, ibuf, buf, evec, scfh
c Scr: paah, aaah, ppah, phah
c Out: auxaa, auxpp, auxph, auxhh, 
c      hfxaa, hfxpp, hfxph, hfxhh, 
c      ppph, pphh, phph, phhh
c
c Igor Schweigert, Jan 2004
c $Id: oep2ints.FPP,v 1.2 2008/06/03 15:57:36 taube Exp $
c     
      implicit none
c
c     Arguments
c
      logical isaux
c     
      integer
     &     nbas, nocc, naux, luint, lnbuf, ibuf (lnbuf)
c
      character*(*) axptype
c     
      double precision
     &     buf (lnbuf), evec (nbas, nbas), scfh (*), 
     &     auxaa (*), auxpp (*), auxph (*), auxhh (*),
     &     hfxaa (*), hfxpp (*), hfxph (*), hfxhh (*),
     &     paah (*), aaah (*), ppah (*), phah (*),
     &     ppph (*), pphh (*), phph (*), phhh (*)
c
c     Local variables
c
      logical
     &     i1i2, i3i4, i1i2i3i4
c           
      integer
     &     ind, nut, iand, ishft,
     &     n1, n2, i1, i2, i3, i4, i5
c     
      character*80
     &     fname
      character*1 asv_null
c
      double precision x
c
c     Define the two-electron integral offset operators and the
c     associated variables.
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
c     Read the auxiliary AOMEs from the OEPINT file and transform one AO
c     index into a MO index: <> If the auxiliary integrals are not
c     requested, skip this part. <> Initialize the array being
c     calculated (the PPPH array will be used as a scratch). <> Open the
c     OEPINT file and locate the auxiliary integrals record. <> Read the
c     integrals and transform the first AO index into a MO index. <>
c     Close the OEPINT file. <> Transform the second AO index into a MO
c     index and store the resulting auxiliary MOMEs in to the
c     corresponding arrays.
c
      asv_null=achar(0)
      if (.not.isaux) goto 100
c     
      call zero (auxaa, n_aa * naux)
      call zero (auxpp, n_pp * naux)
      call zero (auxph, n_ph * naux)
      call zero (auxhh, n_hh * naux)
c     
      open(
     &     unit = luint, file = 'OEPINT',
     &     form = 'UNFORMATTED', access = 'SEQUENTIAL')
      call locate (luint, 'AUX3CNTR')
c     
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read (luint) buf, ibuf, nut
         do ind = 1, nut
            auxaa (ibuf (ind)) = buf(ind)
         enddo
      enddo
c     
      close (luint)
c     
      do n1 = 1, naux
         call oepao2mo (
     &        nbas, nocc, evec,
     &        auxaa (n_aa*(n1-1)+1), paah,
     &        auxpp (n_pp*(n1-1)+1), auxph (n_ph*(n1-1)+1),
     &        auxhh (n_hh*(n1-1)+1))
      enddo
c
 100  continue
c
c     Read the 2e integrals in AO basis, and transform one AO index into
c     a hole index: <> Initialize the arrays being calculated. <> Get
c     the name of the VMOL file from JOBARC, open the file, and shift
c     the pointer to the begining of the 2e integral record. <> Read the
c     integrals and transform one AO index to a hole index: <<>> Read
c     the value and indices of the current integral. <<>> Define the
c     symmetry of the integral. <<>> Loop over the hole index, check the
c     symmetry of the integral and update the corresponding integrals in
c     the AAAH array. <> Close the VMOL file.
c
      call zero (aaah, n_aaah)
c     
      call gfname('IIII    ',fname,n1)
      open(unit=luint,file=fname(1:n1),form='UNFORMATTED',
     &     access='SEQUENTIAL')
      call locate(luint,'TWOELSUP')
c     
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do n1=1,nut
c
            x=buf(n1)
            i1=iupki(ibuf(n1))
            i2=iupkj(ibuf(n1))
            i3=iupkk(ibuf(n1))
            i4=iupkl(ibuf(n1))
c
            i1i2 = i1.ne.i2
            i3i4 = i3.ne.i4
            i1i2i3i4 = i_aa (i1,i2) .ne. i_aa (i3,i4)
c     
            do i5=1, nocc
               aaah (i_aaah (i1, i2, i3, i5)) =
     &              aaah (i_aaah (i1, i2, i3, i5)) + evec (i4,i5)*x
               if (i3i4)
     &              aaah (i_aaah (i1, i2, i4, i5)) =
     &              aaah (i_aaah (i1, i2, i4, i5)) + evec (i3,i5)*x
               if (i1i2i3i4)
     &              aaah (i_aaah (i3, i4, i1, i5)) =
     &              aaah (i_aaah (i3, i4, i1, i5)) + evec (i2,i5)*x
               if (i1i2i3i4.and.i1i2)
     &              aaah (i_aaah (i3, i4, i2, i5)) =
     &              aaah (i_aaah (i3, i4, i2, i5)) + evec (i1,i5)*x
            enddo
         enddo
      enddo
c     
      close(unit=luint,status='KEEP')
c
c     Calculate the AOMEs of the Coulomb and non-local exchange
c     operators: <> Initialize the scratch and target arrays to
c     zeros. <> Calculate the Coulomb and non-local exchange
c     AOMEs. Since we might need the Coulomb AOMEs to correct
c     asymptotic, store them temporarily in PAAH. <> Add the Coulomb
c     AOMEs to the SCF Hamiltonian with the Coulomb. <> If requested,
c     obtain the AOMEs of the auxiliary exchange potential, add them to
c     the SCF Hamiltonian and subtract them from the non-local exchange
c     AOMEs. <> Transform non-local exchange AOMEs to the MOs.
c     
      call zero (hfxaa, n_aa)
      call zero (paah,  n_aa)
      do n1 = 1, n_aaah
         i4 = i4_aaah (n1)
         i3 = i3_aaah (n1, i4)
         i2 = i2_aaah (n1, i3, i4)
         i1 = i1_aaah (n1, i2, i3, i4)
         n2 = i_aa (i1, i2)
         paah (n2) = paah (n2) + 2.d0 * aaah (n1) * evec (i3, i4) 
         if (i1.le.i3) then
            n2 = i_aa (i1, i3)
            hfxaa (n2) = hfxaa (n2) - aaah (n1) * evec (i2, i4)
         endif
         if (i1.ne.i2 .and. i2.le.i3) then
            n2 = i_aa (i2, i3)
            hfxaa (n2) = hfxaa (n2) - aaah (n1) * evec (i1, i4)
         endif
      enddo
c
      call daxpy (n_aa, 1.d0, paah, 1, scfh, 1)
c     
      if     (axptype (1:3) .eq. 'FA'//asv_null) then
         call dscal (n_aa,-1.d0/dble (2*nocc), paah, 1)
         call daxpy (n_aa, 1.d0, paah, 1, scfh, 1)
         call daxpy (n_aa,-1.d0, paah, 1, hfxaa, 1)
      elseif (axptype (1:7) .eq. 'SLATER'//asv_null) then
          call getrec(1,'JOBARC','SLAT51AO',iintfp*n_aa,paah)
c         open(
c     &        unit = luint, file = 'OEPINT',
c     &        form = 'UNFORMATTED', access = 'SEQUENTIAL')
c         call locate (luint, 'SLAT51AO')
c         nut = lnbuf
c         do while (nut.eq.lnbuf)
c            read(luint) buf, ibuf, nut
c            do ind = 1, nut
c               paah (ibuf (ind)) = buf(ind)
c            enddo
c         enddo
c         close (luint)
         call daxpy (n_aa, 1.d0, paah, 1, scfh, 1)
         call daxpy (n_aa,-1.d0, paah, 1, hfxaa, 1)
      endif
c
      call zero (hfxpp, n_pp)
      call zero (hfxph, n_ph)
      call zero (hfxhh, n_hh)
      call zero (paah, n_ma)
      call oepao2mo (
     &     nbas, nocc, evec,
     &     hfxaa, paah, hfxpp, hfxph, hfxhh)
c
c     AAAH->PAAH. Convert the first AO index to a particle index: <>
c     Initialize the target array. <> Loop over the new particle index
c     and remaining indices, and update the target array.
c     
      call zero (paah, n_paah)
c     
      do n1 = 1, n_aaah
         i4 = i4_aaah (n1)
         i3 = i3_aaah (n1, i4)
         i2 = i2_aaah (n1, i3, i4)
         i1 = i1_aaah (n1, i2, i3, i4)
         i1i2 = i1.ne.i2
         do i5 = 1, nbas-nocc
            n2 = i_paah (i5, i2, i3, i4)
            paah (n2) = paah (n2) + aaah (n1) * evec (i1, i5+nocc)
            if (i1.ne.i2) then
               n2 = i_paah (i5, i1, i3, i4)
               paah (n2) = paah (n2) + aaah (n1) * evec (i2, i5+nocc)
            endif
         enddo
      enddo
c
c     PAAH->PPAH,PHAH. Convert the second AO index to a particle and a
c     hole indices: <> Initialize the target arrays. <> Loop over the
c     source array and construct the target arrays.
c     
      call zero (ppah, n_ppah)
      call zero (phah, n_phah)
c     
      do n1 = 1, n_paah
         i4 = i4_paah (n1)
         i3 = i3_paah (n1, i4)
         i2 = i2_paah (n1, i3, i4)
         i1 = i1_paah (n1, i2, i3, i4)
         do i5 = i1, nbas-nocc
            n2 = i_ppah (i1, i5, i3, i4)
            ppah (n2) = ppah (n2) + paah (n1) * evec (i2, i5+nocc)
         enddo
         do i5 = 1, nocc
            n2 = i_phah (i1, i5, i3, i4)
            phah (n2) = phah (n2) + paah (n1) * evec (i2, i5)
         enddo
      enddo
c
c     PPAH->PPPH, PPHH. Convert the third AO index to a particle and a
c     hole indices: <> Initialize the target arrays. <> Loop over the
c     source array and construct the target arrays.
c     
      call zero (ppph, n_ppph)
      call zero (pphh, n_pphh)
c     
      do n1 = 1, n_ppah
         i4 = i4_ppah (n1)
         i3 = i3_ppah (n1, i4)
         i2 = i2_ppah (n1, i3, i4)
         i1 = i1_ppah (n1, i2, i3, i4)
         do i5 = 1, nbas-nocc
            n2 = i_ppph (i1, i2, i5, i4)
            ppph (n2) = ppph (n2) + ppah (n1) * evec (i3, i5+nocc)
         enddo
         do i5 = 1, i4
            n2 = i_pphh (i1, i2, i5, i4)
            pphh (n2) = pphh (n2) + ppah (n1) * evec (i3, i5)
         enddo
      enddo
c
c     PHAH->PHPH, PHHH. Convert the third AO index.: <> Initialize the
c     target arrays. <> Loop over the new indices and remaining indices,
c     and update the target arrays.
c     
      call zero (phph, n_phph)
      call zero (phhh, n_phhh)
c     
      do n1 = 1, n_phah
         i4 = i4_phah (n1)
         i3 = i3_phah (n1, i4)
         i2 = i2_phah (n1, i3, i4)
         i1 = i1_phah (n1, i2, i3, i4)
         do i5 = 1, nbas-nocc
            if (i_ph (i1, i2) .le. i_ph (i5, i4)) then
               n2 = i_phph (i1, i2, i5, i4)
               phph (n2) = phph (n2) + phah (n1) * evec (i3, i5+nocc)
            endif
         enddo
         do i5 = 1, i4
            n2 = i_phhh (i1, i2, i5, i4)
            phhh (n2) = phhh (n2) + phah (n1) * evec (i3, i5)
         enddo
      enddo
c     
      return
      end

