










      subroutine oep1ints (
     &     nbas, nocc, naux, luint, lnbuf, axptype, ibuf, 
     &     buf, evec, dens, scfh, scrah,
     &     auxaa, auxph, hfxaa, hfxph)
c
c This routine takes care of all the OEP1-related integrals. It reads
c (1) AO two-electron integrals from the VMOL integral file, (2) AO MEs
c of the auxiliary potential if requested, and (3) the AO MEs of the
c auxiliary functions. Then the Coulomb AOMEs are added to the SCF
c hamiltonian matrix (supplied via SCFH), the auxiliary and non-local
c exchange AOMEs are tranformed to the MO basis.
c
c If the auxiliary potential has been requested then its AOMEs are added
c to the SCF Hamiltonian matrix as well and subtracted from the
c non-local exchange AOMEs.
c
c Note that unlike OEP2INTS, here the auxiliary integrals are read in
c after two-electron integrals, because we use AUXAA as a scratch to
c store the Coulomb AOMEs (and auxiliary potential AOMEs if needed.)
c
c In: nbas, nocc, luint, lnbuf, ibuf, buf, evec, scfh
c Scr: scrah
c Out: auxaa, auxph, hfxaa, hfxhh
c
c Igor Schweigert, Jan 2004
c $Id: oep1ints.FPP,v 1.2 2008/06/03 15:57:36 taube Exp $
c     
      implicit none
c
c     Arguments
c     
      integer
     &     nbas, nocc, naux, luint, lnbuf, ibuf (lnbuf)
c
      character*(*)
     &     axptype
c     
      double precision
     &     buf (lnbuf), evec (nbas, nbas), dens (*),  scfh (*),
     &     scrah (*), auxaa (*), auxph (*), hfxaa (*), hfxph (*)
c
c     Local variables
      logical oepint_exsist
      integer
     &     ind, nut, iand, ishft,
     &     n1, n2, n3, i1, i2, i3, i4, 
     &     n, ij, kl, ik, jl, il, jk
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
c     Read the 2e integrals in AO basis, and calculate the Coulomb and
c     non-local exchange AOMEs: <> Initialize the arrays being
c     calculated. <> Get the name of the VMOL file from JOBARC, open the
c     file, and shift the pointer to the begining of the 2e integral
c     record. <> Read the integrals and calculate the AOMEs of the
c     Coulomb and non-local exchange operators (Coulomb AOMEs will be
c     stored in SCRMA): <<>> Read the value and indices of the current
c     integral. <<>> Define the symmetry of the integral. <<>> Update
c     the corresponding integrals. <> Close the VMOL file. <> Scale the
c     diagonal elements by factor of two. <> Add the Coulomb AOMEs to
c     the SCF Hamiltonian. 
c
      asv_null=achar(0)
      call zero (auxaa, n_aa)
      call zero (hfxaa, n_aa)
      call zero (hfxph, n_ph)
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
            x = buf(n1)
            i = iupki ( ibuf(n1))
            j = iupkj ( ibuf(n1))
            k = iupkk ( ibuf(n1))
            l = iupkl ( ibuf(n1))
c     
            ij = i_aa_cr (j, i)
            kl = i_aa_cr (l, k)
            ik = i_aa_cr (k, i)
            jl = i_aa_cr (l, j)
            il = i_aa_cr (l, i)
            jk = i_aa_cr (j, k)
c     
            if (i.eq.j) x = .5d0 * x
            if (k.eq.l) x = .5d0 * x
            if (ij.eq.kl) x = .5d0 * x
c     
            auxaa (ij) = auxaa (ij) + 2.d0 * dens (kl) * x
            auxaa (kl) = auxaa (kl) + 2.d0 * dens (ij) * x
c     
            hfxaa (ik) = hfxaa (ik) - 0.5d0 * dens (jl) * x
            hfxaa (jl) = hfxaa (jl) - 0.5d0 * dens (ik) * x
            hfxaa (il) = hfxaa (il) - 0.5d0 * dens (jk) * x
            hfxaa (jk) = hfxaa (jk) - 0.5d0 * dens (il) * x
         enddo
      enddo
c     
      close(unit=luint,status='KEEP')
c     
      do i = 1, n_a
         n = i_aa (i, i)
         auxaa (n) = 2.d0 * auxaa (n)
         hfxaa (n) = 2.d0 * hfxaa (n)
      enddo
c 
      call daxpy (n_aa, 1.d0, auxaa, 1, scfh, 1)
c
c     If requested, obtain the AOMEs of the auxiliary potential used to
c     correct the asymptotic behavior of the exchange OEP. <> If the
c     scaled Coulomb is used as the auxiliary potential, scale the
c     Coulomb AOMEs, update the SCF Hamiltonian, and non-local exchange
c     AOMEs. <> If the Slater potential is used as the auxiliary
c     potential, get the Slater AOMEs from the OEP integral file, update
c     the SCF Hamiltonian and non-local exchange AOMEs.
c
c
      if     (axptype (1:3) .eq. 'FA'//asv_null) then
c
         call dscal (n_aa, 1.d0/dble (2*nocc), auxaa, 1)
         call daxpy (n_aa,-1.d0, auxaa, 1, scfh, 1)
         call daxpy (n_aa, 1.d0, auxaa, 1, hfxaa, 1)
c     
      elseif (axptype (1:7) .eq. 'SLATER'//asv_null) then
       call getrec(1,'JOBARC','SLAT51AO',iintfp*n_aa,auxaa)

c         open(
c     &        unit = luint, file = 'SLAT51AB',
c     &        form = 'UNFORMATTED', access = 'SEQUENTIAL')
c         call locate (luint, 'SLAT51AB')
c     
c         nut = lnbuf
c         do while (nut.eq.lnbuf)
c            read(luint) buf, ibuf, nut
c            do ind = 1, nut
c               auxaa (ibuf (ind)) = buf(ind)
c            enddo
c         enddo
c     
         close (luint)
c
         call daxpy (n_aa, 1.d0, auxaa, 1, scfh, 1)
         call daxpy (n_aa,-1.d0, auxaa, 1, hfxaa, 1)
c
      endif
c
c     Calculate the PH part of non-local exchange MOME: <> Transform the
c     second AO index into a hole index. <> Transform the first AO index
c     into a particle index.
c     
      call zero (scrah, n_ah)
      do n1 = 1, n_ah
         i2 = i2_ah (n1)
         i1 = i1_ah (n1, i2)
         do i3 = 1, n_a
            n2 = i_aa_cr (i1, i3)
            scrah (n1) = scrah (n1) + hfxaa (n2) * evec (i3, i2)
         enddo
      enddo
c     
      do n1 = 1, n_ph
         i2 = i2_ph (n1)
         i1 = i1_ph (n1, i2) + n_h
         do i3 = 1, n_a
            n2 = i_ah (i3, i2)
            hfxph (n1) = hfxph (n1) + scrah (n2) * evec (i3, i1)
         enddo
      enddo
c
c     Read the auxiliary AOMEs from the OEPINT file and transform them
c     to MOs: <> Initialize the array being calculated.<> Open the
c     OEPINT file and locate the auxiliary integrals record. <> Read the
c     integrals. <> Close the OEPINT file. <> For each auxiliary
c     function, transform the AOMEs into PH part of the MOMEs.
c
      call zero (auxaa, n_aa * naux)
      call zero (auxph, n_ph * naux)
      inquire (file = 'OEPINT', exist = oepint_exsist)
         if (.not. oepint_exsist) then
            write (6,*)    'crappppp'
            call errex
         endif    
      open(
     &     unit = luint, file = 'OEPINT',
     &     form = 'UNFORMATTED', access = 'SEQUENTIAL')
      call locate (luint, 'AUX3CNTR')
c     
      nut = lnbuf
      do while (nut.eq.lnbuf)
         read(luint) buf, ibuf, nut
         do ind = 1, nut
            auxaa (ibuf (ind)) = buf(ind)
         enddo
      enddo
c     
      close (luint)
c  
      do i4 = 1, naux
         call zero (scrah, n_ah)
         do n1 = 1, n_ah
            i2 = i2_ah (n1)
            i1 = i1_ah (n1, i2)
            do i3 = 1, n_a
               n2 = i_aa_cr (i1, i3) + n_aa * (i4-1)
               scrah (n1) = scrah (n1) + auxaa (n2) * evec (i3, i2)
            enddo
         enddo
         do n1 = 1, n_ph
            i2 = i2_ph (n1)
            i1 = i1_ph (n1, i2) + n_h
            n3 = n1 + n_ph * (i4-1)
            do i3 = 1, n_a
               n2 = i_ah (i3, i2)
               auxph (n3) = auxph (n3) + scrah (n2) * evec (i3, i1)
            enddo
         enddo
      enddo
c     
      return
      end
