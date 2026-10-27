





























































































































































































































      Subroutine Erg_dens(Iuhf,Lcorr,N,B,Ha,Hb,Xa,Xb,Dt,W,C,Occ,Eps,
     &                    Scr)
C
C Densities over ERD functions, D = B D(SO) B(T) (Erg_bso). Per spin
C
C    H   = C n C(T)        the SCF (reference) density, n = ORB_OCC*
C
C and for an SCF gradient W = sum_spin C n e C(T) (e = SCFEVL*). For
C a correlated gradient (Lcorr) the relaxed density R (RELDENS*, MO
C basis) and the intermediates I (IINTERM*) give, as in xvdint's
C Onedso and Dsofso,
C
C    X   = C R C(T) - H    the correlation part of the density
C    W   = - C sym(I) C(T)
C
C Once the correlated programs have run (Lcord: a correlated or a
C second derivative calculation), SCFEVC* and SCFEVL* are in the
C correlated MO order, all occupied orbitals first, like R
C and I (xvdint's Symcext and Symd turn all of them back to the SCF
C order instead). ORB_OCC* stays in the SCF order, so in that case
C the occupied orbitals are the first ones, as many as OCCUPY*0 says
C (SYMPOPO* would leave out dropped core orbitals, which are still
C there).
C For RHF, RELDENSA and IINTERMA are the totals (both spins): Xa = Xb
C is half of C R C(T) - 2 H. Ha, Hb, Xa, Xb are per spin, Dt is the
C total one-particle density (H + X over both spins), W the total
C energy weighted density.
C
      Implicit Double Precision (A-H,O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c erdints.com : begin
c
c Basis set, geometry and symmetry data for xerdints (integrals from
c Flocke's OED/ERD, symmetry adapted to the ACES II SO basis without
c VMOL or VMOL2JA). Filled by Erd_rdmol and Erd_symgrp.
c
c Atoms are those of MOL_ERDOED (all real centers, ZMAT order).
c   NATM          : number of atoms
c   ZATM(i)       : nuclear charge
c   XYZATM(3,i)   : coordinates (bohr, computational orientation)
c
c A shell is one angular momentum block on one atom, possibly with
c several (general) contractions, as ERD wants it.
c   NSHL          : number of shells
c   ISHATM(s)     : atom of shell s
c   LSHL(s)       : angular momentum
c   NPRSH(s)      : number of primitives
c   NCTSH(s)      : number of contractions
c   IOFALP(s)     : ALPHA(IOFALP(s)+1 ...) are the exponents
c   IOFCC(s)      : CCOEF(IOFCC(s)+1 ...) are the coefficients,
c                   NPRSH x NCTSH, column major
c   IOFCB(s)      : ICCBEG/ICCEND(IOFCB(s)+1 ...) are the first/last
c                   nonzero primitive of each contraction
c   NFNSH(s)      : number of functions (components x contractions)
c   IOFFN(s)      : offset of the first function of s in ERD order;
c                   within a shell the components run fastest
c                   (contraction major)
c   IOFCRT(s)     : offset of the first Cartesian function of s in
c                   the ZMAT AO order (= ERD Cartesian order)
c   NFNERD        : total number of functions in ERD order
c   NCART         : total number of Cartesian functions
c   SPHER         : spherical harmonics (ERD in spherical mode)
c
c Symmetry (D2h and subgroups; every operation is a sign triple).
c   NGEN          : number of generators (0-3)
c   NOPS          : number of operations, 2**NGEN
c   IOPSGN(k,r)   : sign of coordinate k under operation r (r=0 is E);
c                   bit g of r set means generator g+1 is included
c   IATMAP(i,r)   : image of atom i under operation r
c   ISHMAP(s,r)   : image of shell s under operation r
c
      integer mxshl, mxprmt, mxcct, mxctrt
      parameter (mxshl=5000, mxprmt=50000, mxcct=500000,
     &           mxctrt=20000)

      integer natm, nshl, ishatm, lshl, nprsh, nctsh, iofalp, iofcc,
     &        iofcb, nfnsh, ioffn, iofcrt, nfnerd, ncart, iccbeg,
     &        iccend
      double precision zatm, xyzatm, alpha, ccoef
      logical spher
      common /erdbas/ zatm(mxatms), xyzatm(3,mxatms),
     &                alpha(mxprmt), ccoef(mxcct),
     &                natm, nshl, ishatm(mxshl), lshl(mxshl),
     &                nprsh(mxshl), nctsh(mxshl), iofalp(mxshl),
     &                iofcc(mxshl), iofcb(mxshl), nfnsh(mxshl),
     &                ioffn(mxshl), iofcrt(mxshl), nfnerd, ncart,
     &                iccbeg(mxctrt), iccend(mxctrt), spher

c   CCRAW(i)      : the coefficients as read (MOL_ERDOED, i.e. the MOL
c                   file's), same layout as CCOEF; VMOL writes these
c                   for vprops (Erd_wrvp)
      double precision ccraw
      common /erdraw/ ccraw(mxcct)

      integer ngen, nops, iopsgn, iatmap, ishmap
      common /erdsym/ ngen, nops, iopsgn(3,0:7), iatmap(mxatms,0:7),
     &                ishmap(mxshl,0:7)

c SO basis (Erd_sobld). A shell orbit is a shell and its images; each
c SO is a +-1 combination of one function (same component and
c contraction) on every member shell of an orbit, i.e. VMOL's
c unnormalized SOs. SOs are numbered irrep by irrep.
c   NIRR          : number of irreps (= NOPS)
c   NSO           : number of SOs (= NFNERD)
c   NBFIRR(g)     : SOs in irrep g (1 based, ACES irrep g)
c   NORB          : number of shell orbits
c   NORBM(k)      : members of orbit k (orbit size, 1 2 4 or 8)
c   IORBSH(m,k)   : member shells, IORBSH(1,k) is the representative
c   IORBOP(m,k)   : an operation taking IORBSH(1,k) to IORBSH(m,k)
c   KORBSH(s)     : orbit of shell s
c   IORBSO(k)     : LORBSO(IORBSO(k)+1 ...) are the SOs of orbit k,
c                   NORBM(k)*NFNSH(rep) of them
c   ISOIRR(p)     : irrep of SO p (1 based)
c   ISOORB(p)     : orbit of SO p
c   ISOFN(p)      : function (component x contraction, ERD order)
c                   within the orbit's shells
c   ISOSGN(m,p)   : sign of member m in SO p
c   FNFAC(i)      : ERD function i times FNFAC(i) is the ACES (VMOL)
c                   AO; for Cartesian components x^a y^b z^c this is
c                   sqrt((2a-1)!!(2b-1)!!(2c-1)!!), 1 for spherical ones
c   WVSPH(m,j,L)  : spherical basis, L >= 2: VMOL harmonic j =
c                   sum_m WVSPH(m,j,L) ERD spherical function m
c                   (Erd_vsph). The basis functions SOs are made of are
c                   VMOL's: Erd_fntr turns ERD integrals into them.
      integer mxfn, mxorb
      parameter (mxfn=20000, mxorb=mxshl)
      integer nirr, nso, nbfirr, norb, norbm, iorbsh, iorbop, korbsh,
     &        iorbso, lorbso, isoirr, isoorb, isofn, isosgn
      common /erdso/ nirr, nso, nbfirr(8), norb, norbm(mxorb),
     &               iorbsh(8,mxorb), iorbop(8,mxorb), korbsh(mxshl),
     &               iorbso(mxorb), lorbso(mxfn), isoirr(mxfn),
     &               isoorb(mxfn), isofn(mxfn), isosgn(8,mxfn)
      double precision fnfac, wvsph
      common /erdsof/ fnfac(mxfn), wvsph(13,13,2:6)

c erdints.com : end


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



c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      Logical Lcorr, Lcord
      Dimension B(N,N), Ha(N,N), Hb(N,N), Xa(N,N), Xb(N,N), Dt(N,N),
     &          W(N,N), C(N,N), Occ(N), Eps(N), Scr(N,N,2)
      Integer Iocc(8,2)
      Character*8 Labc(2), Labe(2), Labo(2), Labd(2), Labi(2)
      Data Labc /'SCFEVCA0','SCFEVCB0'/
      Data Labe /'SCFEVLA0','SCFEVLB0'/
      Data Labo /'ORB_OCCA','ORB_OCCB'/
      Data Labd /'RELDENSA','RELDENSB'/
      Data Labi /'IINTERMA','IINTERMB'/

      Lcord = Lcorr .Or. Iflags(3) .Eq. 2
      If (Lcord) Then
         Call Getrec(20,'JOBARC','OCCUPYA0',Nirr,Iocc(1,1))
         If (Iuhf .Eq. 1) Then
            Call Getrec(20,'JOBARC','OCCUPYB0',Nirr,Iocc(1,2))
         Endif
      Endif
      Call Zero(W,N*N)
      Call Zero(Xa,N*N)
      Call Zero(Xb,N*N)
      Fac = 1.0D0
      If (Iuhf .Eq. 0) Fac = 2.0D0
      Do Ispin = 1, Iuhf+1
         Call Getrec(20,'JOBARC',Labc(Ispin),N*N*Iintfp,C)
         Call Getrec(20,'JOBARC',Labe(Ispin),N*Iintfp,Eps)
         If (Lcord) Then
            Nocc = 0
            Do Ir = 1, Nirr
               Nocc = Nocc + Iocc(Ir,Ispin)
            Enddo
            Do J = 1, N
               Occ(J) = 0.0D0
               If (J .Le. Nocc) Occ(J) = 1.0D0
            Enddo
         Else
            Call Getrec(20,'JOBARC',Labo(Ispin),N*Iintfp,Occ)
         Endif

C H = B C n C(T) B(T).

         Do J = 1, N
            Do I = 1, N
               Scr(I,J,2) = Occ(J)*C(I,J)
            Enddo
         Enddo
         Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,C,N,0.0D0,
     &              Scr(1,1,1),N)
         If (Ispin .Eq. 1) Then
            Call Erg_btr(N,B,Scr(1,1,1),Ha,Scr(1,1,2),0.0D0)
         Else
            Call Erg_btr(N,B,Scr(1,1,1),Hb,Scr(1,1,2),0.0D0)
         Endif

         If (.Not. Lcorr) Then

C W += B C (fac n e) C(T) B(T).

            Do J = 1, N
               Do I = 1, N
                  Scr(I,J,2) = Fac*Occ(J)*Eps(J)*C(I,J)
               Enddo
            Enddo
            Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,C,N,0.0D0,
     &                 Scr(1,1,1),N)
            Call Erg_btr(N,B,Scr(1,1,1),W,Scr(1,1,2),1.0D0)
         Else

C X = B C R C(T) B(T) / fac - H, W -= B C sym(I) C(T) B(T).

            Do Ipass = 1, 2
               If (Ipass .Eq. 1) Then
                  Call Getrec(20,'JOBARC',Labd(Ispin),N*N*Iintfp,
     &                        Scr(1,1,1))
               Else
                  Call Getrec(20,'JOBARC',Labi(Ispin),N*N*Iintfp,
     &                        Scr(1,1,1))
                  Do J = 1, N
                     Do I = 1, J-1
                        X = 0.5D0*(Scr(I,J,1) + Scr(J,I,1))
                        Scr(I,J,1) = X
                        Scr(J,I,1) = X
                     Enddo
                  Enddo
               Endif
               Call Dgemm('N','N',N,N,N,1.0D0,C,N,Scr(1,1,1),N,
     &                    0.0D0,Scr(1,1,2),N)
               Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,C,N,
     &                    0.0D0,Scr(1,1,1),N)
               If (Ipass .Eq. 1) Then
                  Call Dscal(N*N,1.0D0/Fac,Scr(1,1,1),1)
                  If (Ispin .Eq. 1) Then
                     Call Erg_btr(N,B,Scr(1,1,1),Xa,Scr(1,1,2),0.0D0)
                     Call Daxpy(N*N,-1.0D0,Ha,1,Xa,1)
                  Else
                     Call Erg_btr(N,B,Scr(1,1,1),Xb,Scr(1,1,2),0.0D0)
                     Call Daxpy(N*N,-1.0D0,Hb,1,Xb,1)
                  Endif
               Else
                  Call Dscal(N*N,-1.0D0,Scr(1,1,1),1)
                  Call Erg_btr(N,B,Scr(1,1,1),W,Scr(1,1,2),1.0D0)
               Endif
            Enddo
         Endif
      Enddo
      If (Iuhf .Eq. 0) Then
         Call Dcopy(N*N,Ha,1,Hb,1)
         Call Dcopy(N*N,Xa,1,Xb,1)
      Endif

      Call Dcopy(N*N,Ha,1,Dt,1)
      Call Daxpy(N*N,1.0D0,Hb,1,Dt,1)
      Call Daxpy(N*N,1.0D0,Xa,1,Dt,1)
      Call Daxpy(N*N,1.0D0,Xb,1,Dt,1)

      Return
      End
