










      Subroutine Erd_sobld
C
C Build the SO basis from the shell orbits (fills the erdso common).
C
C For a function mu on the representative shell A of an orbit, the
C projector onto irrep G gives
C
C    P_G mu_A = sum_r chi_G(r) O_r mu_A = sum_r chi_G(r) p_mu(r) mu_r(A)
C
C where p_mu(r) = +-1 is the parity of mu under r (Erd_cmppar) and
C chi_G(r) = (-1)**popcount(iand(G,r)). It vanishes unless
C chi_G(h) = p_mu(h) for every h of the stabilizer H of A, which
C leaves exactly |orbit| irreps. For those, the coefficient of the
C member B = r_B(A) is |H| chi_G(r_B) p_mu(r_B); dropping |H| gives
C VMOL's unnormalized SOs, sum over members of +-mu_B.
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

      Integer Ipar(3,64), Ipsh(0:7,64), Ichr(0:7,0:7), Nfill(8),
     &        Iso0(8)
      Dimension Fac(64), Cv(36,36)
      Integer Iparv(3,36)
      Logical Compat

      Nirr = Nops
      Nso  = Nfnerd
      If (Nso .Gt. Mxfn) Then
         Write(6,"(a,i6)") " @Erd_sobld: too many functions, max ",
     &                     Mxfn
         Call Errex
      Endif

C Characters.

      Do Ig = 0, Nirr-1
         Do Ir = 0, Nops-1
            Ichr(Ir,Ig) = 1
            If (Mod(Popcnt(Iand(Ig,Ir)),2) .Eq. 1) Ichr(Ir,Ig) = -1
         Enddo
      Enddo

C Shell orbits; members in the order of the operation generating
C them first.

      Do Is = 1, Nshl
         Korbsh(Is) = 0
      Enddo
      Norb = 0
      Do Is = 1, Nshl
         If (Korbsh(Is) .Eq. 0) Then
            Norb = Norb + 1
            Norbm(Norb) = 0
            Do Ir = 0, Nops-1
               Js = Ishmap(Is,Ir)
               If (Korbsh(Js) .Eq. 0) Then
                  Norbm(Norb) = Norbm(Norb) + 1
                  Iorbsh(Norbm(Norb),Norb) = Js
                  Iorbop(Norbm(Norb),Norb) = Ir
                  Korbsh(Js) = Norb
               Endif
            Enddo
         Endif
      Enddo

C First pass: count the SOs of each irrep. Second pass: assign them,
C irrep by irrep, in orbit order and then function order.

      Do Ipass = 1, 2
         Do Ig = 1, Nirr
            Nfill(Ig) = 0
         Enddo
         If (Ipass .Eq. 2) Then
            Iso0(1) = 0
            Do Ig = 2, Nirr
               Iso0(Ig) = Iso0(Ig-1) + Nbfirr(Ig-1)
            Enddo
         Endif
         Nlist = 0

         Do Iorb = 1, Norb
            Irep = Iorbsh(1,Iorb)
            Nm   = Norbm(Iorb)
            L    = Lshl(Irep)
            Nct  = Nctsh(Irep)
            Call Erd_cmppar(L,Ipar,Ncmp)
            Do Ic = 1, Ncmp
               Do Ir = 0, Nops-1
                  Ipsh(Ir,Ic) = 1
                  Do K = 1, 3
                     If (Ipar(K,Ic) .Eq. 1)
     &                  Ipsh(Ir,Ic) = Ipsh(Ir,Ic)*Iopsgn(K,Ir)
                  Enddo
               Enddo
            Enddo
            If (Ipass .Eq. 2) Iorbso(Iorb) = Nlist

C Within a shell ERD runs the components fastest (contraction major).

            Do If = 1, Nfnsh(Irep)
               Ic = Mod(If-1,Ncmp) + 1
               Ncomp = 0
               Do Ig = 0, Nirr-1
                  Compat = .True.
                  Do Ir = 0, Nops-1
                     If (Ishmap(Irep,Ir) .Eq. Irep .And.
     &                   Ichr(Ir,Ig) .Ne. Ipsh(Ir,Ic)) Compat = .False.
                  Enddo
                  If (Compat) Then
                     Ncomp = Ncomp + 1
                     Nfill(Ig+1) = Nfill(Ig+1) + 1
                     If (Ipass .Eq. 2) Then
                        Ip = Iso0(Ig+1) + Nfill(Ig+1)
                        Isoirr(Ip) = Ig + 1
                        Isoorb(Ip) = Iorb
                        Isofn(Ip)  = If
                        Do M = 1, Nm
                           Ir = Iorbop(M,Iorb)
                           Isosgn(M,Ip) = Ichr(Ir,Ig)*Ipsh(Ir,Ic)
                        Enddo
                        Nlist = Nlist + 1
                        Lorbso(Nlist) = Ip
                     Endif
                  Endif
               Enddo
               If (Ncomp .Ne. Nm) Then
                  Write(6,"(a,3i5)") " @Erd_sobld: irrep count error ",
     &                               Iorb, If, Ncomp
                  Call Errex
               Endif
            Enddo
         Enddo

         If (Ipass .Eq. 1) Then
            Do Ig = 1, Nirr
               Nbfirr(Ig) = Nfill(Ig)
            Enddo
         Endif
      Enddo

C Scale factors between ERD's Cartesian functions and VMOL's, and for
C a spherical basis the ERD -> VMOL harmonic transformations.

      If (Spher) Then
         Lmx = 0
         Do Is = 1, Nshl
            Lmx = Max(Lmx,Lshl(Is))
         Enddo
         Do L = 2, Lmx
            Call Erd_vsph(L,Cv,36,Wvsph(1,1,L),13,Iparv)
         Enddo
      Endif

      Do Is = 1, Nshl
         L = Lshl(Is)
         Call Erd_cmpfac(L,Spher,Fac,Ncmp)
         Do If = 1, Nfnsh(Is)
            Fnfac(Ioffn(Is)+If) = Fac(Mod(If-1,Ncmp)+1)
         Enddo
      Enddo

      If (Nlist .Ne. Nso) Then
         Write(6,"(a,2i6)") " @Erd_sobld: SO count error ", Nlist, Nso
         Call Errex
      Endif

      Return
      End
