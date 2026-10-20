










      Subroutine Pccd_rfo_step(Work,Maxcor,Nbas,Nocc,Iuhf,Grdcnv,Tol,
     +                         Etol,Gmax,Grms,Ncycle)
C
C OO_METHOD=RFO: trust-region rational-function (augmented Hessian)
C orbital step over all orbital pairs p<q (oo, ov and vv rotations,
C fully coupled). One call per macro iteration (xvscf/xvtran/xintprc/
C xpccd); the state between calls is kept in JOBARC (RFO_* records).
C Uphill steps are rejected: the previous orbitals are restored and the
C trust radius is reduced. At every convergence the lowest eigenvalues
C of the exact orbital Hessian are found (Pccd_rsp_stab; without it by
C Davidson iterations with finite-difference Hessian-vector products of
C the analytic gradient, one per macro iteration, Pccd_rfo_stab). For a
C negative eigenvalue the optimization continues along that mode (a
C gradient-only method can never leave e.g. a symmetric saddle point).
C
C Parametrization: C(new) = C * exp(K), K(p,q) = x(pq) = -K(q,p), p<q,
C dE/dx(pq) = -2*(G(p,q)-G(q,p)) with G = OBRTGRDA (verified against
C finite-difference energies, 10/2026).
C
C With SYMMETRY=ON only rotations between MOs of the same irrep are
C parameters (Irp(p) = irrep of MO p; SCFEVCA0 columns and OBRTGRDA
C are both ordered occupied by irrep, then virtual by irrep). With
C SYMMETRY=OFF all Irp(p) = 1.
C
C Hessian (file pccd_rfo_h0.in, see Pccd_rfo_h0mode): default (EXACT)
C at every accepted point the exact Hessian, the frozen-amplitude
C Hessian plus the amplitude response (Pccd_rsp; pCCD only), dense up
C to Npar = 1500, above that applied in the Davidson iterations of the
C step; without it (OO-CCD/LCCD, memory, failed checks) as FRZ. FRZ:
C at every accepted point the analytic frozen-amplitude Hessian
C of that point (Pccd_ohess, OHESFRZ; pCCD only) plus a correction
C (the amplitude response) built by the SR1 updates; FRZ0: the frozen
C Hessian at the first point only; CONST: 4|df| (ov) and a constant
C (oo, vv) at the first point. Without OHESFRZ (OO-CCD/LCCD, too
C little memory) the constant start is used. Benchmark 10/2026 (H2O,
C N2 cc-pVDZ 1.1-2.5 A): FRZ needs 2-12x fewer gradients than CONST,
C same minima except N2 2.5 A (another local minimum).
C
C RHF, no dropped MOs.
C
      Implicit Double Precision(A-H,O-Z)
      Logical Grdcnv,Symmetry
      Dimension Work(Maxcor)
      Common/Symm/Symmetry


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



C Mxd: Davidson subspace size (Pccd_rfo_dav).
      Parameter (Mxd=40)

      If (Iuhf .Ne. 0) Then
         Write(6,"(2a)") " @-Pccd_rfo_step: OO_METHOD=RFO requires",
     +                   " REF=RHF."
         Call Errex
      Endif
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nao)
      If (Nao .Ne. Nbas) Then
         Write(6,"(2a)") " @-Pccd_rfo_step: number of MOs and AOs",
     +                   " differ (dropped/linear dependent MOs)."
         Call Errex
      Endif

      Nsq  = Nbas*Nbas
      Irpo = 1
      Iirp = Irpo + 8
      Ibeg = Iirp + Nbas
      If (Ibeg .Gt. Maxcor) Call Insmem("Pccd_rfo_step",Ibeg,Maxcor)
      Call Pccd_rfo_irp(Work(Irpo),Work(Iirp),Nbas,Symmetry)
      Call Pccd_rfo_npar(Work(Iirp),Nbas,Npar)
      Np1  = Npar + 1
      Mxhv = Min(Npar,60)
      Lwork = Max(64*Np1,3*Np1)
C Large Npar (Davidson) with the exact Hessian and the frozen Hessian
C available (Lhes = 0): no model Hessian H and no dense augmented
C matrix Ah, which then only holds the Davidson workspace; the frozen
C Hessian is the only Npar**2 array.
      Call Pccd_rfo_h0mode(Ih0)
      Call Pccd_rfo_eigmode(Npar,Idav)
      Call Getrec(0,"JOBARC","OHESOK  ",Length,Junk)
      Iohs = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","OHESOK  ",1,Iohs)
      Lhes = 1
      If (Idav .Eq. 1 .And. Ih0 .Eq. 3 .And. Iohs .Eq. 1) Lhes = 0
      Lh  = Npar*Npar
      Lah = Np1*Np1
      If (Lhes .Eq. 0) Then
         Lh  = 1
         Lah = 2*Np1*Mxd + 2*Np1 + 2*Mxd*Mxd + 65*Mxd + 1
      Endif
      I000 = Ibeg
      I010 = I000 + Nsq
      I020 = I010 + Nsq
      I030 = I020 + Nsq
      I040 = I030 + Nsq
      I050 = I040 + Nsq
      I060 = I050 + Nsq
      I070 = I060 + Npar
      I080 = I070 + Npar
      I090 = I080 + Npar
      I100 = I090 + Npar
      I110 = I100 + Lh
      I120 = I110 + Lah
      I130 = I120 + Np1
      I140 = I130 + Lwork
      I150 = I140 + Npar*Mxhv
      I160 = I150 + Npar*Mxhv
      I170 = I160 + Mxhv*Mxhv
      I180 = I170 + Mxhv
      I190 = I180 + Nsq
      I200 = I190 + Npar
      I210 = I200 + Npar
      I220 = I210 + Npar*Npar
      Iend = I220
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_rfo_step",Iend,Maxcor)
C Workspace of the exact-Hessian stability check (Pccd_rsp), if it fits.
      Call Pccd_rsp_size(Nbas,Nocc,Npar,Lrsp)
      If (I220+Lrsp .Gt. Maxcor) Lrsp = 0

      Call Pccd_rfo_do(Work(I000),Work(I010),Work(I020),Work(I030),
     +                 Work(I040),Work(I050),Work(I060),Work(I070),
     +                 Work(I080),Work(I090),Work(I100),Work(I110),
     +                 Work(I120),Work(I130),Lwork,Work(I140),
     +                 Work(I150),Work(I160),Work(I170),Work(I180),
     +                 Work(I190),Work(I200),Work(I210),Work(Iirp),
     +                 Ih0,Mxhv,Nbas,Nocc,Npar,Grdcnv,Tol,Etol,Gmax,
     +                 Grms,Ncycle,Work(I220),Lrsp,Lah,Lhes)
C The T2 guess of the next macro iterations follows accepted points only.
      Call Getrec(20,"JOBARC","RFO_ACPT",1,Iacpt)
      If (Iacpt .Eq. 1) Call Pccd_rfo_tguess(Work,Maxcor,"PROMO")
      Return
      End


      Subroutine Pccd_rfo_irp(Nocirr,Irp,Nbas,Symmetry)
C Irrep of each MO in the SCFEVCA0/OBRTGRDA order: occupied by irrep,
C then virtual by irrep. All 1 with SYMMETRY=OFF.
      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry
      Dimension Nocirr(8),Irp(Nbas),Nbfirr(8)
      If (.Not. Symmetry) Then
         Do I = 1, Nbas
            Irp(I) = 1
         Enddo
         Return
      Endif
C OCCUPYA may be stored with 8 entries; NUMBASIR has one per irrep.
      Call Getrec(0,"JOBARC","NUMBASIR",Nirrep,Junk)
      Call Getrec(0,"JOBARC","OCCUPYA ",Length,Junk)
      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Nbfirr)
      Call Getrec(20,"JOBARC","OCCUPYA ",Length,Nocirr)
      K = 0
      Do Irr = 1, Nirrep
         Do I = 1, Nocirr(Irr)
            K = K + 1
            Irp(K) = Irr
         Enddo
      Enddo
      Do Irr = 1, Nirrep
         Do I = 1, Nbfirr(Irr) - Nocirr(Irr)
            K = K + 1
            Irp(K) = Irr
         Enddo
      Enddo
      If (K .Ne. Nbas) Then
         Write(6,"(a,2i6)") " @-Pccd_rfo_irp: MO count mismatch",K,Nbas
         Call Errex
      Endif
      Return
      End

      Subroutine Pccd_rfo_npar(Irp,Nbas,Npar)
C Number of orbital rotation parameters: pairs p<q of the same irrep.
      Dimension Irp(Nbas)
      Npar = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            If (Irp(Ip) .Eq. Irp(Iq)) Npar = Npar + 1
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rfo_h0mode(Ih0)
C Model Hessian of OO_METHOD=RFO from the first word of the file
C pccd_rfo_h0.in: CONST -> 0, FRZ0 -> 1, FRZ -> 2; EXACT or no file
C -> 3.
      Character*8 Word
      Logical Exist
      Ih0 = 3
      Inquire(File="pccd_rfo_h0.in",Exist=Exist)
      If (.Not. Exist) Return
      Word = " "
      Open(Unit=44,File="pccd_rfo_h0.in",Status="Old")
      Read(44,*,End=10,Err=10) Word
 10   Close(44)
      If (Word .Eq. "CONST") Ih0 = 0
      If (Word .Eq. "FRZ0") Ih0 = 1
      If (Word .Eq. "FRZ") Ih0 = 2
      If (Word .Eq. "EXACT") Ih0 = 3
      Return
      End

      Subroutine Pccd_rfo_eigmode(Npar,Idav)
C Eigenproblems of OO_METHOD=RFO (step, stability start vectors): dense
C DSYEV (Idav = 0) up to Npar = Ndense, Davidson (Idav = 1) above. The
C first word of the file pccd_rfo_eig.in, DENSE or DAVIDSON, forces one.
      Parameter (Ndense=1500)
      Character*8 Word
      Logical Exist
      Idav = 0
      If (Npar .Gt. Ndense) Idav = 1
      Inquire(File="pccd_rfo_eig.in",Exist=Exist)
      If (.Not. Exist) Return
      Word = " "
      Open(Unit=44,File="pccd_rfo_eig.in",Status="Old")
      Read(44,*,End=10,Err=10) Word
 10   Close(44)
      If (Word .Eq. "DENSE") Idav = 0
      If (Word .Eq. "DAVIDSON") Idav = 1
      Return
      End
