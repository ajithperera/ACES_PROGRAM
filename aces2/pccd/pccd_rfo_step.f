










      Subroutine Pccd_rfo_step(Work,Maxcor,Nbas,Nocc,Iuhf,Grdcnv,Tol,
     +                         Gmax,Grms,Ncycle)
C
C OO_METHOD=RFO: trust-region rational-function (augmented Hessian)
C orbital step with an SR1-updated model Hessian over all orbital pairs
C p<q (oo, ov and vv rotations, fully coupled). One call per macro
C iteration (xvscf/xvtran/xintprc/xpccd); the state between calls is
C kept in JOBARC (RFO_* records). Uphill steps are rejected: the
C previous orbitals are restored and the trust radius is reduced.
C
C Parametrization: C(new) = C * exp(K), K(p,q) = x(pq) = -K(q,p), p<q,
C dE/dx(pq) = -2*(G(p,q)-G(q,p)) with G = OBRTGRDA (verified against
C finite-difference energies, 10/2026).
C
C RHF, SYMMETRY=OFF, no dropped MOs.
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




      If (Iuhf .Ne. 0 .Or. Symmetry) Then
         Write(6,"(2a)") " @-Pccd_rfo_step: OO_METHOD=RFO requires",
     +                   " REF=RHF and SYMMETRY=OFF."
         Call Errex
      Endif
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nao)
      If (Nao .Ne. Nbas) Then
         Write(6,"(2a)") " @-Pccd_rfo_step: number of MOs and AOs",
     +                   " differ (dropped/linear dependent MOs)."
         Call Errex
      Endif

      Npar = Nbas*(Nbas-1)/2
      Nsq  = Nbas*Nbas
      Np1  = Npar + 1
      Lwork = Max(64*Np1,3*Np1)
      I000 = 1
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
      I110 = I100 + Npar*Npar
      I120 = I110 + Np1*Np1
      I130 = I120 + Np1
      Iend = I130 + Lwork
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_rfo_step",Iend,Maxcor)

      Call Pccd_rfo_do(Work(I000),Work(I010),Work(I020),Work(I030),
     +                 Work(I040),Work(I050),Work(I060),Work(I070),
     +                 Work(I080),Work(I090),Work(I100),Work(I110),
     +                 Work(I120),Work(I130),Lwork,Nbas,Nocc,Npar,
     +                 Grdcnv,Tol,Gmax,Grms,Ncycle)
      Return
      End

