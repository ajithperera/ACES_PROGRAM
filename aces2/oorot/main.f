














































































































































































































      Program Oorot
C
C Pilot utility for the OO-pCCD exact-Newton (finite-difference Hessian)
C test. The mode is read from the first line of oorot.in:
C   SAVE  : SCFEVCA0 (AO x MO)            -> c0.dat
C   APPLY : C = C0*U, C0 from c0.dat,
C           U (MO x MO) from u.dat         -> SCFEVCA0, SCFEVCB0
C   DUMP  : OBRTGRDA (MO x MO), TOTENERG   -> grad.dat
C   SAVEAB: SCFEVCA0, SCFEVCB0            -> ca.dat, cb.dat
C   LOADAB: ca.dat, cb.dat -> SCFEVCA0/B0, SCFEVCAS/BS (UHF guess)
C   LOAD  : c0.dat -> SCFEVCA0/B0 + non-HF flags (OO-pCCD guess)
C   DUMPD : DENSOO, DENSVV (correlation 1-RDM blocks)  -> dens.dat
C   DUMPT2: max|t| and norm of RFO_T2C (last accepted T2) -> stdout
C All text files are column-major with a "nrow ncol" header line.
C
      Implicit Double Precision (A-H,O-Z)
      Character*8 Mode


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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



c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
      Call Aces_init(icore, i0, icrsiz, iuhf, .true.)
      Open(Unit=31,File="oorot.in",Status="Old")
      Read(31,"(a)") Mode
      Close(31)
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nao)
      Nmo  = Nocco(1) + Nvrto(1)
      I000 = i0
      I010 = I000 + Nao*Nmo*Iintfp
      I020 = I010 + Nmo*Nmo*Iintfp
      Iend = I020 + Nao*Nmo*Iintfp
      If (Iend-i0 .Gt. icrsiz) Call Insmem("oorot",Iend-i0,icrsiz)
      If (Mode(1:6) .Eq. "DUMPT2") Then
C Largest |t| and norm of the RFO continuation amplitudes RFO_T2C
         Call Getrec(0,"JOBARC","RFO_T2C ",Length,Junk)
         If (Length .Le. 0) Then
            Write(6,"(a)") " @-oorot: no RFO_T2C record."
            Call Errex
         Endif
         If (Length .Gt. icrsiz) Call Insmem("oorot",Length,icrsiz)
         Call Getrec(20,"JOBARC","RFO_T2C ",Length,icore(i0))
         Call Oorot_t2st(icore(i0),Length/Iintfp)
         Call Aces_fin
         Stop
      Endif
      Call Oorot_do(icore(I000),icore(I010),icore(I020),Nao,Nmo,Mode)
      Call Aces_fin
      Stop
      End


      Subroutine Oorot_t2st(T,N)
      Implicit Double Precision (A-H,O-Z)
      Dimension T(N)
      Tmax = 0.0D0
      Tnrm = 0.0D0
      Do I = 1, N
         Tmax = Max(Tmax,Dabs(T(I)))
         Tnrm = Tnrm + T(I)*T(I)
      Enddo
      Write(6,"(a,i8,a,F12.6,a,F12.6)") " @-oorot T2: n =",N,
     +      "  max|t| =",Tmax,"  |T2| =",Dsqrt(Tnrm)
      Return
      End
