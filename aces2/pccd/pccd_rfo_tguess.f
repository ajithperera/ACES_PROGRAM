















































































































































































































      Subroutine Pccd_rfo_tguess(Work,Maxcor,Mode)
C
C OO_METHOD=RFO amplitude continuation. Away from equilibrium the pCCD
C equations have several solutions at the same orbitals and, starting
C each macro iteration from the MP2 guess, DIIS lands on one or another
C of them depending on round-off (N2 cc-pVDZ, R=2.5 A: two branches
C 3.4 mEh apart, picked at random for orbital steps of 1e-17). The
C orbital optimizer then sees a discontinuous energy and the trust
C radius collapses. Starting T2 from the amplitudes of the last
C accepted orbitals keeps the solver on one branch. Mode:
C   "LOAD ": if RFO_T2C exists, overwrite T2(AB) (list 46) with it.
C            Called before the T2 resorts in Pccd_rdriver.
C   "SAVE ": converged T2(AB) -> RFO_T2N (after Pccd_rdriver).
C   "PROMO": RFO_T2N -> RFO_T2C (the RFO step accepted this point).
C No-op unless OO_METHOD=RFO. RHF (list 46 only).
C
      Implicit Double Precision(A-H,O-Z)
      Character*5 Mode
      Dimension Work(Maxcor)


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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end

      If (Iflags2(178) .Ne. 6) Return
      Nsize = Isymsz(Isytyp(1,46),Isytyp(2,46))
      If (Nsize .Gt. Maxcor) Call Insmem("Pccd_rfo_tguess",Nsize,Maxcor)

      If (Mode .Eq. "LOAD ") Then
         Call Getrec(0,"JOBARC","RFO_T2C ",Length,Junk)
         If (Length .Le. 0) Return
         If (Length .Ne. Nsize*Iintfp) Then
            Write(6,"(a)") " @-Pccd_rfo_tguess: RFO_T2C size mismatch."
            Call Errex
         Endif
         Call Getrec(20,"JOBARC","RFO_T2C ",Nsize*Iintfp,Work)
         Call Putall(Work,Nsize,1,46)
         Write(6,"(2a)") " @-RFO T2 guess: amplitudes of the last",
     +                   " accepted orbitals."
      Elseif (Mode .Eq. "SAVE ") Then
         Call Getall(Work,Nsize,1,46)
         Call Putrec(20,"JOBARC","RFO_T2N ",Nsize*Iintfp,Work)
      Elseif (Mode .Eq. "PROMO") Then
         Call Getrec(0,"JOBARC","RFO_T2N ",Length,Junk)
         If (Length .Le. 0) Return
         Call Getrec(20,"JOBARC","RFO_T2N ",Nsize*Iintfp,Work)
         Call Putrec(20,"JOBARC","RFO_T2C ",Nsize*Iintfp,Work)
      Endif
      Return
      End
