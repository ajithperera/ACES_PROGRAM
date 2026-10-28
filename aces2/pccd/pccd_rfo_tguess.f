

























































































































































































































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
C   "LOAD ": if RFO_T2C exists, overwrite T2(AB) (list 46) with it;
C            otherwise (first macro iteration) with T2 from a TGUESS
C            file of an earlier calculation, if present and of the
C            right size (geometry continuation: copy TGUESS with the
C            orbitals, oorot LOAD). Called before the T2 resorts in
C            Pccd_rdriver.
C   "SAVE ": converged T2(AB) -> RFO_T2N (after Pccd_rdriver).
C   "PROMO": RFO_T2N -> RFO_T2C (the RFO step accepted this point) and
C            to the TGUESS file (DMPTGSS layout, RHF: T1 = 0, list 44
C            slot = 0, list 46), so that the file holds the T2 of the
C            last accepted point (the minimum at the end).
C No-op unless OO_METHOD=RFO. RHF (list 46 only).
C
      Implicit Double Precision(A-H,O-Z)
      Character*5 Mode
      Character*80 Fulnam
      Logical Exist
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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
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
C TGUESS layout (Pccd_ddmptgss/DMPTGSS, RHF): T1, T2(A<B,I<J), T2(Ab,Ij)
      N44  = Isymsz(1,3)
      Ltg  = Nt(1) + N44 + Nsize
      If (Nsize+Ltg .Gt. Maxcor) Call Insmem("Pccd_rfo_tguess",
     +                                       Nsize+Ltg,Maxcor)

      If (Mode .Eq. "LOAD ") Then
         Call Getrec(0,"JOBARC","RFO_T2C ",Length,Junk)
         If (Length .Le. 0) Then
            Call Gfname("TGUESS  ",Fulnam,Namlen)
            Inquire(File=Fulnam(1:Namlen),Exist=Exist,Size=Isz)
            If (.Not. Exist) Return
C Unformatted sequential: one record, 4-byte markers on both sides.
            If (Isz .Ne. 8*Ltg+8) Then
               Write(6,"(2a)") " @-RFO TGUESS file ignored: size does",
     +              " not match this calculation."
               Return
            Endif
            Call Rdtgss(Work,Ltg,"TGUESS  ",Exist)
            Call Putall(Work(Nt(1)+N44+1),Nsize,1,46)
            Write(6,"(2a)") " @-RFO T2 guess: amplitudes from the",
     +                      " TGUESS file."
            Return
         Endif
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
         Call Dzero(Work(Nsize+1),Nt(1)+N44)
         Call Dcopy(Nsize,Work,1,Work(Nsize+Nt(1)+N44+1),1)
         Call Dmptgss(Work(Nsize+1),Ltg,"TGUESS  ")
      Endif
      Return
      End
