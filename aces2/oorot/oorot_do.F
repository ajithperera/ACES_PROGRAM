#include "flags.h"
      Subroutine Oorot_do(C0,U,C,Nao,Nmo,Mode)
      Implicit Double Precision (A-H,O-Z)
      Character*8 Mode
      Dimension C0(Nao,Nmo),U(Nmo,Nmo),C(Nao,Nmo)
#include "machsp.com"
#include "flags.com"
      If (Mode(1:6) .Eq. "SAVEAB") Then
C Alpha and beta MOs (e.g. of a broken-symmetry UHF) -> ca.dat, cb.dat
         Call Getrec(20,"JOBARC","SCFEVCA0",Nao*Nmo*Iintfp,C0)
         Call Getrec(20,"JOBARC","SCFEVCB0",Nao*Nmo*Iintfp,C)
         Open(Unit=32,File="ca.dat",Status="Unknown")
         Write(32,*) Nao, Nmo
         Write(32,"(4(1x,E25.17))") ((C0(I,J),I=1,Nao),J=1,Nmo)
         Close(32)
         Open(Unit=32,File="cb.dat",Status="Unknown")
         Write(32,*) Nao, Nmo
         Write(32,"(4(1x,E25.17))") ((C(I,J),I=1,Nao),J=1,Nmo)
         Close(32)
      Elseif (Mode .Eq. "SAVE    ") Then
         Call Getrec(20,"JOBARC","SCFEVCA0",Nao*Nmo*Iintfp,C0)
         Open(Unit=32,File="c0.dat",Status="Unknown")
         Write(32,*) Nao, Nmo
         Write(32,"(4(1x,E25.17))") ((C0(I,J),I=1,Nao),J=1,Nmo)
         Close(32)
      Elseif (Mode(1:5) .Eq. "APPLY") Then
         Open(Unit=32,File="c0.dat",Status="Old")
         Read(32,*) N1, N2
         If (N1 .Ne. Nao .Or. N2 .Ne. Nmo) Call Errex
         Read(32,*) ((C0(I,J),I=1,Nao),J=1,Nmo)
         Close(32)
         Open(Unit=33,File="u.dat",Status="Old")
         Read(33,*) N1, N2
         If (N1 .Ne. Nmo .Or. N2 .Ne. Nmo) Call Errex
         Read(33,*) ((U(I,J),I=1,Nmo),J=1,Nmo)
         Close(33)
         Call Dgemm("N","N",Nao,Nmo,Nmo,1.0D0,C0,Nao,U,Nmo,0.0D0,
     +              C,Nao)
         Call Putrec(20,"JOBARC","SCFEVCA0",Nao*Nmo*Iintfp,C)
         Call Putrec(20,"JOBARC","SCFEVCB0",Nao*Nmo*Iintfp,C)
      Elseif (Mode(1:5) .Eq. "DUMP ") Then
         Call Getrec(20,"JOBARC","OBRTGRDA",Nmo*Nmo*Iintfp,U)
         Call Getrec(20,"JOBARC","TOTENERG",Iintfp,E)
         Open(Unit=34,File="grad.dat",Status="Unknown")
         Write(34,"(E25.17)") E
         Write(34,*) Nmo, Nmo
         Write(34,"(4(1x,E25.17))") ((U(I,J),I=1,Nmo),J=1,Nmo)
         Close(34)
      Elseif (Mode(1:4) .Eq. "NHF0") Then
C Reset the NONHF flag in the JOBARC copy of IFLAGS (set to 1 by xpccd's
C pccd_do_orbrots after its first cycle).
         Call Getrec(20,"JOBARC","IFLAGS  ",100,Iflags)
         Write(6,"(a,i3)") " @-oorot: NONHF flag was ",
     +                     Iflags(h_IFLAGS_nonhf)
         Iflags(h_IFLAGS_nonhf) = 0
         Call Putrec(20,"JOBARC","IFLAGS  ",100,Iflags)
      Elseif (Mode(1:6) .Eq. "LOADAB") Then
C UHF guess: ca.dat -> SCFEVCA0, cb.dat -> SCFEVCB0, GUESS=JOBARC, so
C the next xvscf iterates from them (broken-symmetry UHF guesses).
         Open(Unit=32,File="ca.dat",Status="Old")
         Read(32,*) N1, N2
         If (N1 .Ne. Nao .Or. N2 .Ne. Nmo) Call Errex
         Read(32,*) ((C0(I,J),I=1,Nao),J=1,Nmo)
         Close(32)
         Open(Unit=32,File="cb.dat",Status="Old")
         Read(32,*) N1, N2
         If (N1 .Ne. Nao .Or. N2 .Ne. Nmo) Call Errex
         Read(32,*) ((C(I,J),I=1,Nao),J=1,Nmo)
         Close(32)
C xvscf (Eigvectors_frm_jobarc) starts from SCFEVCAS/SCFEVCBS unless
C Brueckner, HFSTABILITY=FOLLOW or OO_METHOD is set: write both pairs.
         Call Putrec(20,"JOBARC","SCFEVCA0",Nao*Nmo*Iintfp,C0)
         Call Putrec(20,"JOBARC","SCFEVCB0",Nao*Nmo*Iintfp,C)
         Call Putrec(20,"JOBARC","SCFEVCAS",Nao*Nmo*Iintfp,C0)
         Call Putrec(20,"JOBARC","SCFEVCBS",Nao*Nmo*Iintfp,C)
         Call Getrec(20,"JOBARC","IFLAGS  ",100,Iflags)
         Iflags(h_IFLAGS_guess) = 0
         Call Putrec(20,"JOBARC","IFLAGS  ",100,Iflags)
         Write(6,"(a)") " @-oorot: UHF guess loaded from ca/cb.dat."
      Elseif (Mode(1:4) .Eq. "LOAD") Then
C Restricted guess orbitals from c0.dat -> SCFEVCA0, SCFEVCB0, and the
C flags Pccd_do_orbrots sets for non-HF orbitals, so that the next
C xvscf (no iterations, MOs from JOBARC) builds the matching Fock.
         Open(Unit=32,File="c0.dat",Status="Old")
         Read(32,*) N1, N2
         If (N1 .Ne. Nao .Or. N2 .Ne. Nmo) Call Errex
         Read(32,*) ((C0(I,J),I=1,Nao),J=1,Nmo)
         Close(32)
         Call Putrec(20,"JOBARC","SCFEVCA0",Nao*Nmo*Iintfp,C0)
         Call Putrec(20,"JOBARC","SCFEVCB0",Nao*Nmo*Iintfp,C0)
         Call Getrec(20,"JOBARC","IFLAGS  ",100,Iflags)
         Iflags(h_IFLAGS_scf_maxcyc) = 0
         Iflags(h_IFLAGS_nonhf)      = 1
         Iflags(h_IFLAGS_qrhf_orb)   = 0
         Iflags(h_IFLAGS_guess)      = 0
         Call Putrec(20,"JOBARC","IFLAGS  ",100,Iflags)
         Write(6,"(a)") " @-oorot: guess MOs loaded from c0.dat."
      Elseif (Mode(1:5) .Eq. "DUMPF") Then
C AO Fock matrix (Nao x Nao, Nao = Nmo assumed) and current MOs
         If (Nao .Ne. Nmo) Call Errex
         Call Getrec(20,"JOBARC","FOCKA   ",Nao*Nao*Iintfp,U)
         Call Getrec(20,"JOBARC","SCFEVCA0",Nao*Nmo*Iintfp,C)
         Open(Unit=35,File="fock.dat",Status="Unknown")
         Write(35,*) Nao, Nao
         Write(35,"(4(1x,E25.17))") ((U(I,J),I=1,Nao),J=1,Nao)
         Write(35,*) Nao, Nmo
         Write(35,"(4(1x,E25.17))") ((C(I,J),I=1,Nao),J=1,Nmo)
         Close(35)
      Else
         Write(6,"(2a)") " @-oorot: unknown mode ", Mode
         Call Errex
      Endif
      Return
      End
