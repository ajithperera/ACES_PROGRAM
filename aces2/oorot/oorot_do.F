#include "flags.h"
      Subroutine Oorot_do(C0,U,C,Nao,Nmo,Mode)
      Implicit Double Precision (A-H,O-Z)
      Character*8 Mode
      Dimension C0(Nao,Nmo),U(Nmo,Nmo),C(Nao,Nmo)
#include "machsp.com"
#include "flags.com"
      If (Mode(1:4) .Eq. "SAVE") Then
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
