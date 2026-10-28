#include "flags.h"
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
C   DUMPREC: the double precision JOBARC record named on the second
C           line of oorot.in (all of it)    -> rec_<name>.dat (n, values)
C   FCSET : DROPMO -> frozen pairs of OO_METHOD=RFO pCCD (as xrunpccd
C           does after the cycle-1 SCF; for drivers running the modules)
C All text files are column-major with a "nrow ncol" header line.
C
      Implicit Double Precision (A-H,O-Z)
      Character*8 Mode,Recnam
      Character*24 Fname
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "info.com"
      Call Aces_init(icore, i0, icrsiz, iuhf, .true.)
      Open(Unit=31,File="oorot.in",Status="Old")
      Read(31,"(a)") Mode
      Recnam = " "
      If (Mode(1:7) .Eq. "DUMPREC") Read(31,"(a)") Recnam
      Close(31)
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nao)
      Nmo  = Nocco(1) + Nvrto(1)
      I000 = i0
      I010 = I000 + Nao*Nmo*Iintfp
      I020 = I010 + Nmo*Nmo*Iintfp
      Iend = I020 + Nao*Nmo*Iintfp
      If (Iend-i0 .Gt. icrsiz) Call Insmem("oorot",Iend-i0,icrsiz)
      If (Mode(1:7) .Eq. "DUMPREC") Then
         Call Getrec(0,"JOBARC",Recnam,Length,Junk)
         If (Length .Le. 0) Then
            Write(6,"(2a)") " @-oorot: no record ",Recnam
            Call Errex
         Endif
         If (Length .Gt. icrsiz) Call Insmem("oorot",Length,icrsiz)
         Call Getrec(20,"JOBARC",Recnam,Length,icore(i0))
         Write(Fname,"(a4,a,a4)") "rec_",Recnam(1:Len_trim(Recnam)),
     +                            ".dat"
         Call Oorot_wrec(icore(i0),Length/Iintfp,Fname)
         Call Aces_fin
         Stop
      Endif
      If (Mode(1:5) .Eq. "FCSET") Then
         Call Pccd_fc_setup(.True.,.True.)
         Call Aces_fin
         Stop
      Endif
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


      Subroutine Oorot_wrec(X,N,Fname)
      Implicit Double Precision (A-H,O-Z)
      Character*24 Fname
      Dimension X(N)
      Open(Unit=32,File=Fname,Status="Unknown")
      Write(32,*) N
      Write(32,"(4(1x,E25.17))") (X(I),I=1,N)
      Close(32)
      Return
      End
