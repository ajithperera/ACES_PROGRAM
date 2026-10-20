#include "flags.h"
      Subroutine Pccd_fc_setup(Xpccd,Force)
C
C Frozen (dropped) MOs with OO_METHOD=RFO pCCD: uncorrelated but
C optimized. The DROPMO list (MODROPA, energy-ordered indices from
C xjoda) is turned into the number of frozen occupied and frozen
C virtual MOs per irrep, PCCDNFZO and PCCDNFZV, and NUMDROPA is set to
C zero so that xvtran/xintprc transform and xpccd works in the full MO
C space. xpccd (Pccd_fcmask) keeps the amplitudes of pairs with a
C frozen index at zero, so the energy is that of pCCD in the active
C space with the frozen MOs in the reference, and the orbital gradient
C and Hessian of the full space include the frozen-active and frozen-
C virtual rotations.
C
C JOBARC must be open (Aces_ja_init by the caller): do_brueckner, and
C xoorot FCSET for drivers that run the modules themselves (oonewt.py;
C Force: also with another OO_METHOD, e.g. ROTGRAD for gradient checks).
C Called after the cycle-1 SCF (HF orbital energies, SCF order). The
C frozen MOs are then the lowest occupied/highest virtual MOs of their
C irreps; the RFO steps keep the column order within an irrep, so they
C stay the first NFZO(irrep) occupied and the last NFZV(irrep) virtual
C MOs of each irrep.
C
      Implicit Double Precision(A-H,O-Z)
      Logical Xpccd,Force
#include "maxbasfn.par"
      Dimension Eval(MAXBASFN),Isrt(MAXBASFN),Idrop(MAXBASFN),
     +          Irp(MAXBASFN),Ipos(MAXBASFN),Nbfirr(8),Nocirr(8),
     +          Nfzo(8),Nfzv(8)
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

      If (.Not. Xpccd) Return
      If (Iflags2(h_IFLAGS2_oo_method).Ne.6 .And. .Not.Force) Return
      Call Getrec(0,"JOBARC","NUMDROPA",Length,Junk)
      Ndrop = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","NUMDROPA",1,Ndrop)
      If (Ndrop .Eq. 0) Then
C No DROPMO in this job: clear frozen counts left in a reused JOBARC.
         Call Getrec(0,"JOBARC","PCCDNFZO",Length,Junk)
         If (Length .Gt. 0) Then
            Do Irr = 1, 8
               Nfzo(Irr) = 0
            Enddo
            Call Putrec(20,"JOBARC","PCCDNFZO",8,Nfzo)
            Call Putrec(20,"JOBARC","PCCDNFZV",8,Nfzo)
         Endif
         Return
      Endif
      If (Iflags(h_IFLAGS_calc) .Ne. 52) Then
         Write(6,"(2a)") " @-Pccd_fc_setup: dropped MOs with",
     +         " OO_METHOD=RFO are only implemented for CALC=PCCD."
         Call Errex
      Endif
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nbas)
      Call Getrec(0,"JOBARC","NUMBASIR",Nirrep,Junk)
      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Nbfirr)
      Call Getrec(0,"JOBARC","OCCUPYA ",Length,Junk)
      Call Getrec(20,"JOBARC","OCCUPYA ",Length,Nocirr)
      Call Getrec(20,"JOBARC","SCFEVLA0",Nbas*Iintfp,Eval)
      Call Getrec(20,"JOBARC","MODROPA ",Ndrop,Idrop)

C Irrep and position within the irrep of each MO (SCF order).
      K = 0
      Do Irr = 1, Nirrep
         Do I = 1, Nbfirr(Irr)
            K = K + 1
            Irp(K)  = Irr
            Ipos(K) = I
         Enddo
      Enddo
C Energy order -> SCF index (as vtran/dropvc.F, stable sort).
      Do I = 1, Nbas
         Isrt(I) = I
      Enddo
      Do I = 2, Nbas
         Ik = Isrt(I)
         J  = I - 1
 20      If (J .Ge. 1) Then
            If (Eval(Isrt(J)) .Gt. Eval(Ik)) Then
               Isrt(J+1) = Isrt(J)
               J = J - 1
               Go to 20
            Endif
         Endif
         Isrt(J+1) = Ik
      Enddo

      Do Irr = 1, 8
         Nfzo(Irr) = 0
         Nfzv(Irr) = 0
      Enddo
      Do I = 1, Ndrop
         If (Idrop(I) .Lt. 1 .Or. Idrop(I) .Gt. Nbas) Then
            Write(6,"(a,i6)") " @-Pccd_fc_setup: bad DROPMO index",
     +                        Idrop(I)
            Call Errex
         Endif
         K   = Isrt(Idrop(I))
         Irr = Irp(K)
         If (Ipos(K) .Le. Nocirr(Irr)) Then
            Nfzo(Irr) = Nfzo(Irr) + 1
         Else
            Nfzv(Irr) = Nfzv(Irr) + 1
         Endif
      Enddo
C The frozen MOs must be the lowest occupied and highest virtual MOs
C of each irrep.
      Do I = 1, Ndrop
         K   = Isrt(Idrop(I))
         Irr = Irp(K)
         If (Ipos(K) .Le. Nocirr(Irr)) Then
            Iok = 0
            If (Ipos(K) .Le. Nfzo(Irr)) Iok = 1
         Else
            Iok = 0
            If (Ipos(K) .Gt. Nbfirr(Irr)-Nfzv(Irr)) Iok = 1
         Endif
         If (Iok .Eq. 0) Then
            Write(6,"(2a)") " @-Pccd_fc_setup: only the lowest",
     +           " occupied and the highest virtual MOs of each irrep"
            Write(6,"(a)") " can be frozen with OO_METHOD=RFO."
            Call Errex
         Endif
      Enddo

      Call Putrec(20,"JOBARC","PCCDNFZO",8,Nfzo)
      Call Putrec(20,"JOBARC","PCCDNFZV",8,Nfzv)
      Call Putrec(20,"JOBARC","NUMDROPA",1,0)
      Call Putrec(20,"JOBARC","NUMDROPB",1,0)
      Write(6,"(2a)") " Dropped MOs with OO_METHOD=RFO: uncorrelated",
     +                " but optimized (frozen pairs)."
      Write(6,"(a,8i4)") " Frozen occupied MOs per irrep:",
     +                   (Nfzo(Irr),Irr=1,Nirrep)
      Write(6,"(a,8i4)") " Frozen virtual MOs per irrep :",
     +                   (Nfzv(Irr),Irr=1,Nirrep)
      Return
      End
