      Subroutine Beyond_dipole(Dens,Scr1,Scr2,Nao,Iside)
C
C Ground <-> excited multipole transition moments for the
C beyond-dipole oscillator strength (see multipole.F). Dens is the
C spin-summed AO transition density, scaled exactly as for the
C dipole moments TMRIGHT/TMLEFT (tdens.F). Iside = 1 (right, the
C <0|O|n> moments) or 2 (left, <n|O|0>).
C
C 11/2026 (Ajith Perera): called once per side with the spin-summed
C density (it was called per spin with the unscaled density, so the
C UHF beta call overwrote alpha and the RHF factor 2 and the L
C normalization were missing); the electric dipole computed here is
C saved (DMRIGHT/DMLEFT) so that the cross terms in multipole.F use
C the same integral convention as the other moments.
C
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq
      Dimension Dens(Nao,Nao),Scr1(*),Scr2(Nao,Nao)
      Dimension Qm(6),Om(3),Am(3),Mq(9),Dm(3)
      Character*8 Q_string(2),O_string(2),A_string(2),M_string(2)
      Character*8 D_string(2)

#include "machsp.com"

      Data D_string/"DMRIGHT ","DMLEFT  "/
      Data Q_string/"QMRIGHT ","QMLEFT  "/
      Data O_string/"OMRIGHT ","OMLEFT  "/
      Data A_string/"AMRIGHT ","AMLEFT  "/
      Data M_string/"MQRIGHT ","MQLEFT  "/

      Call Vee_bd_moments(Dens,Scr1,Scr2,Nao,Dm,Qm,Om,Am,Mq)

      Call Putrec(20,"JOBARC",D_string(Iside), 3*Iintfp,Dm)
      Call Putrec(20,"JOBARC",Q_string(Iside), 6*Iintfp,Qm)
      Call Putrec(20,"JOBARC",O_string(Iside), 3*Iintfp,Om)
      Call Putrec(20,"JOBARC",A_string(Iside), 3*Iintfp,Am)
      Call Putrec(20,"JOBARC",M_string(Iside), 9*Iintfp,Mq)

      Return
      End
