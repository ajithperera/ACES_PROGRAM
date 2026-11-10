      Subroutine multipole(Root)
C
C Ground -> excited oscillator strength beyond the electric-dipole
C approximation, from the moments saved by beyond_dipole: the
C right-hand (Iside=1) moments are <0|X|n>, the left-hand ones
C <n|X|0>. See genrlizd_multipole_init/_fin for the formulas.
C
C 11/2026 (Ajith Perera): uses the shared init/fin routines (fixes
C the formulas, see there) and the dipole moments computed with the
C other multipoles (DMRIGHT/DMLEFT); Mq_r/Mq_l were dimensioned 3
C but read with 9 elements.
C
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq_r,Mq_l,Mq_2,M_iso
      Dimension Qm_r(6),Qm_l(6),Om_r(3),Om_l(3)
      Dimension Am_r(3),Am_l(3),Dm_r(3),Dm_l(3)
      Dimension Mq_r(9),Mq_l(9)
      Dimension Qm_2(9),Om_2(3),Am_2(3),Dm_2(3),Mq_2(27)
      Character*8 Q_string(2),O_string(2),A_string(2)
      Character*8 D_string(2),M_string(2)

#include "machsp.com"

      Data D_string/"DMRIGHT ","DMLEFT  "/
      Data Q_string/"QMRIGHT ","QMLEFT  "/
      Data O_string/"OMRIGHT ","OMLEFT  "/
      Data A_string/"AMRIGHT ","AMLEFT  "/
      Data M_string/"MQRIGHT ","MQLEFT  "/

      Call Getrec(20,"JOBARC",D_string(1), 3*Iintfp,Dm_r)
      Call Getrec(20,"JOBARC",D_string(2), 3*Iintfp,Dm_l)
      Call Getrec(20,"JOBARC",Q_string(1), 6*Iintfp,Qm_r)
      Call Getrec(20,"JOBARC",Q_string(2), 6*Iintfp,Qm_l)
      Call Getrec(20,"JOBARC",O_string(1), 3*Iintfp,Om_r)
      Call Getrec(20,"JOBARC",O_string(2), 3*Iintfp,Om_l)
      Call Getrec(20,"JOBARC",A_string(1), 3*Iintfp,Am_r)
      Call Getrec(20,"JOBARC",A_string(2), 3*Iintfp,Am_l)
      Call Getrec(20,"JOBARC",M_string(1), 9*Iintfp,Mq_r)
      Call Getrec(20,"JOBARC",M_string(2), 9*Iintfp,Mq_l)
C
C <0|X|n> (right-hand side) is the "l" argument, <n|X|0> the "r".
C
      Call Genrlizd_multipole_init(Qm_r,Om_r,Am_r,Mq_r,Dm_r,
     +                             Qm_l,Om_l,Am_l,Mq_l,Dm_l,
     +                             Qm_2,Om_2,Am_2,Mq_2,Dm_2)
      Call Genrlizd_multipole_fin(Qm_2,Om_2,Am_2,Dm_2,Mq_2,Root,
     +                            D_iso,A_iso,O_iso,Q_iso,M_iso)

      Return
      End
