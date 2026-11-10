      Subroutine Vee_bd_moments(Dens,Scr1,Scr2,Nao,Dm,Qm,Om,Am,Mq)
C
C Multipole transition moments Sum_pq O(pq) Dens(pq) of a (spin
C summed, AO basis) transition density for the beyond-dipole
C oscillator strengths. All integrals are those of vprops (VPROP
C triangles) and include the electron charge:
C
C   Dm(3)  DIPOLE_x            -r_a
C   Qm(6)  QUAD_xx..yz         -Theta_ab (traceless, Buckingham:
C                              (3r_ar_b-r^2 d_ab)/2; xx,yy,zz,xy,xz,yz)
C   Om(3)  R2R_c               -r^2 r_c   (trace of the 3rd moment;
C                                          the OCTUP* integrals are
C                                          traceless, Sum_a O_aac = 0)
C   Am(3)  AMx                 -l_a,  l = r x nabla  (antisymmetric)
C   Mq(9)  MQab, ab=a+3(b-1)   -(1/2)(r_a l_b + l_b r_a)
C                                                     (antisymmetric)
C
C The two antisymmetric operators are expanded from their packed
C triangles as antisymmetric matrices.
C
C 11/2026 (Ajith Perera): common code of beyond_dipole and
C genrlizd_beyond_dipole. Previously AM and MQ were expanded as
C symmetric matrices, the traceless octupoles were used where the
C r^2 r trace is needed, and the l_b r_a ordering correction to MQ was
C added here with wrong signs for MQXZ/MQYZ; MQ now comes from
C vprops already antisymmetrized.
C
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq
      Dimension Dens(Nao*Nao),Scr1(*),Scr2(Nao,Nao)
      Dimension Dm(3),Qm(6),Om(3),Am(3),Mq(9)
      Character*8 Label_angm(3),Label_dmom(3),Label_octp(3)
      Character*8 Label_quad(6),Label_magd(9)

#include "machsp.com"

      Data Label_dmom/'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Label_quad/'QUAD_XX ','QUAD_YY ','QUAD_ZZ ','QUAD_XY ',
     +                'QUAD_XZ ','QUAD_YZ '/
      Data Label_octp/'R2R_X   ','R2R_Y   ','R2R_Z   '/
      Data Label_angm/'   AMX  ','   AMY  ','   AMZ  '/
      Data Label_magd/'   MQXX ','   MQYX ','   MQZX ','   MQXY ',
     +                '   MQYY ','   MQZY ','   MQXZ ','   MQYZ ',
     +                '   MQZZ '/

      Length = Nao*(Nao+1)/2
      Nao2   = Nao*Nao

      Do Ic = 1, 3
         Call Getrec(20,"JOBARC",Label_dmom(Ic),Length*Iintfp,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Dm(Ic) = Ddot(Nao2,Dens,1,Scr2,1)
      Enddo
      Do Ic = 1, 6
         Call Getrec(20,"JOBARC",Label_quad(Ic),Length*Iintfp,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Qm(Ic) = Ddot(Nao2,Dens,1,Scr2,1)
      Enddo
      Do Ic = 1, 3
         Call Getrec(20,"JOBARC",Label_octp(Ic),Length*Iintfp,Scr1)
         Call Expnd2(Scr1,Scr2,Nao)
         Om(Ic) = Ddot(Nao2,Dens,1,Scr2,1)
      Enddo
      Do Ic = 1, 3
         Call Getrec(20,"JOBARC",Label_angm(Ic),Length*Iintfp,Scr1)
         Call Vee_bd_expnd2a(Scr1,Scr2,Nao)
         Am(Ic) = Ddot(Nao2,Dens,1,Scr2,1)
      Enddo
      Do Ic = 1, 9
         Call Getrec(20,"JOBARC",Label_magd(Ic),Length*Iintfp,Scr1)
         Call Vee_bd_expnd2a(Scr1,Scr2,Nao)
         Mq(Ic) = Ddot(Nao2,Dens,1,Scr2,1)
      Enddo

      Return
      End
