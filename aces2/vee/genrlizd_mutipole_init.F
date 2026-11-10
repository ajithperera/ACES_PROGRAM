      Subroutine Genrlizd_multipole_init(Qml,Oml,Aml,Mql,Dml,Qmr,Omr,
     +                                   Amr,Mqr,Dmr,Qm_2,Om_2,Am_2,
     +                                   Mq_2,Dm_2)
C
C Products of transition moments entering the isotropically averaged
C oscillator strength to second order in the wave vector,
C f = f(mu) + f(m) + f(Q) + f(mu,O) + f(mu,M) (length gauge; e.g.
C Bernadotte, Atkins, Jacob, JCP 137, 204106 (2012)), for a
C transition 0 -> n with omega = E_n - E_0 > 0:
C
C   "l" moments: X_l = <0|X|n>,   "r" moments: X_r = <n|X|0>
C
C (ground state: l = right-hand/TMRIGHT side, r = left/TMLEFT side;
C excited states: 0 is the LOWER state of the pair). The moments are
C those of vee_bd_moments (electron charge included; Am, Mq are of
C the antisymmetric real operators -l and -(1/2){r_a,l_b}, l = r x
C nabla). The products, combined in genrlizd_multipole_fin:
C
C  Dm_2(a)  = Dl(a) Dr(a)
C  Am_2(a)  = Al(a) Ar(a)                  (= -|<n|l_a|0>|^2 exact)
C  Qm_2(ab) = Ql(ab) Qr(ab), traceless Theta (xx,yy,zz,xy,xz,yz);
C             Qm_2(7-9) = Q_aa Q_bb + Q_bb Q_aa (not used by fin)
C  Om_2(c)  = Ol(c) Dr(c) + Dl(c) Or(c), O(c) = r^2 r_c, i.e. the
C             Sum_a O_aac of the raw third moment
C  Mq_2     = e_bad [Dl(d) Mr(ab) - Ml(ab) Dr(d)], the six nonzero
C             e_bad (M(ab) = a+3(b-1), position a, angular momentum b)
C
C 11/2026 (Ajith Perera): Om_2 used the traceless octupoles (whose
C Sum_a O_aac vanishes identically), Dmr for both moments and the y,z
C pairing of Oyyy/Ozzz with mu_x; Mq_2 summed all 27 products
C Mq(i)*Dm(j) (not a rotational invariant) instead of the
C Levi-Civita contraction.
C
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq_2,Mqr,Mql
      Dimension Qmr(6),Omr(3),Amr(3),Dmr(3),Mqr(9)
      Dimension Qml(6),Oml(3),Aml(3),Dml(3),Mql(9)
      Dimension Qm_2(9),Om_2(3),Am_2(3),Dm_2(3),Mq_2(27)

      Do Ic = 1, 3
        Dm_2(Ic) = Dmr(Ic)*Dml(Ic)
        Am_2(Ic) = Amr(Ic)*Aml(Ic)
      Enddo

      Do Ic = 1, 6
        Qm_2(Ic) = Qmr(Ic)*Qml(Ic)
      Enddo
      Qm_2(7) = Qmr(1)*Qml(2)+Qmr(2)*Qml(1)
      Qm_2(8) = Qmr(1)*Qml(3)+Qmr(3)*Qml(1)
      Qm_2(9) = Qmr(2)*Qml(3)+Qmr(3)*Qml(2)

      Do Ic = 1, 3
         Om_2(Ic) = Oml(Ic)*Dmr(Ic) + Dml(Ic)*Omr(Ic)
      Enddo

      Call Dzero(Mq_2,27)
C e_xyz: M(yx) d_z ; e_yzx: M(zy) d_x ; e_zxy: M(xz) d_y  (+)
      Mq_2(1) =  (Dml(3)*Mqr(2) - Mql(2)*Dmr(3))
      Mq_2(2) =  (Dml(1)*Mqr(6) - Mql(6)*Dmr(1))
      Mq_2(3) =  (Dml(2)*Mqr(7) - Mql(7)*Dmr(2))
C e_xzy: M(zx) d_y ; e_yxz: M(xy) d_z ; e_zyx: M(yz) d_x  (-)
      Mq_2(4) = -(Dml(2)*Mqr(3) - Mql(3)*Dmr(2))
      Mq_2(5) = -(Dml(3)*Mqr(4) - Mql(4)*Dmr(3))
      Mq_2(6) = -(Dml(1)*Mqr(8) - Mql(8)*Dmr(1))

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " @-Genrlizd_multipole_init products"
      Write(6,"(a)") " mu*mu"
      Write(6,"(3(1x,F15.10))") (Dm_2(i),i=1,3)
      Write(6,"(a)") " l*l"
      Write(6,"(3(1x,F15.10))") (Am_2(i),i=1,3)
      Write(6,"(a)") " Q*Q"
      Write(6,"(3(1x,F15.10))") (Qm_2(i),i=1,9)
      Write(6,"(a)") " mu*O"
      Write(6,"(3(1x,F15.10))") (Om_2(i),i=1,3)
      Write(6,"(a)") " e*mu*M"
      Write(6,"(3(1x,F15.10))") (Mq_2(i),i=1,6)
#endif

      Return
      End
