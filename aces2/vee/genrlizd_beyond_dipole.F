      Subroutine Genrlizd_Beyond_dipole(Dens,Scr1,Scr2,Qm,Om,Am,
     +                                  Dm,Mq,Nao,
     +                                  Jroot,Irrepy,Iroot,Irrepx)
C
C Multipole moments <L_J|O|R_I> of one generalized (state-to-state)
C transition density (spin summed, AO basis), stored at
C (Jroot,Irrepy,Iroot,Irrepx). See vee_bd_moments for the operators.
C
C 11/2026 (Ajith Perera): uses vee_bd_moments (antisymmetric AM/MQ
C expansion, MQ antisymmetrized in vprops).
C
      Implicit Double Precision(A-H,O-Z)
      Double Precision Mq,Mq1
      Dimension Dens(Nao,Nao),Scr1(*),Scr2(Nao,Nao)
      Dimension Qm(4,8,4,8,6),Om(4,8,4,8,3),Am(4,8,4,8,3),
     +          Mq(4,8,4,8,9),Dm(4,8,4,8,3)
      Dimension Qm1(6),Om1(3),Am1(3),Mq1(9),Dm1(3)

      Call Vee_bd_moments(Dens,Scr1,Scr2,Nao,Dm1,Qm1,Om1,Am1,Mq1)

      Do Ic = 1, 3
         Dm(Jroot,Irrepy,Iroot,Irrepx,Ic) = Dm1(Ic)
         Am(Jroot,Irrepy,Iroot,Irrepx,Ic) = Am1(Ic)
      Enddo
      Do Ic = 1, 6
         Qm(Jroot,Irrepy,Iroot,Irrepx,Ic) = Qm1(Ic)
      Enddo
      Do Ic = 1, 3
         Om(Jroot,Irrepy,Iroot,Irrepx,Ic) = Om1(Ic)
      Enddo
      Do Ic = 1, 9
         Mq(Jroot,Irrepy,Iroot,Irrepx,Ic) = Mq1(Ic)
      Enddo

      Return
      End
