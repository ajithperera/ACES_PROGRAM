      Subroutine Erg_one(G,Dt,W,N,Zscr,Lzscr,Iscr,Liscr)
C
C One-electron part of the gradient,
C
C    dE/dX += sum D(mu,nu) d[T+V](mu,nu)/dX - sum W(mu,nu) dS(mu,nu)/dX,
C
C over the symmetry unique shell pairs only: each pair is weighted
C with the size of its orbit and the result is symmetrized later
C (Erg_symz). Dt is the total density and W the energy weighted
C density, both over ERD functions. S and T depend on the two basis
C centers only (dX on the second = -dX on the first); V also depends
C on the nuclei, so its derivative is taken for every atom, with the
C basis centers on that atom and its nucleus differentiated together.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Parameter (Mxpk = 5000)
      Dimension G(3,Natm), Dt(N,N), W(N,N), Zscr(Lzscr)
      Integer Iscr(Liscr)
      Integer List(2), Icb(Mxpk), Ice(Mxpk), Id(3,3)
      Dimension Alp(Mxpk), Cc(Mxpk), Xn(Mxatms), Yn(Mxatms),
     &          Zn(Mxatms)
      Logical Screen, Lrep
      Data Screen /.False./
      Data Id /1,0,0, 0,1,0, 0,0,1/

      Do I = 1, Natm
         Xn(I) = Xyzatm(1,I)
         Yn(I) = Xyzatm(2,I)
         Zn(I) = Xyzatm(3,I)
      Enddo

      Do Is1 = 1, Nshl
      Do Is2 = 1, Is1
         List(1) = Is1
         List(2) = Is2
         Call Erg_orbit(List,2,Lrep,Nimg)
         If (.Not. Lrep) Goto 100
         Fac = Dble(Nimg)
         If (Is1 .Ne. Is2) Fac = 2.0D0*Fac
         Call Erd_pack(2,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
         Ia1 = Ishatm(Is1)
         Ia2 = Ishatm(Is2)
         X1 = Xyzatm(1,Ia1)
         Y1 = Xyzatm(2,Ia1)
         Z1 = Xyzatm(3,Ia1)
         X2 = Xyzatm(1,Ia2)
         Y2 = Xyzatm(2,Ia2)
         Z2 = Xyzatm(3,Ia2)
         Nc1 = Nctsh(Is1)
         Nc2 = Nctsh(Is2)
         Np1 = Nprsh(Is1)
         Np2 = Nprsh(Is2)
         L1  = Lshl(Is1)
         L2  = Lshl(Is2)
         Nf1 = Nfnsh(Is1)
         Nf2 = Nfnsh(Is2)
         I0  = Ioffn(Is1)
         J0  = Ioffn(Is2)

C Overlap and kinetic energy: differentiate the first center.

         If (Ia1 .Ne. Ia2) Then
         Do K = 1, 3
         Do Iop = 1, 2
            If (Iop .Eq. 1) Then
               Call Oed__memory_ovl_derv_batch(Nalp,Ncc,Nc1,Nc2,
     &              Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,Alp,Cc,Spher,
     &              Imin,Iopt,Jmin,Jopt)
            Else
               Call Oed__memory_kin_derv_batch(Nalp,Ncc,Nc1,Nc2,
     &              Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,Alp,Cc,Spher,
     &              Imin,Iopt,Jmin,Jopt)
            Endif
            Call Erg_scrchk("Erg_one",Imin,Liscr,Jmin,Lzscr)
            Imax = Max(Imin,Min(Iopt,Liscr))
            Jmax = Max(Jmin,Min(Jopt,Lzscr))
            If (Iop .Eq. 1) Then
               Call Oed__gener_ovl_derv_batch(Imax,Jmax,Nalp,Ncc,
     &              Ncsum,Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
            Else
               Call Oed__gener_kin_derv_batch(Imax,Jmax,Nalp,Ncc,
     &              Ncsum,Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Id(1,K),Id(2,K),Id(3,K),0,0,0,
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
            Endif
            If (Nbatch .Gt. 0) Then
               Sum = 0.0D0
               Do J = 1, Nf2
               Do I = 1, Nf1
                  X = Zscr(Nfirst + (J-1)*Nf1 + I - 1)
                  If (Iop .Eq. 1) Then
                     Sum = Sum - W(I0+I,J0+J)*X
                  Else
                     Sum = Sum + Dt(I0+I,J0+J)*X
                  Endif
               Enddo
               Enddo
               G(K,Ia1) = G(K,Ia1) + Fac*Sum
               G(K,Ia2) = G(K,Ia2) - Fac*Sum
            Endif
         Enddo
         Enddo
         Endif

C Nuclear attraction: every atom.

         Do Ic = 1, Natm
         Do K = 1, 3
            Kd1 = 0
            If (Ia1 .Eq. Ic) Kd1 = 1
            Kd2 = 0
            If (Ia2 .Eq. Ic) Kd2 = 1
            Call Oed__memory_nai_derv_batch(Nalp,Ncc,Nc1,Nc2,
     &           Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &           Natm,Xn,Yn,Zn,Ic,
     &           Kd1*Id(1,K),Kd1*Id(2,K),Kd1*Id(3,K),
     &           Kd2*Id(1,K),Kd2*Id(2,K),Kd2*Id(3,K),
     &           Id(1,K),Id(2,K),Id(3,K),Alp,Cc,Spher,
     &           Imin,Iopt,Jmin,Jopt)
            Call Erg_scrchk("Erg_one",Imin,Liscr,Jmin,Lzscr)
            Imax = Max(Imin,Min(Iopt,Liscr))
            Jmax = Max(Jmin,Min(Jopt,Lzscr))
            Call Oed__gener_nai_derv_batch(Imax,Jmax,Nalp,Ncc,
     &           Ncsum,Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &           Natm,Xn,Yn,Zn,Zatm,Ic,
     &           Kd1*Id(1,K),Kd1*Id(2,K),Kd1*Id(3,K),
     &           Kd2*Id(1,K),Kd2*Id(2,K),Kd2*Id(3,K),
     &           Id(1,K),Id(2,K),Id(3,K),
     &           Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &           Nbatch,Nfirst,Zscr)
            If (Nbatch .Gt. 0) Then
               Sum = 0.0D0
               Do J = 1, Nf2
               Do I = 1, Nf1
                  Sum = Sum + Dt(I0+I,J0+J)*
     &                  Zscr(Nfirst + (J-1)*Nf1 + I - 1)
               Enddo
               Enddo
               G(K,Ic) = G(K,Ic) + Fac*Sum
            Endif
         Enddo
         Enddo
  100    Continue
      Enddo
      Enddo

      Return
      End
