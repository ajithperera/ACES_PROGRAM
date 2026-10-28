      Subroutine Erd_ao1e(S,T,V,N,Zscr,Lzscr,Iscr,Liscr)
C
C One-electron AO integrals with OED: overlap S, kinetic energy T and
C nuclear attraction V (all atoms of MOL_ERDOED), as full N x N
C matrices (N = Nfnerd) over the basis functions the SOs are made of
C (VMOL's; see Erd_fntr), in ERD's function order. Within a shell
C pair the OED batch is (f1,f2) with f1 running fastest.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Parameter (Mxpk = 5000)
      Dimension S(N,N), T(N,N), V(N,N), Zscr(Lzscr)
      Integer Iscr(Liscr)
      Integer List(2), Icb(Mxpk), Ice(Mxpk)
      Dimension Alp(Mxpk), Cc(Mxpk), Xn(Mxatms), Yn(Mxatms),
     &          Zn(Mxatms)
      Logical Screen
      Data Screen /.False./

      Do I = 1, Natm
         Xn(I) = Xyzatm(1,I)
         Yn(I) = Xyzatm(2,I)
         Zn(I) = Xyzatm(3,I)
      Enddo

      Do Is1 = 1, Nshl
      Do Is2 = 1, Is1
         List(1) = Is1
         List(2) = Is2
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

         Do Iop = 1, 3
            If (Iop .Eq. 1) Then
               Call Oed__memory_ovl_batch(Nalp,Ncc,Nc1,Nc2,Np1,Np2,
     &              L1,L2,X1,Y1,Z1,X2,Y2,Z2,Alp,Cc,Spher,
     &              Imin,Iopt,Jmin,Jopt)
            Else If (Iop .Eq. 2) Then
               Call Oed__memory_kin_batch(Nalp,Ncc,Nc1,Nc2,Np1,Np2,
     &              L1,L2,X1,Y1,Z1,X2,Y2,Z2,Alp,Cc,Spher,
     &              Imin,Iopt,Jmin,Jopt)
            Else
               Call Oed__memory_nai_batch(Nalp,Ncc,Nc1,Nc2,Np1,Np2,
     &              L1,L2,X1,Y1,Z1,X2,Y2,Z2,Natm,Alp,Cc,Spher,
     &              Imin,Iopt,Jmin,Jopt)
            Endif
            If (Imin .Gt. Liscr .Or. Jmin .Gt. Lzscr) Then
               Write(6,"(a,4i12)") " @Erd_ao1e: OED scratch too small",
     &                             Imin, Liscr, Jmin, Lzscr
               Call Errex
            Endif
            Imax = Min(Iopt,Liscr)
            Imax = Max(Imax,Imin)
            Jmax = Min(Jopt,Lzscr)
            Jmax = Max(Jmax,Jmin)

            If (Iop .Eq. 1) Then
               Call Oed__gener_ovl_batch(Imax,Jmax,Nalp,Ncc,Ncsum,
     &              Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
            Else If (Iop .Eq. 2) Then
               Call Oed__gener_kin_batch(Imax,Jmax,Nalp,Ncc,Ncsum,
     &              Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
            Else
               Call Oed__gener_nai_batch(Imax,Jmax,Nalp,Ncc,Ncsum,
     &              Nc1,Nc2,Np1,Np2,L1,L2,X1,Y1,Z1,X2,Y2,Z2,
     &              Natm,Xn,Yn,Zn,Zatm,
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
            Endif

            Nf1 = Nfnsh(Is1)
            Nf2 = Nfnsh(Is2)
            I0  = Ioffn(Is1)
            J0  = Ioffn(Is2)
            Do J = 1, Nf2
            Do I = 1, Nf1
               If (Nbatch .Gt. 0) Then
                  X = Zscr(Nfirst + (J-1)*Nf1 + I - 1)
               Else
                  X = 0.0D0
               Endif
               If (Iop .Eq. 1) Then
                  S(I0+I,J0+J) = X
                  S(J0+J,I0+I) = X
               Else If (Iop .Eq. 2) Then
                  T(I0+I,J0+J) = X
                  T(J0+J,I0+I) = X
               Else
                  V(I0+I,J0+J) = X
                  V(J0+J,I0+I) = X
               Endif
            Enddo
            Enddo
         Enddo
      Enddo
      Enddo

C ERD functions -> the basis functions of the SOs, rows and columns.

      Do Is = 1, Nshl
         I0 = Ioffn(Is) + 1
         Call Erd_fntr(S(I0,1),1,N,N,Is)
         Call Erd_fntr(T(I0,1),1,N,N,Is)
         Call Erd_fntr(V(I0,1),1,N,N,Is)
         Call Erd_fntr(S(1,I0),N,N,1,Is)
         Call Erd_fntr(T(1,I0),N,N,1,Is)
         Call Erd_fntr(V(1,I0),N,N,1,Is)
      Enddo

      Return
      End
