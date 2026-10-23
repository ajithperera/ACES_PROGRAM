      Subroutine Erd_sphchk(Zscr,Lzscr)
C
C Check of the spherical functions against OED itself: for every
C shell with L >= 2 the spherical-mode overlap of the shell with
C itself, turned into VMOL's harmonics (W^T S W, Erd_vsph), must equal
C Cv^T S' Cv, S' the Cartesian-mode overlap scaled to ACES AOs (Fac)
C and Cv VMOL's harmonics over those AOs; contraction block by block.
C Stops if it does not.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2, Mxpk = 5000,
     &           Mxf = 30*Mxyz)
      Dimension Zscr(Lzscr), Alp(Mxpk), Cc(Mxpk), Cv(Mxyz,Mxyz),
     &          W(2*Lmax+1,2*Lmax+1), Fac(Mxyz), Sc(Mxf*Mxf)
      Integer List(2), Icb(Mxpk), Ice(Mxpk), Iparv(3,Mxyz)
      Logical Lsph, Lcrt

      Lcrt = .False.
      Lsph = .True.
      Errmx = 0.0D0
      Liscr = Lzscr/2
      Lz    = Lzscr - Liscr
      Do Is = 1, Nshl
         L = Lshl(Is)
         If (L .Lt. 2) Goto 100
         List(1) = Is
         List(2) = Is
         Call Erd_pack(2,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
         Ia = Ishatm(Is)
         X = Xyzatm(1,Ia)
         Y = Xyzatm(2,Ia)
         Z = Xyzatm(3,Ia)
         Nc = Nctsh(Is)
         Np = Nprsh(Is)
         Nk = (L+1)*(L+2)/2
         Ns = 2*L + 1
         Call Erd_cmpfac(L,Lcrt,Fac,N)
         Call Erd_vsph(L,Cv,Mxyz,W,2*Lmax+1,Iparv)

         Call Oed__memory_ovl_batch(Nalp,Ncc,Nc,Nc,Np,Np,L,L,X,Y,Z,
     &        X,Y,Z,Alp,Cc,Lcrt,Imin,Iopt,Jmin,Jopt)
         If (Imin .Gt. Liscr .Or. Jmin .Gt. Lz) Call Errex
         Call Oed__gener_ovl_batch(Liscr,Lz,Nalp,Ncc,Ncsum,Nc,Nc,Np,Np,
     &        L,L,X,Y,Z,X,Y,Z,Alp,Cc,Icb,Ice,Lcrt,.False.,
     &        Zscr(Lz+1),Nbc,Nfc,Zscr)
         N1 = Nk*Nc
         Do I = 1, N1*N1
            Sc(I) = Zscr(Nfc+I-1)
         Enddo
         Call Oed__memory_ovl_batch(Nalp,Ncc,Nc,Nc,Np,Np,L,L,X,Y,Z,
     &        X,Y,Z,Alp,Cc,Lsph,Imin,Iopt,Jmin,Jopt)
         If (Imin .Gt. Liscr .Or. Jmin .Gt. Lz) Call Errex
         Call Oed__gener_ovl_batch(Liscr,Lz,Nalp,Ncc,Ncsum,Nc,Nc,Np,Np,
     &        L,L,X,Y,Z,X,Y,Z,Alp,Cc,Icb,Ice,Lsph,.False.,
     &        Zscr(Lz+1),Nbs,Nfs,Zscr)
         N2 = Ns*Nc
         Do It = 1, Nc
         Do Iu = 1, Nc
            Do M = 1, Ns
            Do M2 = 1, Ns
               Sum = 0.0D0
               Do K = 1, Nk
               Do K2 = 1, Nk
                  Sum = Sum + Cv(K,M)*Cv(K2,M2)*Fac(K)*Fac(K2)*
     &                  Sc(K + (It-1)*Nk + N1*(K2-1 + (Iu-1)*Nk))
               Enddo
               Enddo
               Ssph = 0.0D0
               Do J = 1, Ns
               Do J2 = 1, Ns
                  Ssph = Ssph + W(J,M)*W(J2,M2)*Zscr(Nfs + J-1 +
     &                   (It-1)*Ns + N2*(J2-1 + (Iu-1)*Ns))
               Enddo
               Enddo
               Errmx = Max(Errmx,Dabs(Sum-Ssph))
            Enddo
            Enddo
         Enddo
         Enddo
 100     Continue
      Enddo

      Write(6,"(2a,1p,e10.2)") " @Erd_sphchk: spherical functions vs",
     &     " OED, largest overlap error", Errmx
      If (Errmx .Gt. 1.0D-10) Then
         Write(6,"(a)") " @Erd_sphchk: Cartesian -> spherical map wrong"
         Call Errex
      Endif

      Return
      End
