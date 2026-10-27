      Subroutine Erg_two(G,Ha,Hb,Xa,Xb,N,Zscr,Lzscr,Iscr,Liscr,
     &                   Gam,Lgam,Lsep,Cns,Bx,Ng,Irg,Ncol,Icol0,Lcol,
     &                   Ipx,Npr,Ioffg,Gso,W1,W2,Lw)
C
C Two-electron part of the gradient,
C
C    dE/dX += 1/2 sum Gam(mu,nu,la,si) d(mu nu|la si)/dX.
C
C The separable part of Gam (Lsep) is built from the per spin SCF
C densities Hs and correlation corrections Xs (Erg_dens), linear in X
C as in xvdint's Intexp (H, X: sums over spin):
C
C    H(mu,nu) H(la,si) + H(mu,nu) X(la,si) + X(mu,nu) H(la,si)
C    - 1/2 sum_spin [ Hs(mu,la) Hs(nu,si) + Hs(mu,si) Hs(nu,la)
C                   + Hs(mu,la) Xs(nu,si) + Xs(mu,la) Hs(nu,si)
C                   + Hs(mu,si) Xs(nu,la) + Xs(mu,si) Hs(nu,la) ];
C
C for an SCF gradient X = 0. Cns times the nonseparable part from the
C MOINTS lists (Erg_gam, Erg_gnsp) is added when Cns is not zero.
C The sum runs over the symmetry unique canonical shell quartets
C only. Each quartet carries its permutational degeneracy and its
C orbit size; the result is symmetrized later (Erg_symz). The
C derivative is taken with respect to every atom that carries a
C center of the quartet but the last one: all centers on that atom
C are differentiated together, and the last atom gets minus the sum
C (translational invariance). One-center quartets are skipped.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Parameter (Mxpk = 5000)
      Dimension G(3,Natm), Ha(N,N), Hb(N,N), Xa(N,N), Xb(N,N),
     &          Zscr(Lzscr), Gam(Lgam), Bx(N,Ng), Gso(*), W1(Lw),
     &          W2(Lw)
      Integer Iscr(Liscr), Irg(Ng), Ncol(*), Icol0(*), Lcol(*),
     &        Ipx(Ng,Ng), Npr(8), Ioffg(8)
      Integer List(4), Icb(Mxpk), Ice(Mxpk), Iat(4), Iatd(4),
     &        Ider(3,4)
      Dimension Alp(Mxpk), Cc(Mxpk), X(4), Y(4), Z(4), Gq(3)
      Logical Screen, Lrep, Lsep
      Data Screen /.True./

      Do Is1 = 1, Nshl
      Do Is2 = 1, Is1
      Ij = Is1*(Is1-1)/2 + Is2
      Do Is3 = 1, Is1
      Do Is4 = 1, Is3
      Kl = Is3*(Is3-1)/2 + Is4
      If (Kl .Gt. Ij) Goto 100
         List(1) = Is1
         List(2) = Is2
         List(3) = Is3
         List(4) = Is4

C Distinct atoms of the quartet.

         Nd = 0
         Do I = 1, 4
            Iat(I) = Ishatm(List(I))
            Lnew = 1
            Do J = 1, Nd
               If (Iatd(J) .Eq. Iat(I)) Lnew = 0
            Enddo
            If (Lnew .Eq. 1) Then
               Nd = Nd + 1
               Iatd(Nd) = Iat(I)
            Endif
         Enddo
         If (Nd .Eq. 1) Goto 100

         Call Erg_orbit(List,4,Lrep,Nimg)
         If (.Not. Lrep) Goto 100
         Fac = 0.5D0*Dble(Nimg)
         If (Is1 .Ne. Is2) Fac = 2.0D0*Fac
         If (Is3 .Ne. Is4) Fac = 2.0D0*Fac
         If (Ij  .Ne. Kl ) Fac = 2.0D0*Fac

         N1 = Nfnsh(Is1)
         N2 = Nfnsh(Is2)
         N3 = Nfnsh(Is3)
         N4 = Nfnsh(Is4)
         If (N1*N2*N3*N4 .Gt. Lgam) Then
            Write(6,"(a)") " @Erg_two: Gamma buffer too small"
            Call Errex
         Endif
         I0 = Ioffn(Is1)
         J0 = Ioffn(Is2)
         K0 = Ioffn(Is3)
         L0 = Ioffn(Is4)
         Ind = 0
         Do L = L0+1, L0+N4
         Do K = K0+1, K0+N3
         Do J = J0+1, J0+N2
         Do I = I0+1, I0+N1
            Ind = Ind + 1
            If (Lsep) Then
               Hij = Ha(I,J) + Hb(I,J)
               Hkl = Ha(K,L) + Hb(K,L)
               Gam(Ind) = Hij*Hkl + Hij*(Xa(K,L) + Xb(K,L))
     &                  + (Xa(I,J) + Xb(I,J))*Hkl - 0.5D0*(
     &              Ha(I,K)*Ha(J,L) + Ha(I,L)*Ha(J,K)
     &            + Ha(I,K)*Xa(J,L) + Xa(I,K)*Ha(J,L)
     &            + Ha(I,L)*Xa(J,K) + Xa(I,L)*Ha(J,K)
     &            + Hb(I,K)*Hb(J,L) + Hb(I,L)*Hb(J,K)
     &            + Hb(I,K)*Xb(J,L) + Xb(I,K)*Hb(J,L)
     &            + Hb(I,L)*Xb(J,K) + Xb(I,L)*Hb(J,K))
            Else
               Gam(Ind) = 0.0D0
            Endif
         Enddo
         Enddo
         Enddo
         Enddo
         If (Cns .Ne. 0.0D0) Call Erg_gnsp(List,Cns,Gam,Bx,N,Ng,Irg,
     &                    Ncol,Icol0,Lcol,Ipx,Npr,Ioffg,Gso,W1,W2,Lw)

         Call Erd_pack(4,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
         Do I = 1, 4
            X(I) = Xyzatm(1,Iat(I))
            Y(I) = Xyzatm(2,Iat(I))
            Z(I) = Xyzatm(3,Iat(I))
         Enddo

         Do Ia = 1, Nd-1
            Do Kc = 1, 3
               Do I = 1, 4
                  Do K = 1, 3
                     Ider(K,I) = 0
                  Enddo
                  If (Iat(I) .Eq. Iatd(Ia)) Ider(Kc,I) = 1
               Enddo
               Call Erd__memory_eri_derv_batch(Nalp,Ncc,
     &              Nctsh(Is1),Nctsh(Is2),Nctsh(Is3),Nctsh(Is4),
     &              Nprsh(Is1),Nprsh(Is2),Nprsh(Is3),Nprsh(Is4),
     &              Lshl(Is1),Lshl(Is2),Lshl(Is3),Lshl(Is4),
     &              X(1),Y(1),Z(1),X(2),Y(2),Z(2),
     &              X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &              Ider(1,1),Ider(2,1),Ider(3,1),
     &              Ider(1,2),Ider(2,2),Ider(3,2),
     &              Ider(1,3),Ider(2,3),Ider(3,3),
     &              Ider(1,4),Ider(2,4),Ider(3,4),
     &              Alp,Cc,Spher,Imin,Iopt,Jmin,Jopt)
               Call Erg_scrchk("Erg_two",Imin,Liscr,Jmin,Lzscr)
               Imax = Max(Imin,Min(Iopt,Liscr))
               Jmax = Max(Jmin,Min(Jopt,Lzscr))
               Call Erd__gener_eri_derv_batch(Imax,Jmax,Nalp,Ncc,Ncsum,
     &              Nctsh(Is1),Nctsh(Is2),Nctsh(Is3),Nctsh(Is4),
     &              Nprsh(Is1),Nprsh(Is2),Nprsh(Is3),Nprsh(Is4),
     &              Lshl(Is1),Lshl(Is2),Lshl(Is3),Lshl(Is4),
     &              X(1),Y(1),Z(1),X(2),Y(2),Z(2),
     &              X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &              Ider(1,1),Ider(2,1),Ider(3,1),
     &              Ider(1,2),Ider(2,2),Ider(3,2),
     &              Ider(1,3),Ider(2,3),Ider(3,3),
     &              Ider(1,4),Ider(2,4),Ider(3,4),
     &              Alp,Cc,Icb,Ice,Spher,Screen,Iscr,
     &              Nbatch,Nfirst,Zscr)
               Gq(Kc) = 0.0D0
               If (Nbatch .Gt. 0) Gq(Kc) = Ddot(N1*N2*N3*N4,Gam,1,
     &                                         Zscr(Nfirst),1)
            Enddo
            Do Kc = 1, 3
               G(Kc,Iatd(Ia)) = G(Kc,Iatd(Ia)) + Fac*Gq(Kc)
               G(Kc,Iatd(Nd)) = G(Kc,Iatd(Nd)) - Fac*Gq(Kc)
            Enddo
         Enddo
  100 Continue
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
