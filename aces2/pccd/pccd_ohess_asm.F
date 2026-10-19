      Subroutine Pccd_ohess_putg(G,T,Is1,Is2,Is3,Is4,Nocc,Nvrt,Nbas)
C
C Place the dense Dirac list T(p,q,r,s) = <pq|rs> = (pr|qs) (index
C spaces Is = 1 occupied, 2 virtual) into G(p,q,r,s) = (pq|rs) over all
C MOs, with its eight permutational images (real orbitals).
C
      Implicit Double Precision(A-H,O-Z)
      Dimension G(Nbas,Nbas,Nbas,Nbas),N(4),Io(4),Is(4)
      Dimension T(*)

      Is(1) = Is1
      Is(2) = Is2
      Is(3) = Is3
      Is(4) = Is4
      Do K = 1, 4
         N(K)  = Nocc
         Io(K) = 0
         If (Is(K) .Eq. 2) Then
            N(K)  = Nvrt
            Io(K) = Nocc
         Endif
      Enddo
      K = 0
      Do I4 = 1, N(4)
      Is_ = Io(4) + I4
      Do I3 = 1, N(3)
      Ir = Io(3) + I3
      Do I2 = 1, N(2)
      Iq = Io(2) + I2
      Do I1 = 1, N(1)
         Ip = Io(1) + I1
         K = K + 1
         X = T(K)
         G(Ip,Ir,Iq,Is_) = X
         G(Ir,Ip,Iq,Is_) = X
         G(Ip,Ir,Is_,Iq) = X
         G(Ir,Ip,Is_,Iq) = X
         G(Iq,Is_,Ip,Ir) = X
         G(Is_,Iq,Ip,Ir) = X
         G(Iq,Is_,Ir,Ip) = X
         G(Is_,Iq,Ir,Ip) = X
      Enddo
      Enddo
      Enddo
      Enddo
      Return
      End

      Subroutine Pccd_ohess_sen0(T,Is1,Is2,Is3,Is4,Mode,W,A,S,Nocc,
     +                           Nvrt,Nbas,Dmax)
C
C Add W times the dense density list T (slots p,q,r,s; spaces Is) to
C the seniority-zero matrices A(p,q) = Gam(pp|qq) and S(p,q) =
C Gam(pq|qp) + Gam(pq|pq) of the 8-fold symmetrized Gam. Mode 1: the
C element sits at Mulliken (pr|qs) (Dirac partner list); Mode 2: at
C (pq|rs). Dmax collects the largest |W*T| outside that pattern.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension A(Nbas,Nbas),S(Nbas,Nbas),N(4),Io(4),Is(4),M(4)
      Dimension T(*)
      Data Half /0.5D0/

      Is(1) = Is1
      Is(2) = Is2
      Is(3) = Is3
      Is(4) = Is4
      Do K = 1, 4
         N(K)  = Nocc
         Io(K) = 0
         If (Is(K) .Eq. 2) Then
            N(K)  = Nvrt
            Io(K) = Nocc
         Endif
      Enddo
      K = 0
      Do I4 = 1, N(4)
      Do I3 = 1, N(3)
      Do I2 = 1, N(2)
      Do I1 = 1, N(1)
         K = K + 1
         X = W*T(K)
         If (X .Eq. 0.0D0) Go to 10
         M(1) = Io(1) + I1
         If (Mode .Eq. 1) Then
            M(2) = Io(3) + I3
            M(3) = Io(2) + I2
         Else
            M(2) = Io(2) + I2
            M(3) = Io(3) + I3
         Endif
         M(4) = Io(4) + I4
         If (M(1) .Eq. M(2) .And. M(3) .Eq. M(4)) Then
            A(M(1),M(3)) = A(M(1),M(3)) + Half*X
            A(M(3),M(1)) = A(M(3),M(1)) + Half*X
         Elseif ((M(1) .Eq. M(4) .And. M(2) .Eq. M(3)) .Or.
     +           (M(1) .Eq. M(3) .And. M(2) .Eq. M(4))) Then
            S(M(1),M(2)) = S(M(1),M(2)) + Half*X
            S(M(2),M(1)) = S(M(2),M(1)) + Half*X
         Else
            Dmax = Max(Dmax,Dabs(X))
         Endif
 10      Continue
      Enddo
      Enddo
      Enddo
      Enddo
      Return
      End

      Subroutine Pccd_ohess_asm(G,Hh,F,A,S,Occ,Wy,H,Gr,Irp,Nbas,Nocc,
     +                          Npar)
C
C Assemble the frozen-amplitude orbital Hessian H(Npar,Npar) from the
C seniority-zero densities (Occ, A, S) and the integrals G = (pq|rs)
C (formulas in Pccd_ohess). On entry Hh = MO Fock matrix f, on exit
C the core Hamiltonian h = f - (2J - K). Gr(p,q) = F(p,q) - F(q,p),
C the orbital gradient in the OBRTGRDA convention.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension G(Nbas,Nbas,Nbas,Nbas),Hh(Nbas,Nbas),F(Nbas,Nbas)
      Dimension A(Nbas,Nbas),S(Nbas,Nbas),Occ(Nbas),Gr(Nbas,Nbas)
      Dimension Wy(Nbas,Nbas,Nbas),H(Npar,Npar),Irp(Nbas)
      Data Zero,Half,Two /0.0D0,0.5D0,2.0D0/

      Do Iq = 1, Nbas
         Do Ip = 1, Nbas
            X = Zero
            Do J = 1, Nocc
               X = X + Two*G(Ip,Iq,J,J) - G(Ip,J,J,Iq)
            Enddo
            Hh(Ip,Iq) = Hh(Ip,Iq) - X
         Enddo
      Enddo

      Do Iq = 1, Nbas
         Do Ip = 1, Nbas
            X = Occ(Ip)*Hh(Ip,Iq)
            Do Is_ = 1, Nbas
               X = X + A(Ip,Is_)*G(Iq,Ip,Is_,Is_)
     +               + S(Ip,Is_)*G(Iq,Is_,Is_,Ip)
            Enddo
            F(Ip,Iq) = X
         Enddo
      Enddo
      Do Iq = 1, Nbas
         Do Ip = 1, Nbas
            Gr(Ip,Iq) = F(Ip,Iq) - F(Iq,Ip)
         Enddo
      Enddo

C Wy(p,q,s) = sum_m S(p,m)(qm|sm) + A(p,m)(qs|mm)
      Do Is_ = 1, Nbas
         Do Iq = 1, Nbas
            Do Ip = 1, Nbas
               X = Zero
               Do M = 1, Nbas
                  X = X + S(Ip,M)*G(Iq,M,Is_,M) + A(Ip,M)*G(Iq,Is_,M,M)
               Enddo
               Wy(Ip,Iq,Is_) = X
            Enddo
         Enddo
      Enddo

      K = 0
      Do Ip = 1, Nbas
      Do Iq = Ip+1, Nbas
      If (Irp(Ip) .Eq. Irp(Iq)) Then
         K = K + 1
         L = 0
         Do Ir = 1, Nbas
         Do Is_ = Ir+1, Nbas
         If (Irp(Ir) .Eq. Irp(Is_)) Then
            L = L + 1
            H(K,L) = Pccd_ohess_a4(Ip,Iq,Ir,Is_,G,Hh,F,A,S,Occ,Wy,Nbas)
     +             - Pccd_ohess_a4(Iq,Ip,Ir,Is_,G,Hh,F,A,S,Occ,Wy,Nbas)
     +             - Pccd_ohess_a4(Ip,Iq,Is_,Ir,G,Hh,F,A,S,Occ,Wy,Nbas)
     +             + Pccd_ohess_a4(Iq,Ip,Is_,Ir,G,Hh,F,A,S,Occ,Wy,Nbas)
         Endif
         Enddo
         Enddo
      Endif
      Enddo
      Enddo
      Do L = 1, Npar
         Do K = L+1, Npar
            X = Half*(H(K,L) + H(L,K))
            H(K,L) = X
            H(L,K) = X
         Enddo
      Enddo
      Return
      End

      Double Precision Function Pccd_ohess_a4(Ip,Iq,Ir,Is_,G,Hh,F,A,
     +                                        S,Occ,Wy,Nbas)
C
C A4(p,q,r,s) = 2n(p)d(p,r)h(q,s) - (F(p,r)+F(r,p))d(q,s) + 2Y(p,q,r,s)
C Y(p,q,r,s)  = 2A(p,r)(qp|sr) + S(p,r)[(qr|sp)+(qs|pr)] + d(p,r)Wy(p,q,s)
C
      Implicit Double Precision(A-H,O-Z)
      Dimension G(Nbas,Nbas,Nbas,Nbas),Hh(Nbas,Nbas),F(Nbas,Nbas)
      Dimension A(Nbas,Nbas),S(Nbas,Nbas),Occ(Nbas)
      Dimension Wy(Nbas,Nbas,Nbas)
      Data Two /2.0D0/

      Y = Two*A(Ip,Ir)*G(Iq,Ip,Is_,Ir)
     +  + S(Ip,Ir)*(G(Iq,Ir,Is_,Ip) + G(Iq,Is_,Ip,Ir))
      X = Two*Y
      If (Ip .Eq. Ir) X = X + Two*Occ(Ip)*Hh(Iq,Is_)
     +                      + Two*Wy(Ip,Iq,Is_)
      If (Iq .Eq. Is_) X = X - (F(Ip,Ir) + F(Ir,Ip))
      Pccd_ohess_a4 = X
      Return
      End
