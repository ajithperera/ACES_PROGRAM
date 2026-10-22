      Subroutine Pccd_doci(E,W,Ek,Nb,No,Wr,Lwr,Edoci,Ndim,Nit,Rn,
     +                     Iconv)
C
C DOCI (seniority-zero FCI) at the current orbitals from the pair
C parameters of pCCD (Pccd_rsp_par): in the basis of the C(Nb,No)
C closed-shell determinants S (No doubly occupied of Nb orbitals),
C   <S|H|S>  = sum_{p in S} E(p) + sum_{p<q in S} W(p,q),
C   <S'|H|S> = Ek(p,q)  (S' = S - {q} + {p}: one pair moved q -> p),
C the electronic energy (no nuclear repulsion), as pCCD's Eref + Ec.
C Lowest root by Davidson (diagonal preconditioner, restart to the
C current vector after Mxs vectors) from the lowest diagonal element.
C Determinants in colex order, rank = sum_k C(o_k,k), o_1 < ... < o_No
C (0-based orbitals). Wr: Lwr words of scratch, (2*Mxs+3)*Ndim used.
C Iconv = 1: |r| < Rtol; 0: not converged in Mxit; -1: too large (not
C done). Validated against a determinant FCI in the seniority-zero
C space (pilot_doci.py, H2O STO-3G, 9e-14).
C
      Implicit Double Precision(A-H,O-Z)
      Parameter (Mxb=1000, Mxo=64, Mxsp=20, Mxit=300, Rtol=1.0D-7)
      Dimension E(Nb),W(Nb,Nb),Ek(Nb,Nb),Wr(*)
      Dimension Ib(0:Mxb,Mxo),Hs(Mxsp,Mxsp),Hc(Mxsp,Mxsp),Th(Mxsp),
     +          Wk(3*Mxsp)
      Data Zero,One /0.0D0,1.0D0/

      Iconv = -1
      Ndim  = 0
      Nit   = 0
      Rn    = Zero
      Edoci = Zero
      If (Nb .Gt. Mxb .Or. No .Gt. Mxo .Or. No .Lt. 1) Return
C Binomials C(n,k), n < Nb, k <= No; the ranks that occur are < Ndim.
      Do N = 0, Nb
         Bin = Dble(N)
         Ib(N,1) = N
         Do K = 2, No
            Bin = Bin*Dble(N-K+1)/Dble(K)
            If (Bin .Lt. Zero) Bin = Zero
            If (Bin .Gt. 2.0D9) Then
               Ib(N,K) = -1
            Else
               Ib(N,K) = Nint(Bin)
            Endif
         Enddo
      Enddo
      Bin = One
      Do K = 1, No
         Bin = Bin*Dble(Nb-K+1)/Dble(K)
      Enddo
      If (Bin .Gt. 2.0D9) Return
      Ndim = Nint(Bin)
      Mxs  = Min(Mxsp,(Lwr/Ndim - 3)/2)
      If (Mxs .Lt. 3) Return
      Iconv = 0

      Id = 1
      Ix = Id + Ndim
      Is = Ix + Ndim
      Iv = Is + Ndim
      Ia = Iv + Mxs*Ndim
      Call Pccd_doci_dia(E,W,Nb,No,Wr(Id),Ndim)
      Imin = 1
      Do I = 2, Ndim
         If (Wr(Id+I-1) .Lt. Wr(Id+Imin-1)) Imin = I
      Enddo
      Call Dzero(Wr(Iv),Ndim)
      Wr(Iv+Imin-1) = One
      M = 1
      Do It = 1, Mxit
         Nit = It
         Jv = Iv + (M-1)*Ndim
         Ja = Ia + (M-1)*Ndim
         Call Pccd_doci_sig(Ek,Nb,No,Ib,Mxb,Wr(Id),Wr(Jv),Wr(Ja),Ndim)
         Do J = 1, M
            Hs(J,M) = Ddot(Ndim,Wr(Iv+(J-1)*Ndim),1,Wr(Ja),1)
            Hs(M,J) = Hs(J,M)
         Enddo
         Do J = 1, M
            Do I = 1, M
               Hc(I,J) = Hs(I,J)
            Enddo
         Enddo
         Call Dsyev("V","U",M,Hc,Mxsp,Th,Wk,3*Mxsp,Info)
         Edoci = Th(1)
C x = V y, s = (HV) y, r = s - theta x (in Wr(Is) until restart).
         Call Dgemv("N",Ndim,M,One,Wr(Iv),Ndim,Hc(1,1),1,Zero,
     +              Wr(Ix),1)
         Call Dgemv("N",Ndim,M,One,Wr(Ia),Ndim,Hc(1,1),1,Zero,
     +              Wr(Is),1)
         Rn = Zero
         Do I = 0, Ndim-1
            Rn = Rn + (Wr(Is+I) - Edoci*Wr(Ix+I))**2
         Enddo
         Rn = Dsqrt(Rn)
         If (Rn .Lt. Rtol) Then
            Iconv = 1
            Return
         Endif
         If (M .Eq. Mxs) Then
            Call Dcopy(Ndim,Wr(Ix),1,Wr(Iv),1)
            Call Dcopy(Ndim,Wr(Is),1,Wr(Ia),1)
            Hs(1,1) = Edoci
            M = 1
         Endif
C Correction (r)/(theta - D), orthogonalized twice against V.
         Do I = 0, Ndim-1
            Den = Edoci - Wr(Id+I)
            If (Dabs(Den) .Lt. 1.0D-8) Den = Sign(1.0D-8,Den)
            Wr(Ix+I) = (Wr(Is+I) - Edoci*Wr(Ix+I))/Den
         Enddo
         Do Ipass = 1, 2
            Do J = 1, M
               Ov = Ddot(Ndim,Wr(Iv+(J-1)*Ndim),1,Wr(Ix),1)
               Call Daxpy(Ndim,-Ov,Wr(Iv+(J-1)*Ndim),1,Wr(Ix),1)
            Enddo
         Enddo
         Vn = Dnrm2(Ndim,Wr(Ix),1)
         If (Vn .Lt. 1.0D-12) Return
         M = M + 1
         Call Dcopy(Ndim,Wr(Ix),1,Wr(Iv+(M-1)*Ndim),1)
         Call Dscal(Ndim,One/Vn,Wr(Iv+(M-1)*Ndim),1)
      Enddo
      Return
      End

      Subroutine Pccd_doci_dia(E,W,Nb,No,D,Ndim)
C Diagonal of the DOCI Hamiltonian, determinants in colex order.
      Implicit Double Precision(A-H,O-Z)
      Parameter (Mxo=64)
      Dimension E(Nb),W(Nb,Nb),D(Ndim),Iocc(Mxo)
      Do K = 1, No
         Iocc(K) = K - 1
      Enddo
      Do I = 1, Ndim
         Dsum = 0.0D0
         Do K = 1, No
            Ip  = Iocc(K) + 1
            Dsum = Dsum + E(Ip)
            Do L = 1, K-1
               Dsum = Dsum + W(Ip,Iocc(L)+1)
            Enddo
         Enddo
         D(I) = Dsum
         If (I .Lt. Ndim) Call Pccd_doci_nxt(Iocc,No)
      Enddo
      Return
      End

      Subroutine Pccd_doci_nxt(Iocc,No)
C Next combination in colex order (o_1 < ... < o_No, 0-based).
      Dimension Iocc(No)
      J = 1
 10   If (J .Lt. No) Then
         If (Iocc(J)+1 .Eq. Iocc(J+1)) Then
            J = J + 1
            Go To 10
         Endif
      Endif
      Iocc(J) = Iocc(J) + 1
      Do I = 1, J-1
         Iocc(I) = I - 1
      Enddo
      Return
      End

      Subroutine Pccd_doci_sig(Ek,Nb,No,Ib,Ldb,D,X,Y,Ndim)
C Y = H X for the DOCI Hamiltonian: diagonal D, one-pair moves q -> p
C with Ek(p,q); the rank of the moved determinant in O(No).
      Implicit Double Precision(A-H,O-Z)
      Parameter (Mxo=64, Mxb=1000)
      Logical Locc(0:Mxb)
      Dimension Ek(Nb,Nb),Ib(0:Ldb,*),D(Ndim),X(Ndim),Y(Ndim),
     +          Iocc(Mxo)
      Do Ip = 0, Nb-1
         Locc(Ip) = .False.
      Enddo
      Do K = 1, No
         Iocc(K) = K - 1
         Locc(K-1) = .True.
      Enddo
      Do I = 1, Ndim
         Y(I) = D(I)*X(I)
      Enddo
      Do I = 1, Ndim
         Xi = X(I)
         Do Kq = 1, No
            Iq = Iocc(Kq)
            Do Ip = 0, Nb-1
               If (Locc(Ip)) Go To 20
C Rank of S - {q} + {p}.
               Ir  = 0
               K   = 0
               Ins = 0
               Do J = 1, No
                  Io = Iocc(J)
                  If (Io .Eq. Iq) Go To 15
                  If (Ins .Eq. 0 .And. Ip .Lt. Io) Then
                     K   = K + 1
                     Ir  = Ir + Ib(Ip,K)
                     Ins = 1
                  Endif
                  K  = K + 1
                  Ir = Ir + Ib(Io,K)
 15               Continue
               Enddo
               If (Ins .Eq. 0) Ir = Ir + Ib(Ip,No)
               Y(Ir+1) = Y(Ir+1) + Ek(Ip+1,Iq+1)*Xi
 20            Continue
            Enddo
         Enddo
         If (I .Lt. Ndim) Then
            Do K = 1, No
               Locc(Iocc(K)) = .False.
            Enddo
            Call Pccd_doci_nxt(Iocc,No)
            Do K = 1, No
               Locc(Iocc(K)) = .True.
            Enddo
         Endif
      Enddo
      Return
      End
