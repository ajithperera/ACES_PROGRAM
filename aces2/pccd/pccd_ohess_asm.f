










      Subroutine Pccd_ohess_s0el(Ip,Iq,Ir,Is_,Mode,X,A,S,Nbas,Dmax)
C
C Add the density element X (list slots p,q,r,s) to the seniority-zero
C matrices A(p,q) = Gam(pp|qq) and S(p,q) = Gam(pq|qp) + Gam(pq|pq) of
C the 8-fold symmetrized Gam. Mode 1: the element sits at Mulliken
C (pr|qs) (Dirac partner list); Mode 2: at (pq|rs). Dmax collects the
C largest |X| outside that pattern.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension A(Nbas,Nbas),S(Nbas,Nbas)
      Data Half /0.5D0/

      M1 = Ip
      If (Mode .Eq. 1) Then
         M2 = Ir
         M3 = Iq
      Else
         M2 = Iq
         M3 = Ir
      Endif
      M4 = Is_
      If (M1 .Eq. M2 .And. M3 .Eq. M4) Then
         A(M1,M3) = A(M1,M3) + Half*X
         A(M3,M1) = A(M3,M1) + Half*X
      Elseif ((M1 .Eq. M4 .And. M2 .Eq. M3) .Or.
     +        (M1 .Eq. M3 .And. M2 .Eq. M4)) Then
         S(M1,M2) = S(M1,M2) + Half*X
         S(M2,M1) = S(M2,M1) + Half*X
      Else
         Dmax = Max(Dmax,Dabs(X))
      Endif
      Return
      End

      Subroutine Pccd_ohess_scat(Ia,Ib,Ic,Id,X,Ja,Jb,Jc,Jd,A,S,H,Xj,Xk,
     +                           Kp,Nbas,Npar,Ijk)
C
C One integral X = (ab|cd) of a list with Mulliken index spaces
C (Ja Jb|Jc Jd). Each class of permutationally equal integrals is taken
C once, from its member with a >= b when Ja = Jb, c >= d when Jc = Jd
C and pair(a,b) >= pair(c,d) when (Ja,Jb) = (Jc,Jd) (each class lives
C in one list), and every distinct ordered image goes to Pccd_ohess_img.
C Ijk = 1: only the classes Xj, Xk (no H).
C
      Implicit Double Precision(A-H,O-Z)
      Dimension M(4,8)
      Dimension A(*),S(*),H(*),Xj(*),Xk(*),Kp(*)

      If (Ja .Eq. Jb .And. Ia .Lt. Ib) Return
      If (Jc .Eq. Jd .And. Ic .Lt. Id) Return
      If (Ja .Eq. Jc .And. Jb .Eq. Jd) Then
         Iab = Max(Ia,Ib)*(Max(Ia,Ib)-1)/2 + Min(Ia,Ib)
         Icd = Max(Ic,Id)*(Max(Ic,Id)-1)/2 + Min(Ic,Id)
         If (Iab .Lt. Icd) Return
      Endif
C (ab|cd) (ba|cd) (ab|dc) (ba|dc) and the same with the pairs swapped.
      Do K = 1, 4
         I1 = Ia
         I2 = Ib
         I3 = Ic
         I4 = Id
         If (K .Eq. 2 .Or. K .Eq. 4) Then
            I1 = Ib
            I2 = Ia
         Endif
         If (K .Ge. 3) Then
            I3 = Id
            I4 = Ic
         Endif
         M(1,K)   = I1
         M(2,K)   = I2
         M(3,K)   = I3
         M(4,K)   = I4
         M(1,K+4) = I3
         M(2,K+4) = I4
         M(3,K+4) = I1
         M(4,K+4) = I2
      Enddo
      Do K = 1, 8
         Do L = 1, K-1
            If (M(1,L) .Eq. M(1,K) .And. M(2,L) .Eq. M(2,K) .And.
     +          M(3,L) .Eq. M(3,K) .And. M(4,L) .Eq. M(4,K)) Go to 10
         Enddo
         Call Pccd_ohess_img(M(1,K),M(2,K),M(3,K),M(4,K),X,A,S,H,Xj,Xk,
     +                       Kp,Nbas,Npar,Ijk)
 10      Continue
      Enddo
      Return
      End

      Subroutine Pccd_ohess_img(Ip,Iq,Ir,Is_,X,A,S,H,Xj,Xk,Kp,Nbas,
     +                          Npar,Ijk)
C
C The integral G(p,q,r,s) = (pq|rs) = X (one ordered element): the
C integral classes Xj(p,q,y) = (pq|yy), Xk(p,q,y) = (py|qy), and its
C terms in 2Y of the frozen Hessian (Pccd_ohess),
C   Y(p,q,r,s) = 2A(p,r)(qp|sr) + S(p,r)[(qr|sp) + (qs|pr)] + ...,
C added to H through Pccd_ohess_addh.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension A(Nbas,Nbas),S(Nbas,Nbas),Xj(Nbas,Nbas,Nbas),
     +          Xk(Nbas,Nbas,Nbas),H(Npar,Npar),Kp(Nbas,Nbas)
      Data Two,Four /2.0D0,4.0D0/

      If (Ir .Eq. Is_) Xj(Ip,Iq,Ir) = X
      If (Iq .Eq. Is_) Xk(Ip,Ir,Iq) = X
      If (Ijk .Eq. 1) Return
C (qp|sr): q=p', p=q', s=r', r=s'
      Call Pccd_ohess_addh(Iq,Ip,Is_,Ir,Four*A(Iq,Is_)*X,H,Kp,Nbas,
     +                     Npar)
C (qr|sp): q=p', r=q', s=r', p=s'
      Call Pccd_ohess_addh(Is_,Ip,Iq,Ir,Two*S(Is_,Iq)*X,H,Kp,Nbas,
     +                     Npar)
C (qs|pr): q=p', s=q', p=r', r=s'
      Call Pccd_ohess_addh(Ir,Ip,Is_,Iq,Two*S(Ir,Is_)*X,H,Kp,Nbas,
     +                     Npar)
      Return
      End

      Subroutine Pccd_ohess_addh(Iu,Iv,Iw,Ix,V,H,Kp,Nbas,Npar)
C
C V belongs to A4(u,v,w,x); H(pq,rs) = A4(p,q,r,s) - A4(q,p,r,s)
C - A4(p,q,s,r) + A4(q,p,s,r), p<q, r<s. Kp(p,q) = +k (p<q) or -k
C (p>q) for parameter k, 0 if (p,q) is not a parameter.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension H(Npar,Npar),Kp(Nbas,Nbas)
      K = Kp(Iu,Iv)
      L = Kp(Iw,Ix)
      If (K .Eq. 0 .Or. L .Eq. 0) Return
      If ((K .Gt. 0) .Eqv. (L .Gt. 0)) Then
         H(Iabs(K),Iabs(L)) = H(Iabs(K),Iabs(L)) + V
      Else
         H(Iabs(K),Iabs(L)) = H(Iabs(K),Iabs(L)) - V
      Endif
      Return
      End

      Subroutine Pccd_ohess_kp(Irp,Kp,Nbas)
C Kp(p,q) = k, Kp(q,p) = -k for the RFO parameter k of the pair p<q
C (same irrep, row-major order as Pccd_rfo_step); 0 otherwise.
      Dimension Irp(Nbas),Kp(Nbas,Nbas)
      Do Iq = 1, Nbas
         Do Ip = 1, Nbas
            Kp(Ip,Iq) = 0
         Enddo
      Enddo
      K = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               Kp(Ip,Iq) = K
               Kp(Iq,Ip) = -K
            Endif
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_ohess_fin(Hh,F,A,S,Occ,Xj,Xk,Wy,H,Gr,Irp,Nbas,
     +                          Nocc,Npar)
C
C Complete the frozen-amplitude orbital Hessian (formulas in Pccd_ohess)
C from the integral classes Xj(p,q,y) = (pq|yy), Xk(p,q,y) = (py|qy):
C on entry H holds the 2Y integral terms and Hh = MO Fock matrix f; on
C exit Hh = core Hamiltonian h = f - (2J - K), F the generalized Fock
C matrix, Gr(p,q) = F(p,q) - F(q,p) (OBRTGRDA convention), H complete
C and symmetrized. Wy: Nbas**3 scratch.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Hh(Nbas,Nbas),F(Nbas,Nbas),A(Nbas,Nbas),S(Nbas,Nbas)
      Dimension Occ(Nbas),Gr(Nbas,Nbas),Irp(Nbas),H(Npar,Npar)
      Dimension Xj(Nbas,Nbas,Nbas),Xk(Nbas,Nbas,Nbas),Wy(Nbas,Nbas,Nbas)
      Data Zero,Half,One,Two /0.0D0,0.5D0,1.0D0,2.0D0/

      Call Pccd_ohess_hmo(Hh,Xj,Xk,Nbas,Nocc)
C F(p,q) = n(p)h(p,q) + sum_s A(p,s)(qp|ss) + S(p,s)(qs|sp)
      Do Iq = 1, Nbas
         Do Ip = 1, Nbas
            X = Occ(Ip)*Hh(Ip,Iq)
            Do Is_ = 1, Nbas
               X = X + A(Ip,Is_)*Xj(Iq,Ip,Is_) + S(Ip,Is_)*Xk(Iq,Ip,Is_)
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
      Nsq = Nbas*Nbas
      Call Dgemm("N","T",Nbas,Nsq,Nbas,One,S,Nbas,Xk,Nsq,Zero,Wy,Nbas)
      Call Dgemm("N","T",Nbas,Nsq,Nbas,One,A,Nbas,Xj,Nsq,One,Wy,Nbas)

C Kronecker-delta terms of A4(u,v,w,x):
C   d(u,w)[2n(u)h(v,x) + 2Wy(u,v,x)] - d(v,x)[F(u,w) + F(w,u)]
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
            X = Zero
            If (Ip .Eq. Ir) X = X + Two*(Occ(Ip)*Hh(Iq,Is_)
     +                                 + Wy(Ip,Iq,Is_))
            If (Iq .Eq. Is_) X = X - F(Ip,Ir) - F(Ir,Ip)
            If (Iq .Eq. Ir) X = X - Two*(Occ(Iq)*Hh(Ip,Is_)
     +                                 + Wy(Iq,Ip,Is_))
            If (Ip .Eq. Is_) X = X + F(Iq,Ir) + F(Ir,Iq)
            If (Ip .Eq. Is_) X = X - Two*(Occ(Ip)*Hh(Iq,Ir)
     +                                 + Wy(Ip,Iq,Ir))
            If (Iq .Eq. Ir) X = X + F(Ip,Is_) + F(Is_,Ip)
            If (Iq .Eq. Is_) X = X + Two*(Occ(Iq)*Hh(Ip,Ir)
     +                                 + Wy(Iq,Ip,Ir))
            If (Ip .Eq. Ir) X = X - F(Iq,Is_) - F(Is_,Iq)
            H(K,L) = H(K,L) + X
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

      Subroutine Pccd_ohess_hmo(Hh,Xj,Xk,Nbas,Nocc)
C Hh = f - (2J - K): the MO Fock matrix f in, the core Hamiltonian h out
C (Xj(p,q,y) = (pq|yy), Xk(p,q,y) = (py|qy)).
      Implicit Double Precision(A-H,O-Z)
      Dimension Hh(Nbas,Nbas),Xj(Nbas,Nbas,Nbas),Xk(Nbas,Nbas,Nbas)
      Do Iq = 1, Nbas
         Do Ip = 1, Nbas
            X = 0.0D0
            Do J = 1, Nocc
               X = X + 2.0D0*Xj(Ip,Iq,J) - Xk(Ip,Iq,J)
            Enddo
            Hh(Ip,Iq) = Hh(Ip,Iq) - X
         Enddo
      Enddo
      Return
      End
