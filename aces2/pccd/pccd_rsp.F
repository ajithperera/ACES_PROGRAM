C Exact orbital Hessian-vector products of OO-pCCD (RHF) by amplitude
C response, for the OO_METHOD=RFO stability check.
C
C pCCD is coupled cluster on the seniority-zero part of the Hamiltonian,
C which depends on the MO integrals only through
C   e(p) = 2h(p,p) + (pp|pp),  W(p,q) = 4(pp|qq) - 2(pq|qp),
C   K(p,q) = (pq|pq).
C With pair amplitudes t(i,a) (D: diagonal pair energies),
C   r(i,a) = K(i,a) + D(i,a)t(i,a) + sum_j.ne.i K(i,j)t(j,a)
C          + sum_b.ne.a K(a,b)t(i,b) + (t K^T t)(i,a)
C          - 2t(i,a)(s(i)+s(a)) + 2K(i,a)t(i,a)**2,
C   s(i) = sum_b K(i,b)t(i,b), s(a) = sum_j K(j,a)t(j,a),
C   E = E_ref + sum K(i,a)t(i,a),  L = E + sum lam(i,a) r(i,a).
C Frozen pairs (dropped MOs, PCCDNFZO/NFZV) carry r = t (t = 0).
C With x the RFO rotation parameters and b a direction,
C   H b = L_xx b + d/dx [L_t(P).t1 + r(P).l1],
C   A t1 = -r(t; dP/dx b),  A^T l1 = -(L_tx b + L_tt t1),  A = dr/dt,
C L_xx b from the analytic frozen Hessian (OHESFRZ, Pccd_ohess). The
C relaxation part is the orbital gradient of a functional sum G.P whose
C "densities" G come from Pccd_rsp_gphi; dP/dx b and the gradient need
C h, (pq|yy) and (py|qy) only (O(n**3), saved by Pccd_ohess).
C Python pilot (pccd_oo_survey_1016/pilot_resp*.py, 10/2026): pair
C energy = ACES pCCD to 2e-11; H b vs the full central-difference
C Hessian of ACES gradients: 1e-6 (N2 cc-pVDZ 1.2 A), 1e-5 (2.52 A,
C saddle -0.05635 vs -0.05636).
C
C Workspace layout (Pccd_rsp_init) in Common /Pccdrsp/.

      Subroutine Pccd_rsp_hv(Hf,Wr,Npar,X,Y)
C Y = H X, H the exact orbital Hessian (frozen part Hf + relaxation).
      Implicit Double Precision(A-H,O-Z)
      Dimension Hf(Npar,Npar),Wr(*),X(Npar),Y(Npar)
      Data Zero,One /0.0D0,1.0D0/
      Call Dgemv("N",Npar,Npar,One,Hf,Npar,X,1,Zero,Y,1)
      Call Pccd_rsp_rlx(Wr,Npar,X,Y)
      Return
      End

      Subroutine Pccd_rsp_rlx(Wr,Npar,X,Y)
C Y = Y + (amplitude relaxation part of the exact Hessian) X.
      Implicit Double Precision(A-H,O-Z)
      Dimension Wr(*),X(Npar),Y(Npar)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Common /Pccdrsi/ Irpr(1000)
      Data Zero,Half,One /0.0D0,0.5D0,1.0D0/
      Nsq = Nb*Nb
      Np  = Nb + 2*Nsq
C dP/dx X (the parameter derivatives along X)
      Call Pccd_rsp_bmat(X,Wr(Lbm),Irpr,Nb)
      Call Pccd_rsp_d1(Wr(Lbm),Wr(Lh),Wr(Ljx),Wr(Lkx),Wr(Lde),Wr(Ldw),
     +                 Wr(Ldk),Nb)
C t1 = -A^-1 r(t; d1)   (r is linear in P)
      Call Pccd_rsp_res(Wr(Lt),Wr(Lde),Wr(Ldw),Wr(Ldk),Wr(Lmk),
     +                  Wr(Lt1))
      Do K = 1, Nov
         If (Wr(Lmk+K-1) .Eq. Zero) Wr(Lt1+K-1) = Zero
      Enddo
      Call Dscal(Nov,-One,Wr(Lt1),1)
      Call Dgetrs("N",Nov,1,Wr(La),Nov,Wr(Lpv),Wr(Lt1),Nov,Info)
C w = K(d1)ov + A(t; d1)^T lam + [A(t+t1) - A(t)]^T lam
      Call Pccd_rsp_kov(Wr(Ldk),Wr(Lvw),Nb,No)
      Call Pccd_rsp_vjp(Wr(Llm),Wr(Lt),Wr(Lde),Wr(Ldw),Wr(Ldk),
     +                  Wr(Lmk),Wr(Lr1),1)
      Call Daxpy(Nov,One,Wr(Lr1),1,Wr(Lvw),1)
      Call Dcopy(Nov,Wr(Lt),1,Wr(Lr2),1)
      Call Daxpy(Nov,One,Wr(Lt1),1,Wr(Lr2),1)
      Call Pccd_rsp_vjp(Wr(Llm),Wr(Lr2),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),
     +                  Wr(Lr1),0)
      Call Daxpy(Nov,One,Wr(Lr1),1,Wr(Lvw),1)
      Call Pccd_rsp_vjp(Wr(Llm),Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),
     +                  Wr(Lr1),0)
      Call Daxpy(Nov,-One,Wr(Lr1),1,Wr(Lvw),1)
      Do K = 1, Nov
         Wr(Ll1+K-1) = -Wr(Lvw+K-1)
      Enddo
      Call Dgetrs("T",Nov,1,Wr(La),Nov,Wr(Lpv),Wr(Ll1),Nov,Info)
C G = {K ov: t1} + [Gphi(lam, t+t1) - Gphi(lam, t-t1)]/2 + Gphi(l1, t)
      Call Pccd_rsp_gphi(Wr(Llm),Wr(Lr2),Wr(Lmk),Wr(Lg1),Nb,No)
      Call Dcopy(Nov,Wr(Lt),1,Wr(Lr2),1)
      Call Daxpy(Nov,-One,Wr(Lt1),1,Wr(Lr2),1)
      Call Pccd_rsp_gphi(Wr(Llm),Wr(Lr2),Wr(Lmk),Wr(Lg2),Nb,No)
      Call Daxpy(Np,-One,Wr(Lg2),1,Wr(Lg1),1)
      Call Dscal(Np,Half,Wr(Lg1),1)
      Call Pccd_rsp_gphi(Wr(Ll1),Wr(Lt),Wr(Lmk),Wr(Lg2),Nb,No)
      Call Daxpy(Np,One,Wr(Lg2),1,Wr(Lg1),1)
      Call Pccd_rsp_addk(Wr(Lt1),Wr(Lg1),Nb,No)
      Call Pccd_rsp_zgr(Wr(Lg1),Wr(Lh),Wr(Ljx),Wr(Lkx),Wr(Lz),Irpr,Nb,
     +                  Npar,Wr(Ixs))
      Call Daxpy(Npar,One,Wr(Ixs),1,Y,1)
      Return
      End

      Subroutine Pccd_rsp_par(H,Xj,Xk,E,W,Ek,Nb,Iful)
C Pair parameters from h, Xj(p,q,y) = (pq|yy), Xk(p,q,y) = (py|qy)
C (Iful = 1); with Iful = 0 the arguments are already derivatives and
C only the linear combination is formed (not used).
      Implicit Double Precision(A-H,O-Z)
      Dimension H(Nb,Nb),Xj(Nb,Nb,Nb),Xk(Nb,Nb,Nb),E(Nb),W(Nb,Nb),
     +          Ek(Nb,Nb)
      Do Ip = 1, Nb
         E(Ip) = 2.0D0*H(Ip,Ip) + Xj(Ip,Ip,Ip)
         Do Iq = 1, Nb
            Ek(Ip,Iq) = Xk(Ip,Ip,Iq)
            W(Ip,Iq)  = 4.0D0*Xj(Ip,Ip,Iq) - 2.0D0*Ek(Ip,Iq)
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_mask(Am,Irp,Nfzo,Nfzv)
C Am(i,a) = 0 for frozen pairs, 1 otherwise.
      Implicit Double Precision(A-H,O-Z)
      Dimension Irp(*),Nfzo(8),Nfzv(8),Lo(8),Lv(8),Nvi(8)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Common /Pccdrsi/ Irpr(1000)
      Dimension Am(No,Nv)
      Do Ip = 1, Nb
         Irpr(Ip) = Irp(Ip)
      Enddo
      Do Irr = 1, 8
         Lo(Irr)  = 0
         Lv(Irr)  = 0
         Nvi(Irr) = 0
      Enddo
      Do Ia = 1, Nv
         Nvi(Irp(No+Ia)) = Nvi(Irp(No+Ia)) + 1
      Enddo
      Do Ia = 1, Nv
         Ir = Irp(No+Ia)
         Lv(Ir) = Lv(Ir) + 1
         Do I = 1, No
            Am(I,Ia) = 1.0D0
         Enddo
         If (Lv(Ir) .Gt. Nvi(Ir)-Nfzv(Ir)) Then
            Do I = 1, No
               Am(I,Ia) = 0.0D0
            Enddo
         Endif
      Enddo
      Do I = 1, No
         Ir = Irp(I)
         Lo(Ir) = Lo(Ir) + 1
         If (Lo(Ir) .Le. Nfzo(Ir)) Then
            Do Ia = 1, Nv
               Am(I,Ia) = 0.0D0
            Enddo
         Endif
      Enddo
      Return
      End

      Subroutine Pccd_rsp_kov(Ek,X,Nb,No)
C X(i,a) = K(i,a).
      Implicit Double Precision(A-H,O-Z)
      Dimension Ek(Nb,Nb),X(No,*)
      Do Ia = 1, Nb-No
         Do I = 1, No
            X(I,Ia) = Ek(I,No+Ia)
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_res(T,E,W,Ek,Am,R)
C Pair residual r(t; P) (linear in P = (E, W, Ek)); frozen pairs r = t.
C The sums over j.ne.i, b.ne.a and (t K^T t) by DGEMM (O(o v**2)), the
C excluded diagonal terms K(i,i)t(i,a), K(a,a)t(i,a) taken off after.
      Implicit Double Precision(A-H,O-Z)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Dimension T(No,Nv),E(Nb),W(Nb,Nb),Ek(Nb,Nb),Am(No,Nv),R(No,Nv)
      Dimension Si(1000),Sa(1000),Xm(1000*1000)
      Save Xm
      Data Zero,One /0.0D0,1.0D0/
      Do I = 1, No
         Si(I) = Zero
         Do Ia = 1, Nv
            Si(I) = Si(I) + Ek(I,No+Ia)*T(I,Ia)
         Enddo
      Enddo
      Do Ia = 1, Nv
         Sa(Ia) = Zero
         Do I = 1, No
            Sa(Ia) = Sa(Ia) + Ek(I,No+Ia)*T(I,Ia)
         Enddo
      Enddo
C R = Koo t + t Kvv + (t Kov^T) t
      Call Dgemm("N","N",No,Nv,No,One,Ek,Nb,T,No,Zero,R,No)
      Call Dgemm("N","N",No,Nv,Nv,One,T,No,Ek(No+1,No+1),Nb,One,R,No)
      Call Dgemm("N","T",No,No,Nv,One,T,No,Ek(1,No+1),Nb,Zero,Xm,No)
      Call Dgemm("N","N",No,Nv,No,One,Xm,No,T,No,One,R,No)
      Do Ia = 1, Nv
         Iav = No + Ia
         Do I = 1, No
            X = R(I,Ia) - (Ek(I,I) + Ek(Iav,Iav))*T(I,Ia)
     +        + Ek(I,Iav) + Pccd_rsp_dia(I,Iav,E,W,Nb,No)*T(I,Ia)
     +        - 2.0D0*T(I,Ia)*(Si(I) + Sa(Ia))
     +        + 2.0D0*Ek(I,Iav)*T(I,Ia)*T(I,Ia)
            R(I,Ia) = Am(I,Ia)*X + (One - Am(I,Ia))*T(I,Ia)
         Enddo
      Enddo
      Return
      End

      Double Precision Function Pccd_rsp_dia(I,Iav,E,W,Nb,No)
C D(i,a) = e(a) - e(i) + sum_{k occ, k.ne.i} (W(a,k) - W(i,k)).
      Implicit Double Precision(A-H,O-Z)
      Dimension E(Nb),W(Nb,Nb)
      X = E(Iav) - E(I)
      Do K = 1, No
         If (K .Ne. I) X = X + W(Iav,K) - W(I,K)
      Enddo
      Pccd_rsp_dia = X
      Return
      End

      Subroutine Pccd_rsp_jac(T,E,W,Ek,Am,A,Dt,Rp,Rm)
C A(ia,jb) = dr(ia)/dt(jb): r is quadratic in t, so the central
C difference is exact up to round-off.
      Implicit Double Precision(A-H,O-Z)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Dimension T(Nov),E(Nb),W(Nb,Nb),Ek(Nb,Nb),Am(Nov),A(Nov,Nov),
     +          Dt(Nov),Rp(Nov),Rm(Nov)
      Parameter (Del=1.0D-3)
      Do K = 1, Nov
         Call Dcopy(Nov,T,1,Dt,1)
         Dt(K) = T(K) + Del
         Call Pccd_rsp_res(Dt,E,W,Ek,Am,Rp)
         Dt(K) = T(K) - Del
         Call Pccd_rsp_res(Dt,E,W,Ek,Am,Rm)
         Do L = 1, Nov
            A(L,K) = (Rp(L) - Rm(L))/(2.0D0*Del)
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_vjp(Rl,T,E,W,Ek,Am,V,Ilin)
C V = A(t; P)^T lam (frozen pairs: r = t). Ilin = 1: P is a derivative
C (dP/dx b); the frozen-pair identity block has no P dependence and is
C left out.
      Implicit Double Precision(A-H,O-Z)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Dimension Rl(No,Nv),T(No,Nv),E(Nb),W(Nb,Nb),Ek(Nb,Nb),Am(No,Nv),
     +          V(No,Nv)
      Dimension Si(1000),Sa(1000),Ui(1000),Ua(1000)
      Do I = 1, No
         Si(I) = 0.0D0
         Ui(I) = 0.0D0
         Do Ia = 1, Nv
            Si(I) = Si(I) + Ek(I,No+Ia)*T(I,Ia)
            Ui(I) = Ui(I) + Am(I,Ia)*Rl(I,Ia)*T(I,Ia)
         Enddo
      Enddo
      Do Ia = 1, Nv
         Sa(Ia) = 0.0D0
         Ua(Ia) = 0.0D0
         Do I = 1, No
            Sa(Ia) = Sa(Ia) + Ek(I,No+Ia)*T(I,Ia)
            Ua(Ia) = Ua(Ia) + Am(I,Ia)*Rl(I,Ia)*T(I,Ia)
         Enddo
      Enddo
      Do Ib = 1, Nv
         Ibv = No + Ib
         Do J = 1, No
            Rjb = Am(J,Ib)*Rl(J,Ib)
            X = Rjb*Pccd_rsp_dia(J,Ibv,E,W,Nb,No)
C K(i,j) t(j,a), K(a,b) t(i,b)
            Do I = 1, No
               If (I .Ne. J) X = X + Ek(I,J)*Am(I,Ib)*Rl(I,Ib)
            Enddo
            Do Ia = 1, Nv
               If (Ia .Ne. Ib) X = X + Ek(No+Ia,Ibv)*Am(J,Ia)*Rl(J,Ia)
            Enddo
C (t K^T t): lam (t^T K) and (K t^T lam)
            Y = 0.0D0
            Do Ia = 1, Nv
               Z = 0.0D0
               Do K = 1, No
                  Z = Z + T(K,Ia)*Ek(K,Ibv)
               Enddo
               Y = Y + Am(J,Ia)*Rl(J,Ia)*Z
            Enddo
            Do I = 1, No
               Z = 0.0D0
               Do Ic = 1, Nv
                  Z = Z + T(I,Ic)*Ek(J,No+Ic)
               Enddo
               Y = Y + Z*Am(I,Ib)*Rl(I,Ib)
            Enddo
            X = X + Y - 2.0D0*Rjb*(Si(J) + Sa(Ib))
     +            - 2.0D0*Ek(J,Ibv)*(Ui(J) + Ua(Ib))
     +            + 4.0D0*Rjb*Ek(J,Ibv)*T(J,Ib)
            If (Ilin .Eq. 0) X = X + (1.0D0 - Am(J,Ib))*Rl(J,Ib)
            V(J,Ib) = X
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_gen(T,Am,G,Nb,No)
C G = d/dP of E = E_ref + sum K(i,a) t(i,a); G = (Ge(Nb), Gw, Gk).
      Implicit Double Precision(A-H,O-Z)
      Dimension T(No,*),Am(No,*),G(*)
      Nsq = Nb*Nb
      Call Dzero(G,Nb+2*Nsq)
      Iw = Nb
      Ik = Nb + Nsq
      Do I = 1, No
         G(I) = 1.0D0
         Do J = 1, No
            If (J .Ne. I) G(Iw+(J-1)*Nb+I) = 0.5D0
         Enddo
      Enddo
      Do Ia = 1, Nb-No
         Do I = 1, No
            G(Ik+(No+Ia-1)*Nb+I) = T(I,Ia)
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_addk(T,G,Nb,No)
C G(K ov) += T.
      Implicit Double Precision(A-H,O-Z)
      Dimension T(No,*),G(*)
      Ik = Nb + Nb*Nb
      Do Ia = 1, Nb-No
         Do I = 1, No
            G(Ik+(No+Ia-1)*Nb+I) = G(Ik+(No+Ia-1)*Nb+I) + T(I,Ia)
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_gphi(Rl,T,Am,G,Nb,No)
C G = d/dP of Phi = sum_ia lam(i,a) r(i,a)(t; P) over unfrozen pairs,
C G = (Ge(Nb), Gw(Nb,Nb), Gk(Nb,Nb)), entries as used by Pccd_rsp_res.
      Implicit Double Precision(A-H,O-Z)
      Dimension Rl(No,*),T(No,*),Am(No,*),G(*)
      Dimension Ui(1000),Ua(1000)
      Nv  = Nb - No
      Nsq = Nb*Nb
      Iw  = Nb
      Ik  = Nb + Nsq
      Call Dzero(G,Nb+2*Nsq)
      Do I = 1, No
         Ui(I) = 0.0D0
      Enddo
      Do Ia = 1, Nv
         Ua(Ia) = 0.0D0
         Do I = 1, No
            U = Am(I,Ia)*Rl(I,Ia)*T(I,Ia)
            Ui(I)  = Ui(I) + U
            Ua(Ia) = Ua(Ia) + U
         Enddo
      Enddo
C D(i,a): e(a) - e(i) + sum_{k.ne.i} W(a,k) - W(i,k)
      Do Ia = 1, Nv
         G(No+Ia) = G(No+Ia) + Ua(Ia)
      Enddo
      Do I = 1, No
         G(I) = G(I) - Ui(I)
      Enddo
      Do K = 1, No
         Do Ia = 1, Nv
            U = Am(K,Ia)*Rl(K,Ia)*T(K,Ia)
            G(Iw+(K-1)*Nb+No+Ia) = G(Iw+(K-1)*Nb+No+Ia) + Ua(Ia) - U
         Enddo
         Do I = 1, No
            If (I .Ne. K) G(Iw+(K-1)*Nb+I) = G(Iw+(K-1)*Nb+I) - Ui(I)
         Enddo
      Enddo
C K(i,a): lam, (t lam^T t), -2 Ui t(i,c), -2 Ua t(k,a), 2 lam t**2
      Do Ic = 1, Nv
         Do K = 1, No
            X = Am(K,Ic)*Rl(K,Ic)*(1.0D0 + 2.0D0*T(K,Ic)*T(K,Ic))
     +        - 2.0D0*Ui(K)*T(K,Ic) - 2.0D0*Ua(Ic)*T(K,Ic)
            Y = 0.0D0
            Do Ia = 1, Nv
               Z = 0.0D0
               Do I = 1, No
                  Z = Z + Am(I,Ia)*Rl(I,Ia)*T(I,Ic)
               Enddo
               Y = Y + T(K,Ia)*Z
            Enddo
            G(Ik+(No+Ic-1)*Nb+K) = G(Ik+(No+Ic-1)*Nb+K) + X + Y
         Enddo
      Enddo
C K(i,j) t(j,a) -> Gk(i,j) = sum_a lam(i,a) t(j,a), i.ne.j
      Do J = 1, No
         Do I = 1, No
            If (I .Ne. J) Then
               X = 0.0D0
               Do Ia = 1, Nv
                  X = X + Am(I,Ia)*Rl(I,Ia)*T(J,Ia)
               Enddo
               G(Ik+(J-1)*Nb+I) = G(Ik+(J-1)*Nb+I) + X
            Endif
         Enddo
      Enddo
C K(b,a) t(i,b) -> Gk(b,a) = sum_i lam(i,a) t(i,b), a.ne.b
      Do Ia = 1, Nv
         Do Ib = 1, Nv
            If (Ia .Ne. Ib) Then
               X = 0.0D0
               Do I = 1, No
                  X = X + Am(I,Ia)*Rl(I,Ia)*T(I,Ib)
               Enddo
               G(Ik+(No+Ia-1)*Nb+No+Ib) = G(Ik+(No+Ia-1)*Nb+No+Ib) + X
            Endif
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_bmat(X,B,Irp,Nb)
C Antisymmetric B(p,q) = x(pq), B(q,p) = -x(pq), pairs p<q of one irrep.
      Implicit Double Precision(A-H,O-Z)
      Dimension X(*),B(Nb,Nb),Irp(Nb)
      Call Dzero(B,Nb*Nb)
      K = 0
      Do Ip = 1, Nb
         Do Iq = Ip+1, Nb
            If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               B(Ip,Iq) = X(K)
               B(Iq,Ip) = -X(K)
            Endif
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_d1(B,H,Xj,Xk,De,Dw,Dk,Nb)
C Parameter derivatives along x = e B (C = C exp(eB)):
C   de(r)  = 4 sum_t h(r,t)B(t,r) + 4 sum_t B(t,r)(tr|rr)
C   dJ(x,y) = 2 sum_t B(t,x)(tx|yy) + 2 sum_t B(t,y)(ty|xx)
C   dK(x,y) = 2 sum_t B(t,x)(ty|xy) + 2 sum_t B(t,y)(tx|yx)
C   dW = 4 dJ - 2 dK; (tx|yy) = Xj(t,x,y), (ty|xy) = Xk(t,x,y).
      Implicit Double Precision(A-H,O-Z)
      Dimension B(Nb,Nb),H(Nb,Nb),Xj(Nb,Nb,Nb),Xk(Nb,Nb,Nb),De(Nb),
     +          Dw(Nb,Nb),Dk(Nb,Nb)
      Do Ir = 1, Nb
         X = 0.0D0
         Do It = 1, Nb
            X = X + B(It,Ir)*(H(Ir,It) + Xj(It,Ir,Ir))
         Enddo
         De(Ir) = 4.0D0*X
      Enddo
      Do Iy = 1, Nb
         Do Ix = 1, Nb
            Xjx = 0.0D0
            Xkx = 0.0D0
            Do It = 1, Nb
               Xjx = Xjx + B(It,Ix)*Xj(It,Ix,Iy) + B(It,Iy)*Xj(It,Iy,Ix)
               Xkx = Xkx + B(It,Ix)*Xk(It,Ix,Iy) + B(It,Iy)*Xk(It,Iy,Ix)
            Enddo
            Dk(Ix,Iy) = 2.0D0*Xkx
            Dw(Ix,Iy) = 8.0D0*Xjx - 4.0D0*Xkx
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_zgr(G,H,Xj,Xk,Z,Irp,Nb,Npar,Y)
C Orbital gradient (pairs p<q of one irrep) of F = sum G.P:
C Y(pq) = Z(p,q) - Z(q,p), Z(p,q) = occurrences of orbital q in the
C integrals of F replaced by orbital p.
      Implicit Double Precision(A-H,O-Z)
      Dimension G(*),H(Nb,Nb),Xj(Nb,Nb,Nb),Xk(Nb,Nb,Nb),Z(Nb,Nb),
     +          Irp(Nb),Y(Npar)
      Nsq = Nb*Nb
      Iw  = Nb
      Ik  = Nb + Nsq
      Do Iq = 1, Nb
         Do Ip = 1, Nb
            X = G(Iq)*(4.0D0*H(Ip,Iq) + 4.0D0*Xj(Ip,Iq,Iq))
            Do Iy = 1, Nb
               Gwq = G(Iw+(Iy-1)*Nb+Iq)
               Gwy = G(Iw+(Iq-1)*Nb+Iy)
               Gcq = G(Ik+(Iy-1)*Nb+Iq) - 2.0D0*Gwq
               Gcy = G(Ik+(Iq-1)*Nb+Iy) - 2.0D0*Gwy
               X = X + 8.0D0*(Gwq + Gwy)*Xj(Ip,Iq,Iy)
     +               + 2.0D0*(Gcq + Gcy)*Xk(Ip,Iq,Iy)
            Enddo
            Z(Ip,Iq) = X
         Enddo
      Enddo
      K = 0
      Do Ip = 1, Nb
         Do Iq = Ip+1, Nb
            If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               Y(K) = Z(Ip,Iq) - Z(Iq,Ip)
            Endif
         Enddo
      Enddo
      Return
      End

      Subroutine Pccd_rsp_stab(Hf,Gv,Wr,Lwr,X,Ev,Irp,Nbas,Nocc,Npar,
     +                         Iact,Alam,Vec,Init)
C Stability check with the exact Hessian (Pccd_rsp_hv) at a stationary
C point, all within this xpccd run: the lowest Nroot eigenpairs by
C Davidson. Exact invariance zero modes (|lambda| < Xzero) are skipped.
C Iact = 1: minimum (Alam = lowest nonzero eigenvalue); 2: negative
C eigenvalue Alam < Xlamth with eigenvector Vec; -1: not available
C (use the finite-difference check). X(Npar,Nroot), Ev(Nroot): scratch.
C Init = 0: Pccd_rsp_init was already done at this point.
      Implicit Double Precision(A-H,O-Z)
      Dimension Hf(Npar,Npar),Gv(Npar),Wr(Lwr),X(Npar,*),Ev(*),
     +          Irp(Nbas),Vec(Npar)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Logical Test
      Parameter (Xzero=1.0D-8, Xlamth=-1.0D-4, Ndns=80)
      Iact = -1
C Davidson: roots, residual tolerance, subspace size; debug override in
C pccd_rsp_dav.in. The subspace is the largest up to Mxsp that fits the
C workspace (Pccd_rsp_size reserves 40). Tests 10/2026 (N2 cc-pVDZ
C 1.4-2.52 A, N2 6-31G, C2H4 cc-pVTZ sym off): subspace 40 -> 200 gave
C the same roots with 2.5-8x fewer products (C2H4 3044 -> 528); Rtol
C 1e-6 lost exact zero modes (skipped eigenvalues) and 1e-5 converged
C C2H4 to a wrong lowest root after 9 products.
      Mxrt = 8
      Rtol = 1.0D-8
      Mxsp = 200
      Inquire(File="pccd_rsp_dav.in",Exist=Test)
      If (Test) Then
         Open(Unit=45,File="pccd_rsp_dav.in",Status="Old")
         Read(45,*,End=5,Err=5) Mxrt,Rtol,Mxsp
 5       Close(45)
         Mxrt = Min(Mxrt,30)
         Write(6,"(a,i3,a,1p,E9.2,a,i4)") " @-RFO stability Davidson:",
     +        Mxrt," roots, Rtol",Rtol,", subspace",Mxsp
      Endif
      Alam = 0.0D0
      If (Nbas .Gt. 1000) Return
      Iok = 1
      If (Init .Eq. 1) Call Pccd_rsp_init(Wr,Lwr,Nbas,Nocc,Npar,Irp,Gv,
     +                                    Iok,1)
      If (Iok .Ne. 1) Then
         Write(6,"(2a)") " @-RFO exact Hessian not available; finite",
     +                   "-difference stability check."
         Return
      Endif
C Debug (file pccd_rsp_test.in, Npar values b): H b split into the
C frozen part and the relaxation, to pccd_rsp_test.out.
      Inquire(File="pccd_rsp_test.in",Exist=Test)
      If (Test) Then
         Open(Unit=45,File="pccd_rsp_test.in",Status="Old")
         Read(45,*) (Vec(K),K=1,Npar)
         Close(45)
         Call Pccd_rsp_hv(Hf,Wr,Npar,Vec,X(1,1))
         Call Dgemv("N",Npar,Npar,1.0D0,Hf,Npar,Vec,1,0.0D0,X(1,2),1)
         Open(Unit=45,File="pccd_rsp_test.out",Status="Unknown")
         Write(45,"(1p,3E25.16)") (X(K,2),X(K,1)-X(K,2),X(K,1),K=1,Npar)
         Close(45)
      Endif
      Nroot = Min(Mxrt,Npar)
 7    Continue
      If (2*Npar*Mxsp + 2*Npar + 2*Mxsp*Mxsp + 65*Mxsp .Gt. Ldvl .And.
     +    Mxsp .Gt. 2*Nroot+1) Then
         Mxsp = Mxsp - 10
         Go to 7
      Endif
      If (Npar .Le. Ndns) Then
C Small: the exact Hessian column by column (Npar products), dense
C diagonalization (fits in the Davidson workspace for Npar <= 2*40).
         Call Dcopy(Npar*Npar,Hf,1,Wr(Ldv),1)
         Do K = 1, Npar
            Call Dzero(Vec,Npar)
            Vec(K) = 1.0D0
            Call Pccd_rsp_rlx(Wr,Npar,Vec,Wr(Ldv+(K-1)*Npar))
         Enddo
         Call Pccd_rsp_sym(Wr(Ldv),Npar)
         Lwk = Ldvl - Npar*Npar
         Call Dsyev("V","U",Npar,Wr(Ldv),Npar,Wr(Ldv+Npar*Npar),
     +              Wr(Ldv+Npar*Npar+Npar),Lwk-Npar,Info)
         Iconv = 0
         If (Info .Eq. 0) Iconv = 1
         Nmv = Npar
         Do Ir = 1, Nroot
            Ev(Ir) = Wr(Ldv+Npar*Npar+Ir-1)
            Call Dcopy(Npar,Wr(Ldv+(Ir-1)*Npar),1,X(1,Ir),1)
         Enddo
      Else
      Call Pccd_rfo_dav(Hf,Gv,Wr,Npar,2,Nroot,X,Ev,Rtol,Wr(Ldv),Ldvl,
     +                  Mxsp,Iconv,Nmv)
      Endif
      If (Iconv .Ne. 1) Then
         Write(6,"(2a,i6,a)") " @-RFO exact-Hessian Davidson not",
     +        " converged after",Nmv," Hv; finite-difference check."
         Return
      Endif
      Nzero = 0
      Iphys = 0
      Do Ir = 1, Nroot
         If (Dabs(Ev(Ir)) .Lt. Xzero) Then
            Nzero = Nzero + 1
         Elseif (Iphys .Eq. 0) Then
            Iphys = Ir
         Endif
      Enddo
      Write(6,"(a,i5,a,i2,a,1p,8E11.3)") " @-RFO exact Hessian (",Nmv,
     +     " Hv,",Nzero," zero): lowest roots",(Ev(Ir),Ir=1,Nroot)
      If (Ev(1) .Lt. Xlamth) Then
         Iact = 2
         Alam = Ev(1)
         Call Dcopy(Npar,X(1,1),1,Vec,1)
      Else
         Iact = 1
         If (Iphys .Gt. 0) Alam = Ev(Iphys)
         Write(6,"(a,1p,E12.4,a,i2,a)")
     +      " @-RFO lowest Hessian eigenvalue",Alam,": minimum (",Nzero,
     +      " zero modes skipped)."
      Endif
      Return
      End

      Subroutine Pccd_rsp_sym(A,N)
C A = (A + A^T)/2.
      Implicit Double Precision(A-H,O-Z)
      Dimension A(N,N)
      Do J = 1, N
         Do I = J+1, N
            X = 0.5D0*(A(I,J) + A(J,I))
            A(I,J) = X
            A(J,I) = X
         Enddo
      Enddo
      Return
      End
