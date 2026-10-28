










      Subroutine Pccd_pair_prm(Wr,Irp)
C Read h, (pq|yy), (py|qy) (RSP_HMO, RSP_JX, RSP_KX of Pccd_ohess) into
C the /Pccdrsp/ workspace, form the pair parameters e, W, K and the
C frozen-pair mask (PCCDNFZO/PCCDNFZV: the first Nfzo occupied and last
C Nfzv virtual MOs of each irrep, occupied/virtual by irrep as Irp).
      Implicit Double Precision(A-H,O-Z)
      Dimension Wr(*),Irp(*),Nfzo(8),Nfzv(8)
      Common /MACHSP/ Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Nsq = Nb*Nb
      Call Getrec(20,"JOBARC","RSP_HMO ",Nsq*Iintfp,Wr(Lh))
      Call Getrec(20,"JOBARC","RSP_JX  ",Nsq*Nb*Iintfp,Wr(Ljx))
      Call Getrec(20,"JOBARC","RSP_KX  ",Nsq*Nb*Iintfp,Wr(Lkx))
      Call Pccd_rsp_par(Wr(Lh),Wr(Ljx),Wr(Lkx),Wr(Le),Wr(Lw),Wr(Lk),
     +                  Nb,1)
      Do Irr = 1, 8
         Nfzo(Irr) = 0
         Nfzv(Irr) = 0
      Enddo
      Call Getrec(0,"JOBARC","PCCDNFZO",Length,Junk)
      If (Length .Gt. 0) Then
         Call Getrec(20,"JOBARC","PCCDNFZO",8,Nfzo)
         Call Getrec(20,"JOBARC","PCCDNFZV",8,Nfzv)
      Endif
      Call Pccd_rsp_mask(Wr(Lmk),Irp,Nfzo,Nfzv)
      Return
      End

      Subroutine Pccd_pair_slv(Wr,Irp,Npar,Rtol,Mxit,Ec,Eref,Alm,Nit,
     +                         Rn,Iconv)
C
C Pair-form pCCD (RHF) in place of the CCD-shaped T2/L2 solvers: with
C the parameters of Pccd_pair_prm and the start t in Wr(Lt),
C   r(t) = 0 by Newton steps (A = dr/dt, Pccd_rsp_jac; a step that does
C            not lower max|r| is halved, at most Mxls times),
C   A^T lam = -K(ov),
C   G = dL/dP (L = E + sum lam r, P = (e, W, K)), Y = orbital gradient
C            of L (Pccd_rfo_do convention, Npar, in Wr(Ixs)),
C   the seniority-zero densities n, A, S (Pccd_pair_den) in Wr(Lg2).
C Ec = sum K t (correlation energy), Eref = sum_i e(i) + 1/2 sum_{i.ne.j}
C W(i,j) (electronic reference energy), Alm = max|lam|, Nit Newton
C steps, Rn = max|r|. Iconv = 1: max|r| < Rtol within Mxit steps and A
C nonsingular. O(Nov**3) per step (dense A); r and the products cost
C O(o v**2).
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Wr(*),Irp(*)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Parameter (Mxls=8)
      Data Zero,Half,One /0.0D0,0.5D0,1.0D0/

      Iconv = 0
      Nsq = Nb*Nb
      Np  = Nb + 2*Nsq
      Do K = 1, Nov
         Wr(Lt+K-1) = Wr(Lt+K-1)*Wr(Lmk+K-1)
      Enddo
      Call Pccd_rsp_res(Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),Wr(Lr1))
      Rn = Pccd_pair_amax(Wr(Lr1),Nov)
      Nit = 0
      Do It = 1, Mxit
         If (Rn .Lt. Rtol) Go to 10
         Nit = It
         Call Pccd_rsp_jac(Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),Wr(La),
     +                     Wr(Lt1),Wr(Lr2),Wr(Lvw))
         Call Dgetrf(Nov,Nov,Wr(La),Nov,Wr(Lpv),Info)
         If (Info .Ne. 0) Return
         Call Dgetrs("N",Nov,1,Wr(La),Nov,Wr(Lpv),Wr(Lr1),Nov,Info)
         Alpha = One
         Do Ils = 1, Mxls
            Call Dcopy(Nov,Wr(Lt),1,Wr(Ll1),1)
            Call Daxpy(Nov,-Alpha,Wr(Lr1),1,Wr(Ll1),1)
            Call Pccd_rsp_res(Wr(Ll1),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),
     +                        Wr(Lr2))
            Rn1 = Pccd_pair_amax(Wr(Lr2),Nov)
            If (Rn1 .Lt. Rn) Go to 5
            Alpha = Half*Alpha
         Enddo
 5       Continue
         Call Dcopy(Nov,Wr(Ll1),1,Wr(Lt),1)
         Call Dcopy(Nov,Wr(Lr2),1,Wr(Lr1),1)
         Rn = Rn1
      Enddo
      If (Rn .Ge. Rtol) Return
 10   Continue

C Lambda at the converged t.
      Call Pccd_rsp_jac(Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),Wr(La),
     +                  Wr(Lt1),Wr(Lr2),Wr(Lvw))
      Call Dgetrf(Nov,Nov,Wr(La),Nov,Wr(Lpv),Info)
      If (Info .Ne. 0) Return
      Call Pccd_rsp_kov(Wr(Lk),Wr(Llm),Nb,No)
      Call Dscal(Nov,-One,Wr(Llm),1)
      Call Dgetrs("T",Nov,1,Wr(La),Nov,Wr(Lpv),Wr(Llm),Nov,Info)
      Alm = Pccd_pair_amax(Wr(Llm),Nov)

C dL/dP, the orbital gradient and the densities.
      Call Pccd_rsp_gen(Wr(Lt),Wr(Lmk),Wr(Lg1),Nb,No)
      Call Pccd_rsp_gphi(Wr(Llm),Wr(Lt),Wr(Lmk),Wr(Lg2),Nb,No)
      Call Daxpy(Np,One,Wr(Lg2),1,Wr(Lg1),1)
      Call Pccd_rsp_zgr(Wr(Lg1),Wr(Lh),Wr(Ljx),Wr(Lkx),Wr(Lz),Irp,Nb,
     +                  Npar,Wr(Ixs))
      Call Pccd_pair_den(Wr(Lg1),Wr(Lg2),Wr(Lg2+Nb),Wr(Lg2+Nb+Nsq),Nb)

C Energies.
      Ec = Zero
      Do Ia = 1, Nv
         Do I = 1, No
            K = (Ia-1)*No + I
            Ec = Ec + Wr(Lmk+K-1)*Wr(Lk+(No+Ia-1)*Nb+I-1)*Wr(Lt+K-1)
         Enddo
      Enddo
      Eref = Zero
      Do I = 1, No
         Eref = Eref + Wr(Le+I-1)
         Do J = 1, No
            If (J .Ne. I) Eref = Eref + Half*Wr(Lw+(J-1)*Nb+I-1)
         Enddo
      Enddo
      Iconv = 1
      Return
      End

      Subroutine Pccd_pair_den(G,Occ,A,S,Nb)
C Seniority-zero densities of Pccd_ohess (E = sum n(p)h(p,p) + 1/2 sum
C [A(p,q)(pp|qq) + S(p,q)(pq|qp)], S(p,p) = 0) from G = dL/dP of the
C pair form, P = (e, W, K): e(p) = 2h(p,p) + (pp|pp), W(p,q) = 4(pp|qq)
C - 2(pq|qp), K(p,q) = (pq|pq).
      Implicit Double Precision(A-H,O-Z)
      Dimension G(*),Occ(Nb),A(Nb,Nb),S(Nb,Nb)
      Nsq = Nb*Nb
      Iw  = Nb
      Ik  = Nb + Nsq
      Do Iq = 1, Nb
         Do Ip = 1, Nb
            If (Ip .Eq. Iq) Then
               A(Ip,Ip) = 2.0D0*G(Ip)
               S(Ip,Ip) = 0.0D0
            Else
               Gw = G(Iw+(Iq-1)*Nb+Ip) + G(Iw+(Ip-1)*Nb+Iq)
               Gk = G(Ik+(Iq-1)*Nb+Ip) + G(Ik+(Ip-1)*Nb+Iq)
               A(Ip,Iq) = 4.0D0*Gw
               S(Ip,Iq) = Gk - 2.0D0*Gw
            Endif
         Enddo
         Occ(Iq) = 2.0D0*G(Iq)
      Enddo
      Return
      End

      Double Precision Function Pccd_pair_amax(X,N)
      Implicit Double Precision(A-H,O-Z)
      Dimension X(N)
      Pccd_pair_amax = 0.0D0
      Do K = 1, N
         Pccd_pair_amax = Max(Pccd_pair_amax,Dabs(X(K)))
      Enddo
      Return
      End

      Subroutine Pccd_pair_grd(Y,Irp,G,Nbas)
C G(p,q) = -Y/4, G(q,p) = Y/4 for the parameters p<q (one irrep), so
C that Gv = -2(G(p,q) - G(q,p)) = Y (Pccd_rfo_do).
      Implicit Double Precision(A-H,O-Z)
      Dimension Y(*),Irp(Nbas),G(Nbas,Nbas)
      Call Dzero(G,Nbas*Nbas)
      K = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               G(Ip,Iq) = -0.25D0*Y(K)
               G(Iq,Ip) =  0.25D0*Y(K)
            Endif
         Enddo
      Enddo
      Return
      End
