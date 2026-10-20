      Subroutine Pccd_ohess(Work,Maxcor,Imode,Iok)
C
C Analytic frozen-amplitude orbital Hessian of pCCD (RHF), in the RFO
C parametrization (Pccd_rfo_step): C(new) = C*exp(K), K(p,q) = x(pq),
C pairs p<q of the same irrep, row-major order. Written to JOBARC
C OHESFRZ (Npar*Npar); the orbital gradient from the same densities
C goes to OHESGRAD (Nbas*Nbas, F(p,q)-F(q,p), compare OBRTGRDA).
C Iok = 1 on success, 0 if skipped (memory).
C
C The total spin-summed pCCD densities are seniority zero: the 1-RDM is
C diagonal (n(p)) and the 2-RDM (Mulliken, Gam(pq|rs)) has only the
C elements A(p,q) = Gam(pp|qq) and S(p,q) = Gam(pq|qp) + Gam(pq|pq)
C (p.ne.q). With h the core Hamiltonian and g(pq|rs) the MO integrals
C   F(p,q)    = n(p) h(p,q) + sum_s A(p,s)(qp|ss) + sum_r S(p,r)(qr|rp)
C   Y(p,q,r,s)= 2A(p,r)(qp|sr) + S(p,r)[(qr|sp)+(qs|pr)]
C             + d(p,r)[sum_m S(p,m)(qm|sm) + sum_m A(p,m)(qs|mm)]
C   A4        = 2n(p)d(p,r)h(q,s) - (F(p,r)+F(r,p))d(q,s) + 2Y(p,q,r,s)
C   E2(pq,rs) = (1-P(pq))(1-P(rs)) A4(p,q,r,s)
C (Helgaker, Jorgensen, Olsen 10.8 for a seniority-zero density). The
C densities: n = 2 n_ref + 2 D(DENSOO/DENSVV); Gam = Gam_ref(n_ref) +
C cross terms(n_ref, D) + 2 G(113) + 2 G(133) + 4 G(116) + 4 G(123)
C + 4 G(125), each Gamma list placed through its partner integral
C (Dirac <pq|rs> = (pr|qs); 123 is (ab|ij) - (ai|bj)). Must run before
C Pccd_form_htau_pq, which reorders lists 123/125 in place.
C Validated (10/2026, Python pilot from these lists) against the exact
C frozen-amplitude finite-difference Hessian, H2O STO-3G: 3.8e-5 (FD
C noise).
C The integral lists are read in column batches; every integral goes
C straight into H (the 2Y terms) and into the classes (pq|yy), (py|qy)
C (Nbas**3 each), which give h, F and the Kronecker-delta terms. Memory:
C Npar**2 + 3 Nbas**3 + one list column (no Nbas**4 array).
C Imode 0: densities from DENSOO/DENSVV and the Gamma lists (above).
C Imode 1: integral classes and h only (RSP_JX, RSP_KX, RSP_HMO; no
C          densities, no H), for the pair-form pCCD solver (Pccd_pair).
C Imode 2: as 0 with n, A and S from the pair-form solver (JOBARC
C          PAIR_OCC, PAIR_DA, PAIR_DS) instead of the ACES densities.
C
      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry
      Dimension Work(Maxcor)
      Common/Symm/Symmetry
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
      Data One,Zero /1.0D0,0.0D0/

      Iok = 0
      Nocc = 0
      Nvrt = 0
      Do Irr = 1, Nirrep
         Nocc = Nocc + Pop(Irr,1)
         Nvrt = Nvrt + Vrt(Irr,1)
      Enddo
      Nbas = Nocc + Nvrt
      Nsq  = Nbas*Nbas
      N3   = Nsq*Nbas

C Longest column (distribution) of the lists read.
      Mxrow = 0
      Do Irr = 1, Nirrep
         Do Lt = 9, 18
            Mxrow = Max(Mxrow,Irpdpd(Irr,Lt))
         Enddo
      Enddo

      Irpo = 1
      Iirp = Irpo + 8
      Ikp  = Iirp + Nbas
      Ic   = Ikp  + Nsq
      Ifao = Ic   + Nsq
      Iscr = Ifao + Nsq
      Ihh  = Iscr + Nsq
      Iff  = Ihh  + Nsq
      Iam  = Iff  + Nsq
      Ism  = Iam  + Nsq
      Ion  = Ism  + Nsq
      Ixj  = Ion  + Nbas
      Ixk  = Ixj  + N3
      Iwy  = Ixk  + N3
      Ihs  = Iwy  + N3
      Call Pccd_rfo_irp(Work(Irpo),Work(Iirp),Nbas,Symmetry)
      Call Pccd_rfo_npar(Work(Iirp),Nbas,Npar)
      Ibuf = Ihs + Npar*Npar
      If (Imode .Eq. 1) Ibuf = Ihs
      Iend = Ibuf + Mxrow
      If (Iend .Gt. Maxcor) Then
         Write(6,"(2a,i12,a,i12)") " @-Pccd_ohess: analytic orbital",
     +        " Hessian skipped, needs",Iend," words, have",Maxcor
         Return
      Endif
C List buffer: the rest of the memory, at most one irrep block.
      Lbuf = Min(Maxcor-Ibuf+1,Mxrow*Mxrow)

C Occupations and the reference/cross parts of A and S, then the
C seniority-zero elements of the correlation Gamma lists. Dmax: the
C largest element of those lists outside the seniority-zero pattern.
      Lb = Ibuf
      Dmax = Zero
      If (Imode .Eq. 1) Go to 20
      If (Imode .Eq. 2) Then
         Call Getrec(20,"JOBARC","PAIR_OCC",Nbas*Iintfp,Work(Ion))
         Call Getrec(20,"JOBARC","PAIR_DA ",Nsq*Iintfp,Work(Iam))
         Call Getrec(20,"JOBARC","PAIR_DS ",Nsq*Iintfp,Work(Ism))
         Go to 20
      Endif
      Call Pccd_ohess_occ(Work(Ion),Work(Iscr),Work(Iam),Work(Ism),
     +                    Work(Iirp),Nbas,Nocc)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,113,1,1,1,1,Nocc,1,2.0D0,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,133,2,2,2,2,Nocc,1,2.0D0,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,116,2,2,1,1,Nocc,1,4.0D0,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,123,2,1,2,1,Nocc,1,4.0D0,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,123,2,1,2,1,Nocc,2,-4.0D0,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,125,2,1,2,1,Nocc,1,4.0D0,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      If (Dmax .Gt. 1.0D-10) Write(6,"(2a,1p,E10.3)")
     +   " @-Pccd_ohess: WARNING, 2-RDM not seniority zero, max",
     +   " dropped element",Dmax

C Spatial MO integrals (pq|rs), MOs occupied by irrep then virtual by
C irrep. Space 1 = occupied, 2 = virtual.
 20   Continue
      Mlst = 3
      If (Imode .Eq. 1) Mlst = 4
      Call Pccd_ohess_kp(Work(Iirp),Work(Ikp),Nbas)
      If (Imode .Ne. 1) Call Dzero(Work(Ihs),Npar*Npar)
      Call Dzero(Work(Ixj),2*N3)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,13,1,1,1,1,Nocc,Mlst,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,233,2,2,2,2,Nocc,Mlst,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,16,2,2,1,1,Nocc,Mlst,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,10,1,1,1,2,Nocc,Mlst,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,30,2,2,2,1,Nocc,Mlst,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,25,2,1,2,1,Nocc,Mlst,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)

C h = f - (2J - K), f = C^T FOCKA C at the current orbitals.
      Call Getrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,Work(Ic))
      Call Getrec(20,"JOBARC","FOCKA   ",Nsq*Iintfp,Work(Ifao))
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(Ifao),Nbas,Work(Ic),
     +           Nbas,Zero,Work(Iscr),Nbas)
      Call Dgemm("T","N",Nbas,Nbas,Nbas,One,Work(Ic),Nbas,Work(Iscr),
     +           Nbas,Zero,Work(Ihh),Nbas)
      If (Imode .Eq. 1) Then
         Call Pccd_ohess_hmo(Work(Ihh),Work(Ixj),Work(Ixk),Nbas,Nocc)
         Go to 30
      Endif
      Call Pccd_ohess_fin(Work(Ihh),Work(Iff),Work(Iam),Work(Ism),
     +                    Work(Ion),Work(Ixj),Work(Ixk),Work(Iwy),
     +                    Work(Ihs),Work(Iscr),Work(Iirp),Nbas,Nocc,
     +                    Npar)
      Call Putrec(20,"JOBARC","OHESFRZ ",Npar*Npar*Iintfp,Work(Ihs))
      Call Putrec(20,"JOBARC","OHESGRAD",Nsq*Iintfp,Work(Iscr))
C h and the integral classes of the exact-Hessian response (Pccd_rsp).
 30   Continue
      Call Putrec(20,"JOBARC","RSP_JX  ",N3*Iintfp,Work(Ixj))
      Call Putrec(20,"JOBARC","RSP_KX  ",N3*Iintfp,Work(Ixk))
      Call Putrec(20,"JOBARC","RSP_HMO ",Nsq*Iintfp,Work(Ihh))
      Iok = 1
      Return
      End
