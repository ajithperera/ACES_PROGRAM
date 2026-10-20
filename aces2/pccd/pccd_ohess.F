      Subroutine Pccd_ohess(Work,Maxcor,Iok)
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
C noise). The integrals are held as a dense Nbas**4 array.
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
      N4   = Nsq*Nsq

C Largest unpacked list (dense, all four index spaces) and the largest
C single-irrep block of the lists read.
      Ltmax = Max(Nocc**4,Nvrt**4,Nocc*Nocc*Nvrt*Nvrt,
     +            Nocc**3*Nvrt,Nvrt**3*Nocc)
      Lbmax = 0
      Do Irr = 1, Nirrep
         Do Lt = 9, 18
            Lbmax = Max(Lbmax,Irpdpd(Irr,Lt))
         Enddo
      Enddo
      Lbmax = Lbmax*Lbmax

      Irpo = 1
      Iirp = Irpo + 8
      Igg  = Iirp + Nbas
      Ic   = Igg  + N4
      Ifao = Ic   + Nsq
      Iscr = Ifao + Nsq
      Ihh  = Iscr + Nsq
      Iff  = Ihh  + Nsq
      Iam  = Iff  + Nsq
      Ism  = Iam  + Nsq
      Ion  = Ism  + Nsq
      Iwy  = Ion  + Nbas
      Itt  = Iwy  + Nsq*Nbas
      Ibuf = Itt  + Ltmax
      Ihs  = Ibuf + Lbmax
      Call Pccd_rfo_irp(Work(Irpo),Work(Iirp),Nbas,Symmetry)
      Call Pccd_rfo_npar(Work(Iirp),Nbas,Npar)
      Iend = Ihs + Npar*Npar
      If (Iend .Gt. Maxcor) Then
         Write(6,"(2a,i12,a,i12)") " @-Pccd_ohess: analytic orbital",
     +        " Hessian skipped, needs",Iend," words, have",Maxcor
         Return
      Endif

C Spatial MO integrals g(p,q,r,s) = (pq|rs), MOs occupied by irrep then
C virtual by irrep. Space 1 = occupied, 2 = virtual.
      Lt = Itt
      Lb = Ibuf
      Call Dzero(Work(Igg),N4)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),13,1,1,1,1,Nocc,Nvrt)
      Call Pccd_ohess_putg(Work(Igg),Work(Lt),1,1,1,1,Nocc,Nvrt,Nbas)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),233,2,2,2,2,Nocc,Nvrt)
      Call Pccd_ohess_putg(Work(Igg),Work(Lt),2,2,2,2,Nocc,Nvrt,Nbas)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),16,2,2,1,1,Nocc,Nvrt)
      Call Pccd_ohess_putg(Work(Igg),Work(Lt),2,2,1,1,Nocc,Nvrt,Nbas)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),10,1,1,1,2,Nocc,Nvrt)
      Call Pccd_ohess_putg(Work(Igg),Work(Lt),1,1,1,2,Nocc,Nvrt,Nbas)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),30,2,2,2,1,Nocc,Nvrt)
      Call Pccd_ohess_putg(Work(Igg),Work(Lt),2,2,2,1,Nocc,Nvrt,Nbas)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),25,2,1,2,1,Nocc,Nvrt)
      Call Pccd_ohess_putg(Work(Igg),Work(Lt),2,1,2,1,Nocc,Nvrt,Nbas)

C Occupations and the reference/cross parts of A and S, then the
C seniority-zero elements of the correlation Gamma lists. Dmax: the
C largest element of those lists outside the seniority-zero pattern.
      Call Pccd_ohess_occ(Work(Ion),Work(Iscr),Work(Iam),Work(Ism),
     +                    Work(Iirp),Nbas,Nocc)
      Dmax = Zero
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),113,1,1,1,1,Nocc,Nvrt)
      Call Pccd_ohess_sen0(Work(Lt),1,1,1,1,1,2.0D0,Work(Iam),
     +                     Work(Ism),Nocc,Nvrt,Nbas,Dmax)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),133,2,2,2,2,Nocc,Nvrt)
      Call Pccd_ohess_sen0(Work(Lt),2,2,2,2,1,2.0D0,Work(Iam),
     +                     Work(Ism),Nocc,Nvrt,Nbas,Dmax)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),116,2,2,1,1,Nocc,Nvrt)
      Call Pccd_ohess_sen0(Work(Lt),2,2,1,1,1,4.0D0,Work(Iam),
     +                     Work(Ism),Nocc,Nvrt,Nbas,Dmax)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),123,2,1,2,1,Nocc,Nvrt)
      Call Pccd_ohess_sen0(Work(Lt),2,1,2,1,1,4.0D0,Work(Iam),
     +                     Work(Ism),Nocc,Nvrt,Nbas,Dmax)
      Call Pccd_ohess_sen0(Work(Lt),2,1,2,1,2,-4.0D0,Work(Iam),
     +                     Work(Ism),Nocc,Nvrt,Nbas,Dmax)
      Call Pccd_ohess_unpk(Work(Lb),Work(Lt),125,2,1,2,1,Nocc,Nvrt)
      Call Pccd_ohess_sen0(Work(Lt),2,1,2,1,1,4.0D0,Work(Iam),
     +                     Work(Ism),Nocc,Nvrt,Nbas,Dmax)
      If (Dmax .Gt. 1.0D-10) Write(6,"(2a,1p,E10.3)")
     +   " @-Pccd_ohess: WARNING, 2-RDM not seniority zero, max",
     +   " dropped element",Dmax

C h = f - (2J - K), f = C^T FOCKA C at the current orbitals.
      Call Getrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,Work(Ic))
      Call Getrec(20,"JOBARC","FOCKA   ",Nsq*Iintfp,Work(Ifao))
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(Ifao),Nbas,Work(Ic),
     +           Nbas,Zero,Work(Iscr),Nbas)
      Call Dgemm("T","N",Nbas,Nbas,Nbas,One,Work(Ic),Nbas,Work(Iscr),
     +           Nbas,Zero,Work(Ihh),Nbas)
      Call Pccd_ohess_asm(Work(Igg),Work(Ihh),Work(Iff),Work(Iam),
     +                    Work(Ism),Work(Ion),Work(Iwy),Work(Ihs),
     +                    Work(Iscr),Work(Iirp),Nbas,Nocc,Npar)
      Call Putrec(20,"JOBARC","OHESFRZ ",Npar*Npar*Iintfp,Work(Ihs))
C h and the integral classes of the exact-Hessian response (Pccd_rsp).
      Call Pccd_rsp_save(Work(Igg),Work(Ihh),Work(Iwy),Nbas)
      Call Putrec(20,"JOBARC","OHESGRAD",Nsq*Iintfp,Work(Iscr))
      Iok = 1
      Return
      End
