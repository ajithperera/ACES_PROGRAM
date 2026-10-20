










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
C noise).
C The integral lists are read in column batches; every integral goes
C straight into H (the 2Y terms) and into the classes (pq|yy), (py|qy)
C (Nbas**3 each), which give h, F and the Kronecker-delta terms. Memory:
C Npar**2 + 3 Nbas**3 + one list column (no Nbas**4 array).
C
      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry
      Dimension Work(Maxcor)
      Common/Symm/Symmetry


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
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
      Call Pccd_ohess_occ(Work(Ion),Work(Iscr),Work(Iam),Work(Ism),
     +                    Work(Iirp),Nbas,Nocc)
      Dmax = Zero
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
      Call Pccd_ohess_kp(Work(Iirp),Work(Ikp),Nbas)
      Call Dzero(Work(Ihs),Npar*Npar)
      Call Dzero(Work(Ixj),2*N3)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,13,1,1,1,1,Nocc,3,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,233,2,2,2,2,Nocc,3,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,16,2,2,1,1,Nocc,3,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,10,1,1,1,2,Nocc,3,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,30,2,2,2,1,Nocc,3,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)
      Call Pccd_ohess_lst(Work(Lb),Lbuf,25,2,1,2,1,Nocc,3,One,
     +                    Work(Iam),Work(Ism),Dmax,Work(Ihs),Work(Ixj),
     +                    Work(Ixk),Work(Ikp),Nbas,Npar)

C h = f - (2J - K), f = C^T FOCKA C at the current orbitals.
      Call Getrec(20,"JOBARC","SCFEVCA0",Nsq*Iintfp,Work(Ic))
      Call Getrec(20,"JOBARC","FOCKA   ",Nsq*Iintfp,Work(Ifao))
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(Ifao),Nbas,Work(Ic),
     +           Nbas,Zero,Work(Iscr),Nbas)
      Call Dgemm("T","N",Nbas,Nbas,Nbas,One,Work(Ic),Nbas,Work(Iscr),
     +           Nbas,Zero,Work(Ihh),Nbas)
      Call Pccd_ohess_fin(Work(Ihh),Work(Iff),Work(Iam),Work(Ism),
     +                    Work(Ion),Work(Ixj),Work(Ixk),Work(Iwy),
     +                    Work(Ihs),Work(Iscr),Work(Iirp),Nbas,Nocc,
     +                    Npar)
      Call Putrec(20,"JOBARC","OHESFRZ ",Npar*Npar*Iintfp,Work(Ihs))
C h and the integral classes of the exact-Hessian response (Pccd_rsp).
      Call Putrec(20,"JOBARC","RSP_JX  ",N3*Iintfp,Work(Ixj))
      Call Putrec(20,"JOBARC","RSP_KX  ",N3*Iintfp,Work(Ixk))
      Call Putrec(20,"JOBARC","RSP_HMO ",Nsq*Iintfp,Work(Ihh))
      Call Putrec(20,"JOBARC","OHESGRAD",Nsq*Iintfp,Work(Iscr))
      Iok = 1
      Return
      End
