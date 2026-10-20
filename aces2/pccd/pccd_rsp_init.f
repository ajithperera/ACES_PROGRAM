










      Subroutine Pccd_rsp_init(Wr,Lwr,Nbas,Nocc,Npar,Irp,Gv,Iok,Iprt)
C Read h, (pq|yy), (py|qy) and T2 (list 46), form the pair parameters,
C polish t by Newton steps on the pair equations, factor A = dr/dt,
C solve A^T lam = -K. Checks: |r| of the ACES amplitudes and the orbital
C gradient of L against Gv (Pccd_rfo_do convention). Iok = 1: ready.
C Iprt = 1: print the checks (otherwise only when they fail).
      Implicit Double Precision(A-H,O-Z)
      Dimension Wr(Lwr),Irp(Nbas),Gv(Npar),Nfzo(8),Nfzv(8)


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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Data Zero,One /0.0D0,1.0D0/

      Iok = 0
      If (Nbas .Gt. 1000) Return
      Nb  = Nbas
      No  = Nocc
      Nv  = Nbas - Nocc
      Nov = No*Nv
      Npr = Npar
      Nsq = Nb*Nb
      Np  = Nb + 2*Nsq
      Call Pccd_rsp_size(Nbas,Nocc,Npar,Lrsp)
      If (Lrsp .Gt. Lwr) Return
      Call Getrec(0,"JOBARC","RSP_JX  ",Length,Junk)
      If (Length .Ne. Nsq*Nb*Iintfp) Return
      Lh  = 1
      Ljx = Lh  + Nsq
      Lkx = Ljx + Nsq*Nb
      Le  = Lkx + Nsq*Nb
      Lw  = Le  + Nb
      Lk  = Lw  + Nsq
      Lde = Lk  + Nsq
      Ldw = Lde + Nb
      Ldk = Ldw + Nsq
      Lg1 = Ldk + Nsq
      Lg2 = Lg1 + Np
      Lg3 = Lg2 + Np
      Lz  = Lg3 + Np
      Lbm = Lz  + Nsq
      Lt  = Lbm + Nsq
      Llm = Lt  + Nov
      Lmk = Llm + Nov
      Lt1 = Lmk + Nov
      Ll1 = Lt1 + Nov
      Lvw = Ll1 + Nov
      Lr1 = Lvw + Nov
      Lr2 = Lr1 + Nov
      La  = Lr2 + Nov
      Lpv = La  + Nov*Nov
      Lsc = Lpv + Nov
      Nblk = Irpdpd(1,14)*Irpdpd(1,15)
C Ixs: Npar words for the relaxation vector of Pccd_rsp_hv.
      Ixs = Lsc + Max(Nblk,Nov*Nov)
      Ldv = Ixs + Npar
      Ldvl = Lwr - Ldv + 1
      Call Getrec(20,"JOBARC","RSP_HMO ",Nsq*Iintfp,Wr(Lh))
      Call Getrec(20,"JOBARC","RSP_JX  ",Nsq*Nb*Iintfp,Wr(Ljx))
      Call Getrec(20,"JOBARC","RSP_KX  ",Nsq*Nb*Iintfp,Wr(Lkx))
      Call Pccd_rsp_par(Wr(Lh),Wr(Ljx),Wr(Lkx),Wr(Le),Wr(Lw),Wr(Lk),
     +                  Nb,1)

C Frozen-pair mask: the first Nfzo occupied and last Nfzv virtual MOs
C of each irrep (occupied/virtual by irrep, as Irp).
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

C T2 pair amplitudes from list 46 (irrep-1 block: rows A-b, columns I-j).
      Ncol = Irpdpd(1,14)
      Nrow = Irpdpd(1,15)
      Call Getlst(Wr(Lsc),1,Ncol,1,1,46)
      Call Pccd_rsp_gett(Wr(Lsc),Nrow,Wr(Lt))
      Do K = 1, Nov
         Wr(Lt+K-1) = Wr(Lt+K-1)*Wr(Lmk+K-1)
      Enddo

C Residual of the ACES amplitudes, then Newton polish.
      Call Pccd_rsp_res(Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),Wr(Lr1))
      R0 = Zero
      Do K = 1, Nov
         R0 = Max(R0,Dabs(Wr(Lr1+K-1)))
      Enddo
      Do It = 1, 5
         Call Pccd_rsp_jac(Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),Wr(La),
     +                     Wr(Lt1),Wr(Lr1),Wr(Lr2))
         Call Pccd_rsp_res(Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),Wr(Lr1))
         Rn = Zero
         Do K = 1, Nov
            Rn = Max(Rn,Dabs(Wr(Lr1+K-1)))
         Enddo
         If (Rn .Lt. 1.0D-12) Go to 10
         Call Dgetrf(Nov,Nov,Wr(La),Nov,Wr(Lpv),Info)
         If (Info .Ne. 0) Return
         Call Dgetrs("N",Nov,1,Wr(La),Nov,Wr(Lpv),Wr(Lr1),Nov,Info)
         Call Daxpy(Nov,-One,Wr(Lr1),1,Wr(Lt),1)
      Enddo
 10   Continue
      If (Iprt .Eq. 1 .Or. Rn .Gt. 1.0D-8)
     +   Write(6,"(2a,1p,E10.3,a,E10.3)") " @-RFO exact Hessian: pair",
     +   " equations |r(T2)| =",R0,", after Newton",Rn
      If (Rn .Gt. 1.0D-8) Return
      Call Pccd_rsp_jac(Wr(Lt),Wr(Le),Wr(Lw),Wr(Lk),Wr(Lmk),Wr(La),
     +                  Wr(Lt1),Wr(Lr1),Wr(Lr2))
      Call Dgetrf(Nov,Nov,Wr(La),Nov,Wr(Lpv),Info)
      If (Info .Ne. 0) Return
C A^T lam = -K(ov) (frozen pairs: lam from r = t, no P dependence)
      Call Pccd_rsp_kov(Wr(Lk),Wr(Llm),Nb,No)
      Call Dscal(Nov,-One,Wr(Llm),1)
      Call Dgetrs("T",Nov,1,Wr(La),Nov,Wr(Lpv),Wr(Llm),Nov,Info)

C Consistency: orbital gradient of L(t, lam) vs the ACES gradient.
      Call Pccd_rsp_gen(Wr(Lt),Wr(Lmk),Wr(Lg1),Nb,No)
      Call Pccd_rsp_gphi(Wr(Llm),Wr(Lt),Wr(Lmk),Wr(Lg2),Nb,No)
      Call Daxpy(Np,One,Wr(Lg2),1,Wr(Lg1),1)
      Call Pccd_rsp_zgr(Wr(Lg1),Wr(Lh),Wr(Ljx),Wr(Lkx),Wr(Lz),Irp,Nb,
     +                  Npar,Wr(Ixs))
      Dg = Zero
      Gm = Zero
      Do K = 1, Npar
         Dg = Max(Dg,Dabs(Wr(Ixs+K-1)-Gv(K)))
         Gm = Max(Gm,Dabs(Gv(K)))
      Enddo
      If (Iprt .Eq. 1 .Or. Dg .Gt. 1.0D-6)
     +   Write(6,"(2a,1p,E10.3,a,E10.3,a)") " @-RFO exact Hessian:",
     +   " pair-form gradient vs RFO gradient",Dg," (max|g|",Gm,")"
      If (Dg .Gt. 1.0D-6) Return
      Iok = 1
      Return
      End
