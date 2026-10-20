#include "flags.h"
      Subroutine Pccd_pair_drv(Work,Maxcor,Iuhf,Iok,Itcv,Ibad)
C
C OO-pCCD (RHF, OO_METHOD=RFO) amplitudes, lambda, energy and orbital
C gradient in pair form (Pccd_pair_slv), in place of Pccd_rdriver,
C Pccd_ldriver, Pccd_dens and Pccd_form_htau_pq: pCCD depends on the MO
C integrals only through e(p) = 2h(p,p) + (pp|pp), W(p,q) = 4(pp|qq)
C - 2(pq|qp) and K(p,q) = (pq|pq), so one pass over the integral lists
C (Pccd_ohess Imode 1) gives everything; each Newton step costs
C O(Nov**3), the CCD-shaped iterations O(o**2 v**4). On exit:
C   list 46: T2(aa,ii) = t(i,a) (the rest zero), for Pccd_rfo_tguess,
C            Pccd_rfo_ampchk and Pccd_rsp_init;
C   TOTENERG = SCFENEG + sum K t, PCCD_TCV = 1;
C   OBRTGRDA: the orbital gradient (only its antisymmetric part is
C            used, Gv = -2(G(p,q) - G(q,p)));
C   PAIR_OCC, PAIR_DA, PAIR_DS: n, A, S for Pccd_ohess Imode 2.
C Iok = 0: not done (memory, Pccd_ohess); the caller uses the CCD-shaped
C solvers. Itcv = 1: converged; Ibad = 1: max|lambda| above the cap of
C Pccd_rfo_ampchk (point rejected).
C
      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry
      Dimension Work(Maxcor)
      Common/Symm/Symmetry
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "flags.com"
      Parameter (Ampcap=5.0D0, Mxit=50)
      Data Zero,One /0.0D0,1.0D0/

      Iok  = 0
      Itcv = 0
      Ibad = 0
      Nocc = 0
      Nvrt = 0
      Do Irr = 1, Nirrep
         Nocc = Nocc + Pop(Irr,1)
         Nvrt = Nvrt + Vrt(Irr,1)
      Enddo
      Nbas = Nocc + Nvrt
      If (Nbas .Gt. 1000) Return
      Nsq  = Nbas*Nbas

      Call Pccd_setmet()
      Call Pccd_setlst(Work,Maxcor,Iuhf,"R")
      Call Pccd_init(Work,Maxcor,Iuhf,"R")
      Call Pccd_rfo_tguess(Work,Maxcor,"LOAD ")
      Call Pccd_ohess(Work,Maxcor,1,Jok)
      If (Jok .Ne. 1) Return

      Irpo = 1
      Iirp = Irpo + 8
      Iwr  = Iirp + Nbas
      Lwr  = Maxcor - Iwr + 1
      Call Pccd_rfo_irp(Work(Irpo),Work(Iirp),Nbas,Symmetry)
      Call Pccd_rfo_npar(Work(Iirp),Nbas,Npar)
      Call Pccd_rsp_size(Nbas,Nocc,Npar,Lrsp)
      If (Lrsp .Gt. Lwr) Return
      Iok = 1
      Call Pccd_rsp_lay(Lwr,Nbas,Nocc,Npar)
      Call Pccd_pair_prm(Work(Iwr),Work(Iirp))
      Ncol = Irpdpd(1,14)
      Nrow = Irpdpd(1,15)
      Call Getlst(Work(Iwr+Lsc-1),1,Ncol,1,1,46)
      Call Pccd_rsp_gett(Work(Iwr+Lsc-1),Nrow,Work(Iwr+Lt-1))

      Rtol = 10.0D0**(-Iflags(h_IFLAGS_cc_conv))
      Call Pccd_pair_slv(Work(Iwr),Work(Iirp),Npar,Rtol,Mxit,Ec,Eref,
     +                   Alm,Nit,Rn,Iconv)
      Write(6,"(a,i3,a,1p,E10.3)") " @-PAIR pCCD equations:",Nit,
     +     " Newton steps, max|r| =",Rn
      If (Iconv .Ne. 1) Then
         Write(6,"(2a)") " @-PAIR the pair-form pCCD equations did",
     +                   " not converge."
         Return
      Endif
      Itcv = 1
      Call Putrec(20,"JOBARC","PCCD_TCV",1,1)

C T2 to list 46 (all irreps zeroed first).
      Call Zersym(Work(Iwr+Ldv-1),46)
      Call Pccd_pair_putt(Work(Iwr+Lt-1),Work(Iwr+Lsc-1),Nrow)

      Call Getrec(20,"JOBARC","SCFENEG ",Iintfp,Escf)
      Call Getrec(20,"JOBARC","NUCREP  ",Iintfp,Enuc)
      Etot = Escf + Ec
      Call Putrec(20,"JOBARC","TOTENERG",Iintfp,Etot)
      Write(6,"(a,F20.12,a,F20.12)") " @-PAIR pCCD energy",Etot,
     +     "  correlation",Ec
      Dref = Eref + Enuc - Escf
      If (Dabs(Dref) .Gt. 1.0D-8) Write(6,"(2a,1p,E10.3)")
     +   " @-PAIR WARNING: pair-form reference energy - SCFENEG =",Dref

      Write(6,"(a,F12.6)") " @-RFO L2 max|amp| (pair lambda) =",Alm
      If (Alm .Gt. Ampcap) Then
         Ibad = 1
         Write(6,"(a,F6.2,a)") " @-RFO L2 amplitudes exceed",Ampcap,
     +         " : unphysical solution, point rejected."
         Return
      Endif

      Call Putrec(20,"JOBARC","PAIR_OCC",Nbas*Iintfp,Work(Iwr+Lg2-1))
      Call Putrec(20,"JOBARC","PAIR_DA ",Nsq*Iintfp,
     +            Work(Iwr+Lg2+Nbas-1))
      Call Putrec(20,"JOBARC","PAIR_DS ",Nsq*Iintfp,
     +            Work(Iwr+Lg2+Nbas+Nsq-1))
      Call Pccd_pair_grd(Work(Iwr+Ixs-1),Work(Iirp),Work(Iwr+Lz-1),
     +                   Nbas)
      Call Putrec(20,"JOBARC","OBRTGRDA",Nsq*Iintfp,Work(Iwr+Lz-1))
      Return
      End
