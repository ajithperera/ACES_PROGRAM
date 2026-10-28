      Subroutine Pccd_rsp_lay(Lwr,Nbas,Nocc,Npar)
C Workspace layout of the pair-form pCCD and its response (Pccd_rsp_init,
C Pccd_pair_drv) in Common /Pccdrsp/; Lwr: length of the workspace.
      Implicit Double Precision(A-H,O-Z)
#include "sympop.com"
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Nb  = Nbas
      No  = Nocc
      Nv  = Nbas - Nocc
      Nov = No*Nv
      Npr = Npar
      Nsq = Nb*Nb
      Np  = Nb + 2*Nsq
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
      Return
      End
