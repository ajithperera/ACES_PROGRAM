










      Subroutine Pccd_rsp_size(Nbas,Nocc,Npar,Lrsp)
C Words needed by Pccd_rsp_init/hv and the Davidson of Pccd_rsp_stab.
      Implicit Double Precision(A-H,O-Z)
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
      Parameter (Mx=40)
      Nvrt = Nbas - Nocc
      Nov  = Nocc*Nvrt
      Nsq  = Nbas*Nbas
      Np   = Nbas + 2*Nsq
      Nblk = Irpdpd(1,14)*Irpdpd(1,15)
      Ldav = 2*Npar*Mx + 2*Npar + 2*Mx*Mx + Mx + 64*Mx
      Lrsp = Nsq + 2*Nsq*Nbas + 5*Np + 2*Nsq + 10*Nov + Nov*Nov
     +     + Nov + Max(Nblk,Nov*Nov) + Ldav + 8*Npar + 8
      Return
      End
