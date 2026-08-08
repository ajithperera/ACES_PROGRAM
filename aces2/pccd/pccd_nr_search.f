













































































































































































































      Subroutine Pccd_nr_search(Grd,Grd_oo,Grd_vv,Grd_vo,Grd_ov,
     +                          Grd_stata,Grd_statb,B_like,Lenoo,
     +                          Lenvv,Lenvo,Nbas,Nocc,Nvrt,Work,
     +                          Maxcor,Iuhf,Ispin,Conv_tol,Icycle)

      Implicit Double Precision(A-H,O-Z)
      Logical B_like, Sym_packed 
      Logical Bfgs
      Logical Symmetry, Apprx_CC_hess, Apprx_HF_hess 
     
      Dimension Grd_stata(6)
      Dimension Grd_statb(6)
      Dimension SGrd_stata(6)
      Dimension SGrd_statb(6)
      Dimension Grd(Nbas,Nbas)
      Dimension Grd_oo(Lenoo)
      Dimension Grd_vv(Lenvv)
      Dimension Grd_vo(Lenvo)
      Dimension Grd_ov(Lenvo)
      Dimension Work(Maxcor)



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
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end

      Common /Symm/Symmetry

      Data Ione,Izero,One,Half,Dnull,Two,Onem/1,0,1.0D0,0.50D0,
     +                                        0.0D0,2.0D0,-1.0D0/
 
      B_like = .False.
      Irrepx = Ione
      Itrln  = Nbas*(Nbas-1)/2
      If (Symmetry) Then
         Nsize_vo = Idsymsz(Irrepx,8+Ispin,8+Ispin)
         Nsize_ov = Nsize_vo
         Nsize_oo = Idsymsz(Irrepx,20+Ispin,20+Ispin)
         Nsize_vv = Idsymsz(Irrepx,18+Ispin,18+Ispin)
      Else
         Nsize_vo = (Nocc*Nvrt)**2
         Nsize_ov = Nsize_vo
         Nsize_oo = (Nocc*Nocc)**2
         Nsize_vv = (Nvrt*Nvrt)**2
      Endif 

      Nsize    = Nbas*Nbas
      Lenoo    = Nocc*Nocc
      Lenvv    = Nvrt*Nvrt
      Lenvo    = Nvrt*Nocc
      Lenov    = Lenvo
  
      IHess_update = Iflags2((7))
      If (IHess_update .Eq. 2) Bfgs = .True. 

      If (Bfgs) Then 
     
      I030 = Ione
      I040 = I030 + Nsize_vo
      I050 = I040 + Nsize_vo
      I060 = I050 + Nsize_oo
      I070 = I060 + Nsize_vv
      I080 = I070 + Max(Nsize_oo,Nsize_vo,Nsize_vv)
      I090 = I080 + Max(Nsize_oo,Nsize_vo,Nsize_vv)
      Iend = I090 + Nsize
      Memleft = Maxcor - Iend 
      If (Iend .Gt. Maxcor) Call Insmem("pccd_nr_search",Iend,Maxcor)
 
      If (Icycle .Eq. Izero) Then
         Call Pccd_hess_ov(Work(I040),Work(Iend),Nsize_vo,Memleft,
     +                     Ispin,Iuhf)
         Call Pccd_hess_oo(Work(I050),Work(Iend),Nsize_oo,Memleft,
     +                     Ispin,Iuhf)
         Call Pccd_hess_vv(Work(I060),Work(Iend),Nsize_vv,Memleft,
     +                     Ispin,Iuhf)
         Call set_diags2_one(Work(I040),Work(I050),Work(I060),Lenoo,
     +                       Lenvv,Lenvo)

C Store the current Hessians and gradients in Jobarc. We may move these
C to a unformatted file if the need arise. 

         Call Putrec(20,"JOBARC","VO_ROT_H",Nsize_vo*Iintfp,Work(I040))
         Call Putrec(20,"JOBARC","OO_ROT_H",Nsize_oo*Iintfp,Work(I050))
         Call Putrec(20,"JOBARC","VV_ROT_H",Nsize_vv*Iintfp,Work(I060))

         Call Putrec(20,"JOBARC","VO_ROT_G",Lenvo*Iintfp,Grd_vo)
         Call Putrec(20,"JOBARC","OV_ROT_G",Lenvo*Iintfp,Grd_ov)
         Call Putrec(20,"JOBARC","OO_ROT_G",Lenoo*Iintfp,Grd_oo)
         Call Putrec(20,"JOBARC","VV_ROT_G",Lenvv*Iintfp,Grd_vv)
      Else
         Lhess = Max(Nsize_vo,Nsize_oo,Nsize_vv)
         Lgrad = Max(Lenvo,Lenoo,Lenvv)
         I000  = Ione
         I010  = I000 + Lhess*Iintfp
         I040  = I010 + Lgrad*Iintfp
         I050  = I040 + Nsize_vo*Iintfp
         I060  = I050 + Nsize_oo*Iintfp
         I070  = I060 + Nsize_vv*Iintfp
         I080  = I070 + Max(Nsize_oo,Nsize_vo,Nsize_vv)
         I090  = I080 + Max(Nsize_oo,Nsize_vo,Nsize_vv)
         Iend  = I090 + Nsize
         Memleft = Maxcor - Iend 
         If (Iend .Gt. Maxcor) Call Insmem("pccd_nr_search",Iend,Maxcor)

         Call Pccd_hess_update(Work(I000),Work(I010),Work(I040),
     +                         Work(I050),Work(I060),Grd_oo,Grd_vv,
     +                         Grd_vo,Grd_ov,Work(Iend),Memleft,
     +                         Lhess,Lgrad,Nocc,Nvrt,Lenoo,Lenvv,
     +                         Lenvo,Nsize_oo,Nsize_vv,Nsize_vo,
     +                         Icycle)
C
C Store the current Hessians and gradients in Jobarc. We may move these
C to a unformatted file if the need arise. 

         Call Putrec(20,"JOBARC","VO_ROT_H",Nsize_vo*Iintfp,Work(I040))
         Call Putrec(20,"JOBARC","OO_ROT_H",Nsize_oo*Iintfp,Work(I050))
         Call Putrec(20,"JOBARC","VV_ROT_H",Nsize_vv*Iintfp,Work(I060))

         Call Putrec(20,"JOBARC","VO_ROT_G",Lenvo*Iintfp,Grd_vo)
         Call Putrec(20,"JOBARC","OV_ROT_G",Lenvo*Iintfp,Grd_ov)
         Call Putrec(20,"JOBARC","OO_ROT_G",Lenoo*Iintfp,Grd_oo)
         Call Putrec(20,"JOBARC","VV_ROT_G",Lenvv*Iintfp,Grd_vv)
      Endif 

      Call output(Work(I040),1,lenvo,1,lenvo,lenvo,lenvo,1)
      Call Dcopy(Nsize_vo,Work(I040),1,Work(I070),1)
      Call Eig(Work(I070),Work(I080),1,Lenvo,1)
      Call Pccd_form_hinv(Work(I040),Work(Iend),Memleft,Nbas,
     +                    Lenvo,Nocc,Nvrt,"OV")
      Write(6,*)
      Write(6,"(a)") " The eigenvalues of the OV/VO Hessian"
      Call output(Work(I070),1,Lenvo,1,Lenvo,Lenvo,Lenvo,1)
      Write(6,*)
      Write(6,"(a)") " The OV/VO orbital rotation Hessian inverse"
      Call output(Work(I040),1,lenvo,1,lenvo,lenvo,lenvo,1)
      Call Dcopy(Nsize_oo,Work(I050),1,Work(I070),1)
      Call Eig(Work(I070),Work(I080),1,Lenoo,1)
      Call Pccd_form_hinv(Work(I050),Work(Iend),Memleft,Nbas,
     +                    Lenoo,Nocc,Nvrt,"OO")
      Write(6,*)
      Write(6,"(a)") " The eigenvalues of the OO Hessian"
      Call output(Work(I070),1,Lenoo,1,Lenoo,Lenoo,Lenoo,1)
      Write(6,*)
      Write(6,"(a)") " The OO orbital rotation Hessian inverse"
      Call output(Work(I050),1,lenoo,1,lenoo,lenoo,lenoo,1)
      Call Dcopy(Nsize_vv,Work(I060),1,Work(I070),1)
      Call Eig(Work(I070),Work(I080),1,Lenvv,1)
      Call Pccd_form_hinv(Work(I060),Work(Iend),Memleft,Nbas,
     +                    Lenvv,Nocc,Nvrt,"VV")
      Write(6,*)
      Write(6,"(a)") " The eigenvalues of the VV Hessian"
      Call output(Work(I070),1,Lenvv,1,Lenvv,Lenvv,Lenvv,1)
      Write(6,*)
      Write(6,"(a)") " The VV orbital rotation Hessian inverse"
      Call output(Work(I060),1,lenvv,1,lenvv,lenvv,lenvv,1)
      Call Pccd_nr_update(Work(I040),Work(I050),Work(I060),Grd_oo,
     +                    Grd_vv,Grd_vo,Grd_ov,Work(Iend),Maxcor,
     +                    Lenoo,Lenvv,Lenvo,Nocc,Nvrt,.True.)
      Sym_packed = .False. 
      Call Pccd_frmful(Work(I090),Grd_oo,Grd_vv,Grd_vo,Grd_ov,
     +                 Work(Iend),Memleft,Nocc,Nvrt,Nbas,"Ov_like",
     +                 Sym_packed)
      Call Dcopy(Nbas*Nbas,Work(I090),1,Grd,1)

      Call Pccd_gossip(Grd,Grd_stata,Grd_statb,Sgrad_stata,
     +                 Sgrad_statb,Nbas,Ispin)

      Else

      I000 = Ione
      I010 = I000 + Nbas
      I020 = I010 + Lenvo
      I030 = I020 + Lenvo
      I040 = I030 + Lenoo
      I050 = I040 + Lenvv
      Iend = I050 + Nsize 
      Memleft = Maxcor - Iend
      If (Iend .Gt. Maxcor) Call Insmem("pccd_nr_search",Iend,Maxcor)

      If (Ispin .Eq. 1) Call Getrec(20,"JOBARC","SCFEVALA",
     +                              Nbas*Iintfp,Work(I000))
      If (Ispin .Eq. 2) Call Getrec(20,"JOBARC","SCFEVALB",
     +                              Nbas*Iintfp,Work(I000))


      Indi   = Izero
      Inda   = Izero 
      Indv   = Izero 
      Indo   = Izero 
      Indv0  = Izero 
      Indo0  = Izero 
      Ijunk  = Izero
      Irrepx = Ione 

      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Nocci = Pop(Irrepr,Ispin)
         Nvrti = Vrt(Irrepl,Ispin)

         Do indxi = 1, Nocci
            Indi = Indi + 1
            Indv = Indv0
            Do indxa = 1, Nvrti
               Inda = Inda + 1
               Indv = Indv + 1
               Work(I010+Inda-1) = (Work(I000+Indv-1+Nocc) - 
     +                              Work(I000+Indi-1))
            Enddo
         Enddo
         Indv0 = Indv0 + Nvrti
      Enddo 

      Inda   = Izero 
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Nocci = Pop(Irrepr,Ispin)
         Noccj = Pop(Irrepl,Ispin)

         Do indxi = 1, Nocci
            Do indxj = 1, Noccj
               Inda = Inda + 1
               Work(I030+Inda-1) = (Work(I000+Indxi-1) -
     +                              Work(I000+Indxj-1))
            Enddo
         Enddo
      Enddo

      Inda   = Izero 
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Nvrta = Vrt(Irrepr,Ispin)
         Nvrtb = Vrt(Irrepl,Ispin)

         Do indxa = 1, Nvrta
            Do indxb = 1, Nvrtb
               Inda = Inda + 1
               Work(I040+Inda-1) = (Work(I000+Nocc+Indxa-1) -
     +                              Work(I000+Nocc+Indxb-1))
            Enddo
         Enddo
      Enddo

C For the time being lets leave occ-occ and vrt-vrt Hessian as zero.

      If (Symmetry) Call Pccd_symexp(Work(I030),Work(I040),Work(I010),
     +                               Work(Iend),Memleft,Nocc,Nvrt,Nbas)

      Call Transp(Work(I010),Work(I020),Nocc,Nvrt)

      Write(6,"(a)") "The OO/VV orbital rotation hessians"
      Call output(work(I030),1,nocc,1,nocc,nocc,nocc,1)
      Call output(work(I040),1,nvrt,1,nvrt,nvrt,nvrt,1)
      Write(6,*)
      Write(6,"(a)") "The VO/OV orbital rotation hessians"
      Call output(work(I010),1,nvrt,1,nocc,nvrt,nocc,1)
      Call output(work(I020),1,nocc,1,nvrt,nocc,nvrt,1)


C -----------------debug end-------------------------------------
C Note that I am using 1/(ea-ei) instead of 1/2(ea-ei) since I work with only
C the alpha (or beta block).  The main loop that rotates orbitals goes over
C only alpha block only (in the RHF context). See pccd_rotg.F.

C -----------------------This is the current working code--------------
      Sym_packed = .False.
      Call Pccd_scale_grads(Grd_oo,Grd_vv,Grd_vo,Grd_ov,Work(I010),
     +                      Work(I020),Work(I030),Work(I040),Lenoo,
     +                      Lenvv,Lenvo,Nocc,Nvrt,Icycle)
      Call Pccd_frmful(Work(I050),Grd_oo,Grd_vv,Grd_vo,Grd_ov,
     +                 Work(Iend),Memleft,Nocc,Nvrt,Nbas,"Ov_like",
     +                 Sym_packed)
      Call Dcopy(Nbas*Nbas,Work(I050),1,Grd,1)

      Call Pccd_gossip(Grd,Grd_stata,Grd_statb,Sgrad_stata,
     +                 Sgrad_statb,Nbas,Ispin)

c ------------Endif for Apprx_CC_hess or Apprx_HF_hess--------------------
      Endif 

C----------------------------------------------------------------------------
   
      Return
      End

