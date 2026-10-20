














































































































































































































      Subroutine Pccd_do_orbrots(Work,Maxcor,Iuhf,IGrad_calc,IHess_calc,
     +                           Scale,OOmicroItr,OOmacroItr,OOtotalItr)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Grdcnv
      Logical pCCD,CCD,LCCD
      Logical OneP_ONLY,Hdump,Onep_save,Ohess

      Common /CALC/pCCD,CCD,LCCD
      Common /ORBR_HESS/ONEP_ONLY

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end


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

      Data Ione,Itwo,Izero,One,Half /1,2,0,1.0D0,0.50D0/

      Grdcnv = .False.
      Ncycle = 0
      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)
C RFO_FAIL: unconverged amplitudes, no densities; Pccd_rfo_step only
C rejects the step.
      Call Getrec(0,"JOBARC","RFO_FAIL",Length,Junk)
      Ifail = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_FAIL",1,Ifail)

      If (Ifail .Eq. 1) Then
         Continue
      Else If (Iuhf .EQ. 0 .AND. Idrop .EQ. 0) Then
C Analytic frozen-amplitude orbital Hessian (pCCD): reads Gamma lists
C 123/125 before Pccd_form_htau_pq reorders them. Debug switch
C pccd_ohess.in (gradient check and pccd_ohess.dat).
         Inquire(File="pccd_ohess.in",Exist=Ohess)
C OO_METHOD=RFO unless pccd_rfo_h0.in says CONST (not at the displaced
C points of the stability check).
         Ih0 = 0
         If (Iflags2(178) .Eq. 6) Then
            Call Pccd_rfo_h0mode(Ih0)
            Call Getrec(0,"JOBARC","RFO_PHAS",Length,Junk)
            Iphase = 0
            If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_PHAS",1,
     +                                     Iphase)
            If (Iphase .Eq. 1) Ih0 = 0
         Endif
C PCCDPAIR = 1 (set by Pccd_main for this run): the amplitudes, lambda,
C energy and gradient come from the pair-form solver (Pccd_pair_drv);
C there are no ACES densities, so the frozen Hessian is formed from the
C pair densities and the gradient lists are skipped.
         Call Getrec(0,"JOBARC","PCCDPAIR",Length,Junk)
         Ipair = 0
         If (Length .Gt. 0) Call Getrec(20,"JOBARC","PCCDPAIR",1,Ipair)
         Iok = 0
         If (Ipair .Eq. 1) Then
            If (Ih0 .Ge. 1) Call Pccd_ohess(Work,Maxcor,2,Iok)
            If (Iok .Eq. 1) Call Pccd_pair_chk(Work,Maxcor)
            Call Putrec(20,"JOBARC","OHESOK  ",1,Iok)
            Go to 10
         Endif
         If (pCCD .And. (Ohess .Or. Ih0 .Ge. 1))
     +      Call Pccd_ohess(Work,Maxcor,0,Iok)
         Call Putrec(20,"JOBARC","OHESOK  ",1,Iok)
C Debug (file pccd_hess.in): dump the two-particle density and all
C other lists 1-250 as they are (GA, AL), then form the gradient and
C the analytic orbital Hessian lists once, in the production order
C (Pccd_form_htau_pq reorders Gamma lists 123/125 in place, so it must
C not run twice), and dump the Hessian lists (FU).
         Inquire(File="pccd_hess.in",Exist=Hdump)
         If (Hdump) Then
            Call Pccd_hess_dump(Work,Maxcor,"GA  ")
            Call Pccd_hess_dump(Work,Maxcor,"AL  ")
            Call Pccd_hess_dump(Work,Maxcor,"ZERO")
            Onep_save = ONEP_only
            ONEP_only = .False.
            Mxc = Maxcor
            Call Pccd_form_htau_pqrs(Work,Mxc,Iuhf)
            Mxc = Maxcor
            Call Pccd_form_htau_pq(Work,Mxc,Iuhf,Igrad_calc,1,Scale)
            Call Pccd_hess_dump(Work,Maxcor,"FU  ")
            ONEP_only = Onep_save
         Else
         If(.Not. ONEP_only) Call Pccd_form_htau_pqrs(Work,Maxcor,Iuhf)
C OO_METHOD=RFO does not use the (unfinished, undebugged) analytic
C Hessian lists that Pccd_form_htau_1D_pqrs builds after the gradient.
         Ihess = IHess_calc
         If (Iflags2(178) .Eq. 6) Ihess = 0
         Call Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          Ihess,Scale)
         Endif
         If (pCCD .And. Ohess) Call Pccd_ohess_chk(Work,Maxcor)
 10      Continue
      Else

C Uhf only or UHF and RHF for forzen core orbitals.

         Call Pccd_form_uhtau(Work,Maxcor,Iuhf)
      Endif 

      Call GetSDInfo(nocca,noccb,nvrta,nvrtb,nbas,nvecDim,Naobfns)

      If (Iflags2(178) .EQ. 1) then ! BB Steep. Desc.
         call Drive_BB_SD(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                    nocca,nvrta,nvecDim)
      else if (Iflags2(178) .EQ. 2) then !BFGS
        print*,'Have not include L-BFGS yet; to do'
      else if (Iflags2(178) .EQ. 3) then !FullNR
        print*,'Have not included FullNR yet; to do'
      else if (Iflags2(178) .EQ. 4) then !AMSgrad
          call Drive_AMSgrad(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                        nocca,nvrta,nvecDim)
      else if (Iflags2(178) .EQ. 6) then !RFO
          Tol  = 10.0D0**(-Iflags(76))
          Etol = 10.0D0**(1-Iflags(4))
          Call Pccd_rfo_step(Work,Maxcor,Nbas,Nocca,Iuhf,Grdcnv,Tol,
     +                       Etol,Grd_max,Grd_rms,Ncycle)
      else if (Iflags2(178) .EQ. 5) then !AMSgrad
          Call Getrec(-20,'JOBARC','ORBOPITR',Ione,Ncycle)
          Call Pccd_rotgrd(Work,Maxcor,Nbas,Nocca,Noccb,Nvrta,Nvrtb,
     +                     Iuhf,Grdcnv,Tol,Grd_max,Grd_rms,Ncycle)
      endif
C 
      If (Grdcnv .And. .Not. Ncycle .Eq. 0) Call Pccd_finalize(Work,
     +                                      Maxcor,Iuhf,Tol,Grd_max,
     +                                      Grd_rms)
      If (IFLAGS(11).GT.0) 
     +IFLAGS(11) = Ione
      IFLAGS(16)= Izero
      IFLAGS(38)     = Ione
C Do not turn on non-HF key-word sinply because this indirectly add singles 
      IFLAGS(34)  =Izero
      IFLAGS(45)     =Izero

      Call Putrec(20,'JOBARC','IFLAGS  ',100,IFLAGS)

      Return
      End 
 
