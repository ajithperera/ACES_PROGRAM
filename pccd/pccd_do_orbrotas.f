













































































































































































































      Subroutine Pccd_do_orbrots(Work,Maxcor,Iuhf,IGrad_calc,IHess_calc,
     +                           Scale,OOmicroItr,OOmacroItr,OOtotalItr)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Grdcnv
      Logical pCCD,CCD,LCCD
      Logical OneP_ONLY

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

      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)
       
      If (Iuhf .EQ. 0 .AND. Idrop .EQ. 0) Then
         If(.Not. ONEP_only) Call Pccd_form_htau_pqrs(Work,Maxcor,Iuhf)
         Call Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          IHess_calc,Scale)
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
 
