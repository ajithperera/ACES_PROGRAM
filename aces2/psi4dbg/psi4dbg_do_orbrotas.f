













































































































































































































      Subroutine Psi4dbg_do_orbrots(Work,Maxcor,Iuhf,IGrad_calc,
     +                              IHess_calc,Scale,OOmicroItr,
     +                              OOmacroItr,OOtotalItr)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)

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

      Data Ione,Itwo,One,Half /1,2,1.0D0,0.50D0/

      Call Psi4dbg_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          IHess_calc,Scale)

      Call GetSDInfo(nocc,nvirt,nbas,nvecDim,Naobfns)


      If (Iflags2(178) .EQ. 0) then ! BB Steep. Desc.
         call Drive_BB_SD(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                    nocc,nvirt,nvecDim)
      else if (Iflags2(178) .EQ. 1) then !BFGS
        print*,'Have not include L-BFGS yet; to do'
      else if (Iflags2(178) .EQ. 2) then !FullNR
        print*,'Have not included FullNR yet; to do'
      else if (Iflags2(178) .EQ. 3) then !AMSgrad
          call Drive_AMSgrad(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                        nocc,nvirt,nvecDim)
      else if (Iflags2(178) .EQ. 4) then !AMSgrad
          Call Psi4dbg_rotgrd(Work,Maxcor,Nbas,Nocc,Nvirt)
      endif

      IFLAGS(16)=0 
C      IFLAGS(h_IFLAGS_non_hf)=1 
      CALL PUTREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
      Return
      End 
 
