






































































































































































































      Subroutine Get_mbpt2_t2s(T2amp,Work,Maxcor,T2ln,T2ln_aa,T2ln_bb,
     +                         T2ln_ab,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab
      Dimension Work(Maxcor),T2amp(T2ln)



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
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Do Ispin = 1, (Iuhf+1)
         Listw = 13 + Ispin
         Listd = 47 + Ispin
         Ibegin = (Ispin-1)*T2ln_aa
         Ioff   =  1
         Do irrep = 1, Nirrep
            Nrow = Irpdpd(Irrep,Ispin)
            Ncol = Irpdpd(Irrep,2+Ispin)
            I000 = 1 
            I010 = I000 + Nrow*Ncol
            Iend = I010 + Nrow*Ncol

            If (Iend .Ge. Maxcor) Call Insmem("get_mbpt_nos",
     +                                         Iend,Maxcor)
            Call Getlst(Work(I000),1,Ncol,2,Irrep,Listw)
            Call Getlst(Work(I010),1,Ncol,2,Irrep,ListD)
            Call vecprd(Work(I000),Work(I010),T2amp(Ibegin+Ioff),
     +                       Nrow*Ncol)
            Ioff = Ioff + Nrow*Ncol
          Enddo
      Enddo 

      Ibegin = T2ln_aa + T2ln_bb
      Ioff   = 1
      Listw  = 16
      ListD  = 50 

      Do irrep = 1, Nirrep
        Nrow = Irpdpd(Irrep,13)
        Ncol = Irpdpd(Irrep,14)
        I000 = 1
        I010 = I000 + Nrow*Ncol
        Iend = I010 + Nrow*Ncol

        If (Iend .Ge. Maxcor) Call Insmem("get_mbpt_nos",
     +                                     Iend,Maxcor)

        Call Getlst(Work(I000),1,Ncol,2,Irrep,Listw)
        Call Getlst(Work(I010),1,Ncol,2,Irrep,ListD)
        Call vecprd(Work(I000),Work(I010),T2amp(Ibegin+Ioff),
     +              Nrow*Ncol)
        Ioff = Ioff + Nrow*Ncol
      Enddo

      Call checksum("T2AA :",T2amp(1),t2ln_aa)
      If (Iuhf .Ne.0) Call checksum("T2BB :",T2amp(1+T2ln_aa),t2ln_bb)
      Call checksum("T2AB :",T2amp(1+T2ln_aa+T2ln_bb),t2ln_ab)
   
      Return
      End 
      
      
