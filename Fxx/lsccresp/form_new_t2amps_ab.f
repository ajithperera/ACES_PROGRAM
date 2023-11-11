










      Subroutine Form_new_t2amps_ab(T2amps,Work,Iwork,Maxcor,Imaxcor,
     +                              List2,T2ln,T2ln_aa,T2ln_bb,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz
      Integer T2ln,T2off,T2ln_aa,T2ln_bb
      
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end


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
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Dimension Work(Maxcor)
      Dimension T2amps(T2ln)
      Dimension Iwork(Imaxcor)

      Data Ione /1/

      Irrepx = Ione
      T2off  = Ione + T2ln_aa + T2ln_bb
      
      Do Irrep = 1, Nirrep

         Irrepr = Irrep
         Irrepl = Dirprd(Irrepr,Irrepx) 
            
          Numdis = Irpdpd(Irrepr,14)
          Dissiz = Irpdpd(Irrepl,13)

          I000 = Ione
          I010 = I000 + Numdis*Dissiz
          I020 = I010 + Numdis*Dissiz
          Iend = I020 + Numdis*Dissiz
          If (Iend .Gt. Maxcor) Call Insmem("form_new_amps_ab",
     +                                       Iend,Maxcor)

C T2(Ab,Ij)

          Call Getlst(Work(I000),1,Numdis,1,Irrepr,List2)
          Call Dcopy(Numdis*Dissiz,T2amps(T2off),1,Work(I010),1)

          Call Form_vvoo_ab(Work(I000),Work(I010),Work(I020),
     +                      Iwork,Imaxcor,Numdis,Dissiz,Irrepr,
     +                      Irrepl)

          Call Putlst(Work(I000),1,Numdis,1,Irrepr,List2)     

          T2off = T2off + Dissiz*Numdis

      Enddo 

      Return
      End
            
