










      Subroutine Form_t2aaaa(T2ab,Work,Maxcor,Lenab,Lenaa)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor),T2ab(Lenab)



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
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

      Data Ione,Inul /1,0/
      
      Ioff   = Ione 
      Joff   = Inul
      Irrepx = Ione

      call getall(T2ab,Lenab,Irrepx,63)

      Call checksum("T2-ABAB :",T2ab,Lenab)
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepx,Irrepr)
 
         Ncol_isqj = Irpdpd(Irrepr,14)
         Nrow_asqb = Irpdpd(Irrepl,15)
         Ncol_itrJ = Irpdpd(Irrepr,3)
         Nrow_atrb = Irpdpd(Irrepl,1)

         I000 = Ione
         Iend = I000 + Lenaa
         If (Iend .Gt. Maxcor) Call Insmem("form_t2aaaa",Iend,
     +                                      Maxcor)
         Call Assym2(Irrepl,Pop(1,1),Nrow_asqb,T2ab(Ioff))

         Call Sqsym(Irrepl,Vrt(1,1),Nrow_atrb,Nrow_asqb,
     +              Ncol_itrJ,Work(I000+joff),T2ab(Ioff))

         Ioff = Ncol_isqj*Nrow_asqb + Ioff
         Joff = Ncol_itrJ*Nrow_atrb + Joff
      Enddo 

      Call Dcopy(Lenaa,Work(I000),1,T2ab,1)

      Write(6,"(a)") " Forming T2AAAA from T2ABAB"
      Call checksum("T2-AAAA :",T2ab,Lenaa)

      Return
      End

