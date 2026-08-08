










      Subroutine Pccd_symexp(Doo,Dvv,Dvo,Work,Maxcor,Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry

      Dimension Doo(Nocc*Nocc)
      Dimension Dvv(Nvrt*Nvrt)
      Dimension Dvo(Nvrt*Nocc)
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

      Common /Symm/ Symmetry 
      
      Data Ione,Itwo,Inull/1,2,0/

      Ioffv   = Ione
      Ioffo   = Ione
      Ioffvo  = Ione
      Irrepx  = Ione

      If (Symmetry) Then
         Nocc2 = Irpdpd(Irrepx,21)
         Nvrt2 = Irpdpd(Irrepx,19)
         Nvovo = Irpdpd(Irrepx,9)
      Else
         Nvrt2   = Nvrt*Nvrt
         Nocc2   = Nocc*Nocc
         Nvovo   = Nocc*Nvrt
      Endif 

      Lenoo = Nocc*Nocc
      Lenvv = Nvrt*Nvrt
      Lenvo = Nvrt*Nocc
      Lenmx = Max(Lenoo,Lenvv,Lenvo)

      Do Ispin = Ione, Ione

         I000 = Ione 
         I010 = I000 + Lenmx 
         I020 = I010 + Lenmx
         Iend = I020 + Lenmx
         If (Iend .Gt. Maxcor) Call Insmem("pccd_symexp",Iend,Maxcor)
         Call Dzero(Work(I000),3*Lenmx)
         Call Pccd_blockcopy(Work(I000),Doo(Ioffo),"OO",Work(I010),
     +                       Work(I020),Nocc,Irpdpd(Irrepx,20+Ispin),
     +                       Nocc2,Ispin)
         Call Dcopy(Lenoo,Work(I000),1,Doo,1)

         I000 = Ione
         I010 = I000 + Lenmx
         I020 = I010 + Lenmx
         Iend = I020 + Lenmx
         If (Iend .Gt. Maxcor) Call Insmem("pccd_symexp",Iend,Maxcor)
         Call Dzero(Work(I000),3*Lenmx)
         Call Pccd_blockcopy(Work(I000),Dvv(Ioffv),"VV",Work(I010),
     +                       Work(I020),Nvrt,Irpdpd(Irrepx,18+Ispin),
     +                       Nvrt2,Ispin)
         Call Dcopy(Lenvv,Work(I000),1,Dvv,1)

         I000 = Ione
         I010 = I000 + Lenmx
         I020 = I010 + Lenmx
         Iend = I020 + Lenmx
         If (Iend .Gt. Maxcor) Call Insmem("pccd_symexp",Iend,Maxcor)
         Call Dzero(Work(I000),3*Lenmx)
         Call Pccd_blockcopy(Work(I000),Dvo(Ioffvo),"VO",Work(I010),
     +                       Work(I020),Nvrt,Irpdpd(Irrepx,8+Ispin),
     +                       Nvovo,Ispin)
         Call Dcopy(Lenvo,Work(I000),1,Dvo,1)

         Ioffo  = Ioffo  + Irpdpd(Irrepx,20+Ispin)
         Ioffv  = Ioffv  + Irpdpd(Irrepx,18+Ispin)
         Ioffvo = Ioffvo + Irpdpd(Irrepx,8+Ispin)
      Enddo

      Write(6,*)
      Write(6,"(a)") "Symmetry expanded OO,VV and VO matrices"
      Call output(Doo,1,Nocc,1,Nocc,Nocc,Nocc,1)
      Call output(Dvv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      Call output(Dvo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
       Return
       End


