










        subroutine GetSDInfo(nocca,noccb,nvrta,nvrtb,nbas,vecDim,
     +                       Naobfns)
        integer nocc,nvirt,nbas,vecDim
        integer Naobfns,Nbfirr
        integer Ndrop,Ione
        Dimension Nbfirr(8)
        
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end


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



c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

        Ione=1
        Nbas=Nocco(1)+Nvrto(1)
        Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)

        If (Idrop .Ne. 0) Then
         Call Getrec(20,"JOBARC","NUMDROPA",1,Ndrop)
         Write(6,"(2a,I4)"),'@SDinfo number of dropped basis', 
     +                      ' fxn:', Ndrop
CSSS         Nbas = (Nbas - Ndrop)
        Endif

! Build and return all relevant quantities

        vecDim=Nbas*Nbas
        nocca=Nocco(1)
        noccb=Nocco(2)
        nvrta=Nvrto(1)
        nvrtb=Nvrto(2)
        Naobfns=nocc+nvrt

        end subroutine
