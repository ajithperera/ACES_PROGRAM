      Subroutine Pccd_antisymmetrize(Htau_pq,Work,Maxcor,Nbas,Nocc,
     +                               Nvrt)

       Implicit Double Precision(A-H,O-Z)
       Logical Sym_packed 

       Dimension Htau_pq(Nbas,Nbas)
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

      Data Ione, Izero /1,0/

      Irrepx  = Ione
      Length  = Nbas*Nbas

      Sym_packed = .False. 
      If (Sym_packed) Then
         Ndim_oo  = Irpdpd(Irrepx,21)
         Ndim_vv  = Irpdpd(Irrepx,19)
         Ndim_vo  = Irpdpd(Irrepx,9)
      Else
         Ndim_oo = Nocc*Nocc
         Ndim_vv = Nvrt*Nvrt
         Ndim_vo = Nocc*Nvrt
      Endif 

      Lenoo    = Ndim_oo
      Lenvv    = Ndim_vv
      Lenvo    = Ndim_vo

      I000 = Ione
      I010 = I000 + Ndim_oo
      I020 = I010 + Ndim_vv
      I030 = I020 + Ndim_vo
      Iend = I030 + Ndim_vo
      Memleft = Maxcor - Iend
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_antisymmetrize",Iend,
     +                                   Maxcor)

      Call Pccd_frmblocks(Htau_pq,Work(I000),Work(I010),Work(I020),
     +                    Work(I030),Work(Iend),Memleft,Lenoo,Lenvv,
     +                    Lenvo,Nocc,Nvrt,Nbas)

C Store the OV and VO gradients befor antsymmetrize them.      

      Call Putrec(20,"JOBARC","VO_GRADS",Ndim_vo*Iintfp,Work(I020))
      Call Putrec(20,"JOBARC","OV_GRADS",Ndim_vo*Iintfp,Work(I030))

      I040 = Iend
      I050 = I040 + Ndim_oo
      I060 = I050 + Ndim_vv
      I070 = I060 + Ndim_vo 
      Iend = I070 + Ndim_vo 
      Memleft = Maxcor - Iend
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_antisymmetrize",Iend,
     +                                   Maxcor)

      Call Pccd_asymmoovv(Work(I000),Work(I010),Work(I040),Work(I050),
     +                    Lenoo,Lenvv,Nocc,Nvrt)
      Call Pccd_asymmovvo(Work(I020),Work(I030),Work(I060),Work(I070),
     +                    Lenvo,Nocc,Nvrt)

      Call Pccd_putblocks(Htau_pq,Work(I000),Work(I010),Work(I020),
     +                    Work(I030),Work(Iend),Memleft,
     +                    Lenoo,Lenvv,Lenvo,Nocc,Nvrt,Nbas,"Vo_like",
     +                    Sym_packed)

      Return 
      End 


       
