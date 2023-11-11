      Subroutine Pccd_form_htau_ov(Doo,Dvv,Dvo,Hov_pq,Hov_qp,Work,
     +                             Maxcor,Nocc,Nvrt,Noccsq,Nvrtsq,
     +                             Nvrtocc,Iuhf,Nonhf)

      Implicit DOuble Precision(A-H,O-Z)
      Logical Nonhf
      Logical Symmetry

      Dimension Doo(Nocc,Nocc),Dvv(Nvrt,Nvrt),Dvo(Nvrt,Nocc)
      Dimension Hov_pq(Nocc*Nvrt)
      Dimension Hov_qp(Nocc*Nvrt)
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

       Common /Sym/Symmetry

       Data Ione /1/

       Call Pccd_ov(Doo,Dvv,Dvo,Hov_pq,Work,Maxcor,Iuhf,Nonhf)
       Call Putrec(20,"JOBARC","REFGRDOV",Nt(1)*IIntfp,Hov_pq)

      Write(6,*)
      write(6,"(a)") "The X(ov) contribution"
      Write(6,"(5(1x,F15.9))") (Hov_pq(i), i=1,Nt(1))
      call checksum("Hov_pq :",Hov_pq,Nocc*Nvrt)

       Return
       End
