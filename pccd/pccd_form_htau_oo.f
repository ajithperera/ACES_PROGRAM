










      Subroutine Pccd_form_htau_oo(Doo,Dvv,Dvo,Hoo_pq,Hoo_qp,Work,
     +                             Maxcor,Nocc,Nvrt,Noccsq,Nvrtsq,
     +                             Nvrtocc,Iuhf,Nonhf)

      Implicit DOuble Precision(A-H,O-Z)
      Logical Non_hf_terms 
      Logical Nonhf
      Logical Symmetry

      Dimension Doo(Nocc,Nocc),Dvv(Nvrt,Nvrt),Dvo(Nvrt,Nocc)
      Dimension Hoo_pq(Nocc*Nocc),Hoo_qp(Nocc*Nocc)
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

       I000 = Ione
       I010 = I000 + Nfmi(1)+Iuhf*Nfmi(2)
       I020 = I010 + Nfea(1)+Iuhf*Nfea(2)
       Iend = I020 + Nt(1)+Iuhf*Nt(2)
       If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_oo",Iend,
     +                                    Maxcor)

       Call Getrec(20,"JOBARC","DENSOO  ",Noccsq*IIntfp,Doo)
       Call Getrec(20,"JOBARC","DENSVV  ",Nvrtsq*IIntfp,Dvv)
       If (Nonhf) Call Getrec(20,"JOBARC","DENSVO  ",
     +                        Nvrtocc*IIntfp,Dvo)

       Call Analyze_Fock(Work(I000),Work(I010),Work(I020),Nfmi(1),
     +                   Nfea(1),Nt(1),Non_hf_terms)

       Call Pccd_oo(Doo,Dvv,Dvo,Hoo_pq,Work,Maxcor,Iuhf,(Non_hf_terms
     +              .or.Nonhf))

       Call Putrec(20,"JOBARC","REFGRDOO",Nfmi(1)*IIntfp,Hoo_pq)

CSSS      call pccd_check_htau("Htau_oo :",Hoo_pq,Nocc,"OO","D")
      Write(6,*)
      write(6,"(a)") "The oo contribution from ref. state"
      Write(6,"(5(1x,F15.9))") (Hoo_pq(i), i=1,Nfmi(1))
      call checksum("Hoo_pq :",Hoo_pq,Nocc*Nocc)

       Return
       End
