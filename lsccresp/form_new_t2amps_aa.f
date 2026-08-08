










      Subroutine Form_new_t2amps_aa(T2amps,Work,Iwork,Maxcor,Imaxcor,
     +                              List2,T2ln,T2ln_aa,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz_trn,Dissiz_sqr
      Integer T2ln,T2ln_aa,T2off
      
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
c active_space.com : begin
      Parameter(Max_xp=100)
      integer Active_oo(8,2),Active_vv(8,2)
      integer Pactive_oo(100,8,2),Pactive_vv(100,8,2)
      integer Ioff_active_oo(8,2),Ioff_active_vv(8,2)
      integer Pioff_active_oo(100,8,2),Pioff_active_vv(100,8,2)
      Double Precision Oo_threshold,Vv_threshold
      common /actvsp_info/Active_oo,Pactive_oo,Active_vv,Pactive_vv,
     +                    Ioff_active_oo,Pioff_active_oo,
     +                    Ioff_active_vv,Pioff_active_vv,
     +                    Oo_threshold,Poo_threshold,
     +                    Vv_threshold,Pvv_threshold,Eta_val(Max_xp),
     +                    E_k(Max_xp),E_ks(Max_xp)
c active_space.com: end

      Dimension Work(Maxcor)
      Dimension IWork(IMaxcor)
      Dimension T2amps(T2ln)

      Data Ione /1/

      Irrepx = Ione
      
      Do Ispin =1, 1+Iuhf

         T2off =  (Ispin-1)*T2ln_aa + Ione

         Do Irrep = 1, Nirrep

            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx) 

            Numdis_trn = Irpdpd(Irrepr,2+Ispin)
            Dissiz_trn = Irpdpd(Irrepl,Ispin)
            Numdis_sqr = Irpdpd(Irrepr,20+Ispin)
            Dissiz_sqr = Irpdpd(Irrepl,18+Ispin)
            
            I000 = Ione
            I010 = I000 + Numdis_sqr*Dissiz_sqr
            I020 = I010 + Numdis_sqr*Dissiz_sqr
            I030 = I020 + Numdis_sqr*Dissiz_sqr
            Iend = I030 + Numdis_sqr*Dissiz_sqr
            If (Iend .Gt. Maxcor) Call Insmem("form_new_amps_aa",
     +                                         Iend,Maxcor)

C T2(A<B,I<J)
            Call Getlst(Work(I000),1,Numdis_trn,1,Irrepr,
     +                  List2+Ispin)
            Call Dcopy(Dissiz_trn*Numdis_trn,T2amps(T2off),Ione,
     +                 Work(I020),Ione)


C T2(A<B,I<J)->T2(A<B,IJ)

            Call Symexp(Irrepr,Pop(1,Ispin),Dissiz_trn,Work(I000))
            Call Symexp(Irrepr,Pop(1,Ispin),Dissiz_trn,Work(I020))

C T2(A<B,IJ)->T2(IJ,A<B)

            Call Transp(Work(I000),Work(I010),Numdis_sqr,
     +                  Dissiz_trn)
            Call Transp(Work(I020),Work(I030),Numdis_sqr,
     +                  Dissiz_trn)

C T2(IJ,A<B)->T2(IJ,AB)

            Call Symexp(Irrepl,Vrt(1,Ispin),Numdis_sqr,Work(I010))
            Call Symexp(Irrepl,Vrt(1,Ispin),Numdis_sqr,Work(I030))

C T2(IJ,AB)->T2(AB,IJ)

            Call Transp(Work(I010),Work(I000),Dissiz_sqr,
     +                  Numdis_sqr)
            Call Transp(Work(I030),Work(I020),Dissiz_sqr,
     +                  Numdis_sqr)

            Call Form_vvoo_aa(Work(I000),Work(I020),Work(I030),
     +                        Iwork,Imaxcor,Numdis_sqr,Dissiz_sqr,
     +                        Irrepr,Irrepl,Ispin)

C T2(AB,IJ)->T2(A<B,IJ)

            Call Sqsym(Irrepl,Vrt(1,Ispin),Dissiz_trn,Dissiz_sqr,
     +                 Numdis_sqr,Work(I010),Work(I000))
            Call Sqsym(Irrepl,Vrt(1,Ispin),Dissiz_trn,Dissiz_sqr,
     +                 Numdis_sqr,Work(I030),Work(I020))

C T2(A<B,IJ)->T2(IJ,A<B)

            Call Transp(Work(I010),Work(I000),Numdis_sqr,
     +                   Dissiz_trn)
            Call Transp(Work(I030),Work(I020),Numdis_sqr,
     +                   Dissiz_trn)

C T2(IJ,AB)->T2(I<J,A<B)

            Call Sqsym(Irrepr,Pop(1,Ispin),Numdis_trn,Numdis_sqr,
     +                 dissiz_trn,Work(I010),Work(I000))
            Call Sqsym(Irrepr,Pop(1,Ispin),Numdis_trn,Numdis_sqr,
     +                 dissiz_trn,Work(I030),Work(I020))

C T2(I<J,A<B)->T2(A<B,I<J)

            Call Transp(Work(I010),Work(I000),Dissiz_trn,
     +                  Numdis_trn)
            Call Transp(Work(I030),Work(I020),Dissiz_trn,
     +                  Numdis_trn)

            Call Putlst(Work(I000),1,Numdis_trn,1,Irrepr,
     +                  List2+Ispin)

            T2off = T2off + Dissiz_trn*Numdis_trn

        Enddo
      Enddo 

      Return
      End
            
