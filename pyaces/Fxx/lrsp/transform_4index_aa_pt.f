










      Subroutine Transform_4index_aa_pt(T2amps,Coo,Cvv,Work,Maxcor,
     +                                  List2_h,List2_t,T2ln,T2ln_aa, 
     +                                  Iuhf,Irrepx,Cct2,Type)
      
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz_trn,Dissiz_sqr
      Integer T2ln,T2ln_aa,T2off
      Character*3 Type
      Logical CCt2
      
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

      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Work(Maxcor)
      Dimension T2amps(T2ln)

      Data Ione /1/

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
            Iend = I010 + Numdis_sqr*Dissiz_sqr
            If (Iend .Gt. Maxcor) Call Insmem("transform_4index_aa",
     +                                         Iend,Maxcor)
C T2(A<B,I<J)
            If (Cct2) Then
                Call Getlst(Work(I000),1,Numdis_trn,1,Irrepr,
     +                      List2_h+Ispin)
            Else
                Call Dcopy(Dissiz_trn*Numdis_trn,T2amps(T2off),Ione,
     +                      Work(I000),Ione)
            Endif

C T2(A<B,I<J)->T2(A<B,IJ)

            Call Symexp(Irrepr,Pop(1,Ispin),Dissiz_trn,Work(I000))

C T2(A<B,IJ)->T2(IJ,A<B)

            Call Transp(Work(I000),Work(I010),Numdis_sqr,
     +                  Dissiz_trn)

            Do I = 1,  Dissiz_trn
               
               Ioff =  (I-1)*Numdis_sqr 

               Call Trans_oo_aa(Work(I010+Ioff),Work(I000+Ioff),Coo,
     +                          Numdis_sqr,Irrepr,Iuhf,Ispin,Type)
            Enddo 

C T2(IJ,A<B)->T2(IJ,AB)

            Call Symexp(Irrepl,Vrt(1,Ispin),Numdis_sqr,Work(I010))

C T2(IJ,AB)->T2(AB.IJ)

            Call Transp(Work(I010),Work(I000),Dissiz_sqr,
     +                  Numdis_sqr)

            Do I = 1, Numdis_sqr

               Joff = (I-1)*Dissiz_sqr 
               Call Trans_vv_aa(Work(I000+Joff),Work(I010+Joff),Cvv,
     +                          Dissiz_sqr,Irrepl,Iuhf,Ispin,Type)
            Enddo

C T2(AB,IJ)->T2(A<B,IJ)
            Call Sqsym(Irrepl,Vrt(1,Ispin),Dissiz_trn,Dissiz_sqr,
     +                 Numdis_sqr,Work(I010),Work(I000))

C T2(A<B,IJ)->T2(IJ,A<B)

            Call Transp(Work(I010),Work(I000),Numdis_sqr,
     +                   Dissiz_trn)

C T2(IJ,AB)->T2(I<J,A<B)

            Call Sqsym(Irrepr,Pop(1,Ispin),Numdis_trn,Numdis_sqr,
     +                 dissiz_trn,Work(I010),Work(I000))

C T2(I<J,A<B)->T2(A<B,I<J)

            Call Transp(Work(I010),Work(I000),Dissiz_trn,
     +                  Numdis_trn)

            If (Cct2) Then
                Call Putlst(Work(I000),1,Numdis_trn,1,Irrepr,
     +                      List2_t+Ispin)
            Else
               Call Dcopy(Dissiz_trn*Numdis_trn,Work(I000),Ione,
     +                      T2amps(T2off),Ione)
            Endif 

            T2off = T2off + Dissiz_trn*Numdis_trn
        Enddo
      Enddo 

      Return
      End
            
