










      Subroutine Transform_2index_pt(T1,Coo,Cvv,Work,Maxcor,List1_h,
     &                               List1_t,T1ln,Iuhf,Irrepx,CCt1,
     &                               Type)
      
      Implicit Double Precision(A-H,O-Z)
      Character*3 Type
      Integer T1ln,T1off
      Logical CCt1
      
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
c symoff.com : begin
      integer  Ioff_oo(8,2),Ioff_vv(8,2)
      common /symoff/ Ioff_oo,Ioff_vv
c symoff.com : end
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
      Dimension T1(T1ln)
      Dimension Work(Maxcor)

      Data Ione,Inull /1,0/

      Do Ispin =1, 1+Iuhf

         Ioff   = Inull
         Ndim   = Irpdpd(Irrepx,8+Ispin)
         T1off  = (Ispin-1)*Irpdpd(Irrepx,9) + Ione

         I000 = Ione
         I010 = I000 + Ndim 
         Iend = I010 + Ndim 

         If (Iend .Gt. Maxcor) Call Insmem("transform_2index",Iend,
     +                                      Maxcor)
         If (CCt1) Then
             Call Getlst(Work(I000),1,1,1,Ispin,List1_h)
         ELse
             Call Dcopy(Ndim,T1(T1off),1,Work(I000),1)
         Endif 

         do Irrep = 1, Nirrep
            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx) 
            
             Nrow = Vrt(Irrepl,Ispin)
             Ncol = Pop(Irrepr,Ispin)

             Joff = Ioff_oo(Irrepr,Ispin)
             Koff = Ioff_vv(Irrepl,Ispin)

             Call Trans(Work(I000+Ioff),Work(I010+Ioff),Coo(Joff),
     +                  Cvv(Koff),Nrow,Ncol,Type)
             
             Ioff = Ioff + Nrow*Ncol
             Joff = Joff + Ncol*Ncol
             Koff = Koff + Nrow*Nrow
         Enddo

         If (CCt1) Then
            Call Putlst(Work(I000),1,1,1,Ispin,List1_t)
         Else
            Call Dcopy(Ndim,Work(I000),1,T1(T1off),1)
         Endif 

      Enddo 

      Return
      End
            
