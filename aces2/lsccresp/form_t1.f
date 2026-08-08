










      Subroutine Form_t1(T1,Work,Maxcor,Pop,Vrt,Ispin,T1ln,T1ln_aa)

      Implicit Double Precision(A-H,O-Z)

      Integer Pop,Vrt 
      Integer I,A,T1ln,T1ln_aa

      Dimension Pop(8),Vrt(8)
      Dimension T1(T1ln)
      Dimension Work(Maxcor)

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




      Data Ione,Izero,Threshold /1,0,1.0D-09/
      Data Dzero /0.0D0/

      Ndim  = Izero 
      Ndimo = Izero 

      Do Irrep = 1, Nirrep 
         Ndim  = Ndim  + Pop(Irrep) + Vrt(Irrep)
         Ndimo = Ndimo + Pop(Irrep) 
      Enddo 
      
      I000 = Ione 
      If (Ispin .Eq. 1) Then
         Call Getrec(20,"JOBARC","SCFEVALA",Ndim*Iintfp,Work(I000))
      Else
         Call Getrec(20,"JOBARC","SCFEVALB",Ndim*Iintfp,Work(I000))
      Endif  
      Ind   = Izero 
      Indi  = Izero
      Inda0 = Izero 

      Do Irrep = 1, Nirrep
         Nrow = Vrt(Irrep)
         Ncol = Pop(Irrep)

         Do I = 1, Ncol
            Indi = Indi + Ione
            Inda = Inda0
            Do A = 1, Nrow 
               Inda = Inda + Ione
               Ind  = Ind + Ione 
               E = Work(Indi) - Work(Ndimo+Inda)
               If (Abs(E) .Gt. Threshold) Then
                  T1(Ind) = T1(Ind)/E
               Else 
                  T1(Ind) = Dzero 
               Endif 
            Enddo
         Enddo
         Inda0 = Inda0 + Nrow 
      Enddo 

      Return
      End


