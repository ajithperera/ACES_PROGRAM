










      Subroutine scrnc_ao2mo(Zao,Scr,Nbfns,Iuhf)
 
      Implicit Double Precision(A-H,O-Z)

      Dimension Zao(Nbfns*Nbfns,Nbfns*Nbfns),Scr(*)



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



    
      Ndim2 =  Nbfns * Nbfns
      Iomo  = 1
      Ivec  = Iomo  + Ndim2 
      Iscr1 = Ivec  + Ndim2
      Iscr2 = Iscr1 + Ndim2 
      Iscr3 = Iscr2 + Ndim2
      Iend  = Iscr3  	

      Call Getrec(20,"JOBARC","CMP2ZMAT",Ndim2*Iintfp,Scr(Iscr1))

      Do Ispin = 1, Iuhf + 1
      Do Ipair = 1, Ndim2 
         Call Ao2mo2(Zao(1,Ipair),Scr(Iomo),Scr(Ivec),Scr(Iend),
     +               Nbfns,Nbfns,Ispin) 

         Write(6,"(a)") "The I-1(pq,xx)"
         Call output(Scr(Iomo),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

C         Call Xgemm("N","N",Nbfns,Nbfns,Nbfns,1.0D0,Scr(Iomo),
C     +               Nbfns,Scr(Iscr1),Nbfns,0.0D0,Scr(Iscr2),
C     +               Nbfns)
C         Call Dzero(Scr(Iomo),Ndim2)
C         Call Xgemm("T","N",Nbfns,Nbfns,Nbfns,1.0D0,Scr(Iscr1),
C     +               Nbfns,Scr(Iscr2),Nbfns,0.0D0,Scr(Iomo),
C     +               Nbfns)
C         Write(6,"(a)") "The I-2(pq,xx)"
C         Call output(Scr(Iomo),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
      Enddo
      Enddo 


      Return
      End 
