










      Subroutine Write_brueckmos(Work,Maxcor,Iuhf,Nbas)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)

      Dimension Nbf4irrep(8), Nocc(8,2)



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




      I000 = 1
      I010 = I000 + Nbas * Nbas 
      I020 = I010 + (Iuhf-1) * Nbas

      Call Getrec(20,"JOBARC","SCFEVCA0",Nbas*Nbas*IINTFP,
     &            Work(I000))
      If (Iuhf .NE. 0)
     &    Call Getrec(20,"JOBARC","SCFEVCA0",Nbas*Nbas*IINTFP,
     &                Work(I010))

      Call Getrec(20, "JOBARC", 'NIRREP  ',   1,      Nirrep)
      Call Getrec(20, "JOBARC", 'NUMBASIR',   Nirrep, Nbf4irrep)
      Call Getrec(20, "JOBARC", 'OCCUPYA0',   Nirrep, Nocc(1,1))
      If (Iuhf .NE. 0) Then
         Call Getrec(20, "JOBARC", 'OCCUPYB0',   Nirrep, Nocc(1,2)) 
      Else
         Call Icopy(8, Nocc(1,1),1,Nocc(1,2),1)   
      Endif 
C
      Call PutBrueckmos(Work(I000), Nbas*Nbas, Iuhf, Nirrep, 
     &                  Nbf4irrep,Nocc)

      Return
      End 
 
