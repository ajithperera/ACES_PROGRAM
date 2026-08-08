      Subroutine Write_brueckmos(Work,Maxcor,Iuhf,Nbas)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)

      Dimension Nbf4irrep(8), Nocc(8,2)

#include "machsp.com"

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
 
