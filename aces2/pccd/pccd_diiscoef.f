










      Subroutine Pccd_diiscoef(Err,Augerr,Tmp,Nphys,Ndim,Info)
C
C DIIS expansion coefficients (as libr/dodiis.F), but a singular DIIS
C matrix is returned as Info.Ne.0 instead of stopping, so that the pCCD
C solver can restart the subspace (seen for stretched N2: duplicated
C error vectors after the extrapolation stagnates).
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Err(Nphys,Nphys),Augerr(Nphys+1,Nphys+2),Tmp(*)
      Parameter (Maxdim=64)
      Dimension Ipiv(Maxdim)

      If (Ndim .Ge. Maxdim) Then
         Write(6,"(a,2i6)") " @-Pccd_diiscoef: ndim too large",Ndim,
     +                      Maxdim
         Call Errex
      Endif
      Augerr(1,1)      =  0.0D0
      Augerr(1,Ndim+2) = -1.0D0
      Do J = 0, Ndim-1
         Augerr(2+J,1) = -1.0D0
         Augerr(1,2+J) = -1.0D0
         Do I = 0, Ndim-1
            Augerr(2+I,2+J) = Err(1+I,1+J)
         Enddo
         Augerr(2+J,Ndim+2) = 0.0D0
      Enddo
      Call Dgesv(Ndim+1,1,Augerr,Nphys+1,Ipiv,Augerr(1,Ndim+2),
     +           Nphys+1,Info)
      If (Info .Ne. 0) Return
      Do I = 1, Ndim
         Tmp(I) = Augerr(1+I,Ndim+2)
      Enddo
      Return
      End
