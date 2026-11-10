      Subroutine Vee_bd_expnd2a(Wpack,Wfull,Ndim)
C
C Expand a row-packed lower triangle, Wpack(i(i-1)/2+j) = W(i,j)
C for j <= i, into the antisymmetric square matrix W (vprops packing
C of the AM and MQ integrals); the diagonal is zero. Same as nmr/expnd2a, kept local to
C vee so that xvee does not depend on libnmr.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Wpack(*),Wfull(Ndim,Ndim)
      Ithru = 0
      Do I = 1, Ndim
         Do J = 1, I
            Ithru = Ithru + 1
            Wfull(I,J) =  Wpack(Ithru)
            Wfull(J,I) = -Wpack(Ithru)
         Enddo
         Wfull(I,I) = 0.0D0
      Enddo
      Return
      End
