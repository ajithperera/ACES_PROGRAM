










      Subroutine Erg_pinv(N,Ng,T,A,P,Ipiv)
C
C P = (T(T) T)**-1 T(T), the pseudoinverse of T (Ng x N, full column
C rank).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension T(Ng,N), A(N,N), P(N,Ng)
      Integer Ipiv(N)

      Call Dgemm('T','N',N,N,Ng,1.0D0,T,Ng,T,Ng,0.0D0,A,N)
      Do J = 1, Ng
         Do I = 1, N
            P(I,J) = T(J,I)
         Enddo
      Enddo
      Call Dgesv(N,Ng,A,N,Ipiv,P,N,Info)
      If (Info .Ne. 0) Then
         Write(6,"(a,i6)") " @Erg_pinv: Dgesv failed, info ", Info
         Call Errex
      Endif

      Return
      End
