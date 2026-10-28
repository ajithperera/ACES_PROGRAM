      Subroutine Erd_minv(A,Ai,N,Ipiv,Work,Lwork)
C Ai = inverse of the N x N matrix A (A is kept), with LAPACK.
      Implicit Double Precision (A-H,O-Z)
      Dimension A(N,N), Ai(N,N), Work(Lwork)
      Integer Ipiv(N)
      Call Dcopy(N*N,A,1,Ai,1)
      Call Dgetrf(N,N,Ai,N,Ipiv,Info)
      If (Info .Eq. 0) Call Dgetri(N,Ai,N,Ipiv,Work,Lwork,Info)
      If (Info .Ne. 0) Then
         Write(6,"(a,i6)") " @Erd_minv: inversion failed, info ", Info
         Call Errex
      Endif
      Return
      End
