      Subroutine Erd_blkcp(A,Lda,M,N,B)
C B(M,N) = the leading M x N block of A(Lda,N).
      Implicit Double Precision (A-H,O-Z)
      Dimension A(Lda,N), B(M,N)
      Do J = 1, N
         Do I = 1, M
            B(I,J) = A(I,J)
         Enddo
      Enddo
      Return
      End
