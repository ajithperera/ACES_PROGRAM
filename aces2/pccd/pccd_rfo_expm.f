










      Subroutine Pccd_rfo_expm(A,U,T,Scr,N)
C U = exp(A) for real antisymmetric A (scaling and squaring, Taylor to
C 12th order), followed by modified Gram-Schmidt on the columns of U.
      Implicit Double Precision(A-H,O-Z)
      Dimension A(N,N),U(N,N),T(N,N),Scr(N,N)
      Data Zero,One /0.0D0,1.0D0/
      Anrm = Zero
      Do J = 1, N
         Do I = 1, N
            Anrm = Max(Anrm,Dabs(A(I,J)))
         Enddo
      Enddo
      Anrm = Anrm*Dble(N)
      Iscl = 0
      Fac  = One
 10   If (Anrm*Fac .Gt. 0.05D0) Then
         Fac  = 0.5D0*Fac
         Iscl = Iscl + 1
         Go to 10
      Endif
      Call Dscal(N*N,Fac,A,1)
C U = I + A + A^2/2! + ... ; T holds the current term.
      Call Dzero(U,N*N)
      Call Dzero(T,N*N)
      Do I = 1, N
         U(I,I) = One
         T(I,I) = One
      Enddo
      Do Iord = 1, 12
         Call Dgemm("N","N",N,N,N,One/Dble(Iord),T,N,A,N,Zero,Scr,N)
         Call Dcopy(N*N,Scr,1,T,1)
         Call Daxpy(N*N,One,T,1,U,1)
      Enddo
      Do Is = 1, Iscl
         Call Dgemm("N","N",N,N,N,One,U,N,U,N,Zero,Scr,N)
         Call Dcopy(N*N,Scr,1,U,1)
      Enddo
      Do J = 1, N
         Do Jp = 1, J-1
            Ov = Ddot(N,U(1,Jp),1,U(1,J),1)
            Call Daxpy(N,-Ov,U(1,Jp),1,U(1,J),1)
         Enddo
         Xn = Dsqrt(Ddot(N,U(1,J),1,U(1,J),1))
         Call Dscal(N,One/Xn,U(1,J),1)
      Enddo
      Return
      End
