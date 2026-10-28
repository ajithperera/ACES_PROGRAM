      Subroutine Erg_btr(N,B,Dso,D,Scr,Beta)
C
C D = B Dso B(T) + Beta D.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension B(N,N), Dso(N,N), D(N,N), Scr(N,N)

      Call Dgemm('N','N',N,N,N,1.0D0,B,N,Dso,N,0.0D0,Scr,N)
      Call Dgemm('N','T',N,N,N,1.0D0,Scr,N,B,N,Beta,D,N)

      Return
      End
