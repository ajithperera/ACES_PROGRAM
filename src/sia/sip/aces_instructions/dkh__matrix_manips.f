C------------------------------------------------------------------------
C  OPERATION   : DKH__MATRIX_MANIPS
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DKH__SQ_PRINT
C                DKH__MATMUL_NN
C                DKH__MATMUL_NT
C                DKH__MATMUL_TN
C                DKH__MATMULD
C                DKH__MATMULM
C                DKH__MATDIVE
C                DKH__MATLFTE
C                DKH__MATRHTE
C                DKH__MATRIXADD
C
C  DESCRIPTION : This is a collection of subroutines used specifically
C                in the Douglas-Kroll-Hess transformations.
C
C                It is essentially DGEMM calls, since all the DKH
C                transformation is is a set of matrix multiplies.
C
C                Rather, it is XGEMM calls, because it was designed
C                with the ACES3 software package in mind.  It is
C                trivial to change though.
C
C
C------------------------------------------------------------------------
C
C
C             ...Print a square matrix!
C
C
         SUBROUTINE  DKH__SQ_PRINT (N, A, string)

         IMPLICIT    NONE

         INTEGER     I,J,K,L,N
         INTEGER     NUMCOL,ISTART
         INTEGER     LASTINDX,FULLLOOPS

         DOUBLE PRECISION  A (1:N,1:N)
         DOUBLE PRECISION  ZERO,TOL

         PARAMETER ( NUMCOL = 5 )
         PARAMETER ( ZERO   = 0.0D0 )
         PARAMETER ( TOL    = 1.0D-08 )

         DOUBLE PRECISION  SCR (1:NUMCOL)

         CHARACTER(*)  STRING

 5000    FORMAT (I17, 4I14)
 5010    FORMAT (I7, 5D14.6)

         WRITE (6,'(/,A)') STRING

         LASTINDX  = MOD (N,NUMCOL)
         FULLLOOPS = (N - LASTINDX) / NUMCOL

         ISTART = 1
         DO I = 1, FULLLOOPS
            WRITE (6,5000) ( (ISTART + J), J = 0, NUMCOL-1 )
            DO K = 1, N
               DO L = ISTART, ISTART+NUMCOL-1
                  IF (DABS (A (K,L)) .GT. TOL) THEN
                     SCR (L-ISTART+1) = A (K,L)
                  ELSE
                     SCR (L-ISTART+1) = ZERO
                  ENDIF
               END DO
               WRITE (6,5010) K, ( SCR (L), L = 1, NUMCOL )
            END DO
            ISTART = ISTART + NUMCOL
         END DO

         IF (LASTINDX .GT. 0) THEN
            WRITE (6,5000) ( (ISTART + J), J = 0, LASTINDX-1)
            DO K = 1, N
               DO L = ISTART, ISTART+LASTINDX-1
                  IF (DABS (A (K,L)) .GT. TOL) THEN
                     SCR (L-ISTART+1) = A (K,L)
                  ELSE
                     SCR (L-ISTART+1) = ZERO
                  ENDIF
               END DO
               WRITE (6,5010) K, ( SCR (L), L=1, LASTINDX )
            END DO
         END IF

         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C             ...Regular matrix multiply!
C
C                   C = alpha A * B + beta * C
C
C
         SUBROUTINE  DKH__MATMUL_NN (N, A, B, C, ALPHA, BETA)

         IMPLICIT    NONE

         INTEGER     N

         DOUBLE PRECISION  ALPHA,BETA
         DOUBLE PRECISION  A (1:N,1:N)
         DOUBLE PRECISION  B (1:N,1:N)
         DOUBLE PRECISION  C (1:N,1:N)
         
         CALL  XGEMM
     +
     +               ( 'N','N',
     +                  N,N,N,ALPHA,
     +                  A,N,B,N,
     +
     +                       BETA,C,N )
     +
     +
         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C              ...Perform a matrix times the transpose of another!
C
C                            T
C                 C = a A * B  + b C
C
C
         SUBROUTINE  DKH__MATMUL_NT (N,A,B,C,ALPHA,BETA)

         IMPLICIT    NONE

         INTEGER     N

         DOUBLE PRECISION  ALPHA,BETA
         DOUBLE PRECISION  A (1:N,1:N)
         DOUBLE PRECISION  B (1:N,1:N)
         DOUBLE PRECISION  C (1:N,1:N)

         CALL  XGEMM
     +
     +               ( 'N','T',
     +                  N,N,N,ALPHA,
     +                  A,N,B,N,
     +
     +                       BETA,C,N )
     +
     +
         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C              ...Perform a transpose matrix times another matrix!
C
C                        T
C                 C = a A * B  + b C 
C
C
         SUBROUTINE  DKH__MATMUL_TN (N,A,B,C,ALPHA,BETA)

         IMPLICIT    NONE

         INTEGER     N

         DOUBLE PRECISION  ALPHA,BETA  
         DOUBLE PRECISION  A (1:N,1:N)
         DOUBLE PRECISION  B (1:N,1:N)
         DOUBLE PRECISION  C (1:N,1:N)

         CALL  XGEMM
     +
     +               ( 'T','N',
     +                  N,N,N,ALPHA,
     +                  A,N,B,N,
     +
     +                       BETA,C,N )
     +
     +
         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C
         SUBROUTINE  DKH__MATMULD (N,A,B,C,SCR,ALPHA,BETA,TT,RR)

         IMPLICIT    NONE

         INTEGER     I,J,N

         DOUBLE PRECISION  ALPHA,BETA
         DOUBLE PRECISION  HALF
         PARAMETER ( HALF = 0.50 D0 )

         DOUBLE PRECISION  TT (1:N)
         DOUBLE PRECISION  RR (1:N)

         DOUBLE PRECISION  A (1:N,1:N)
         DOUBLE PRECISION  B (1:N,1:N)
         DOUBLE PRECISION  C (1:N,1:N)

         DOUBLE PRECISION  SCR (1:N,1:N)

         DO I = 1, N
         DO J = 1, N
            SCR (I,J) = A (I,J) * (HALF / (TT (J) * RR (J) * RR (J)))
         END DO
         END DO

         CALL  XGEMM
     +
     +               ( 'N','N',
     +                  N,N,N,ALPHA,
     +                  SCR,N,B,N,
     +
     +                         BETA,C,N )
     +
     +
         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C
         SUBROUTINE  DKH__MATMULM (N,A,B,C,SCR,ALPHA,BETA,TT,RR)

         IMPLICIT    NONE

         INTEGER     I,J,N

         DOUBLE PRECISION  ALPHA,BETA
         DOUBLE PRECISION  TWO
         PARAMETER  ( TWO = 2.0 D0 )

         DOUBLE PRECISION  TT (1:N)
         DOUBLE PRECISION  RR (1:N)

         DOUBLE PRECISION  A  (1:N,1:N)
         DOUBLE PRECISION  B  (1:N,1:N)
         DOUBLE PRECISION  C  (1:N,1:N)

         DOUBLE PRECISION  SCR (1:N,1:N)

         DO I = 1, N
         DO J = 1, N
            SCR (I,J) = A (I,J) * (TWO * TT (J) * RR (J) * RR (J))
         END DO
         END DO

         CALL  XGEMM
     +
     +               ( 'N','N',
     +                  N,N,N,ALPHA,
     +                  SCR,N,B,N,
     +
     +                         BETA,C,N )
     +
     +
         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
         SUBROUTINE  DKH__MATDIVE (N,A,E)

         IMPLICIT    NONE

         INTEGER     I,J,N

         DOUBLE PRECISION  E (1:N)
         DOUBLE PRECISION  A (1:N,1:N)

         DO I = 1, N
         DO J = 1, N
            A (I,J) = A (I,J) / ( E(I) + E(J) )
         END DO
         END DO

         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C
         SUBROUTINE  DKH__MATLFTE (N,ALPHA,BETA,A,B,SCR,E)

         IMPLICIT    NONE

         INTEGER     I,IJ,N

         DOUBLE PRECISION  ALPHA,BETA

         DOUBLE PRECISION  E (1:N)

         DOUBLE PRECISION  A (1:N,1:N)
         DOUBLE PRECISION  B (1:N,1:N)

         DOUBLE PRECISION  SCR (1:N,1:N)

         CALL  ZERO (SCR, N*N)
         DO I = 1, N
            SCR (I,I) = E (I)
         END DO

         CALL  XGEMM
     +
     +               ( 'N','N',
     +                  N,N,N,ALPHA,
     +                  SCR,N,A,N,
     +
     +                         BETA,B,N )
     +
     +
         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C
         SUBROUTINE  DKH__MATRHTE (N,ALPHA,BETA,A,B,SCR,E)

         IMPLICIT    NONE

         INTEGER     I,J,IJ,N

         DOUBLE PRECISION  ALPHA,BETA 

         DOUBLE PRECISION  E (1:N)

         DOUBLE PRECISION  A (1:N,1:N)
         DOUBLE PRECISION  B (1:N,1:N)

         DOUBLE PRECISION  SCR (1:N,1:N)

         CALL  ZERO (SCR, N*N)
         DO I = 1, N
            SCR (I,I) = E (I)
         END DO

         CALL  XGEMM
     +
     +               ( 'N','N',
     +                  N,N,N,ALPHA,
     +                  A,N,SCR,N,
     +
     +                         BETA,B,N )
     +
     +
         RETURN
         END
C
C
C-----------------------------------------------------------------------
C
C
C
         SUBROUTINE  DKH__MATRIXADD (N,ALPHA,A,BETA,B,C)

         IMPLICIT    NONE

         INTEGER     I,N

         DOUBLE PRECISION  ALPHA,BETA

         DOUBLE PRECISION  A (1:N)
         DOUBLE PRECISION  B (1:N)
         DOUBLE PRECISION  C (1:N)

         DO I = 1, N
            C (I) = ALPHA * A (I) + BETA * B (I)
         END DO

         RETURN
         END
C
C
C             ...ready!
C
C
C-----------------------------------------------------------------------
