         SUBROUTINE  ERD__RS_ATTENUATE
     +
     +                    ( MIJ,MKL,MIJKL,NGQP,
     +                      P,Q,PQPINV,
     +                      DOROOT,
     +
     +                                TVAL,SCALEPQ,RTS )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__RS_ATTENUATE
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : Turns the Rys quadrature of 1/r12 into the one of
C                erf(omega*r12)/r12, omega = RSOMEG (erd__rsomega.inc).
C
C                The Coulomb integral is an integral over t from 0 to
C                1 of exp(-T t**2) times a polynomial in t**2. The
C                attenuated one runs from 0 to sqrt(f) only, with
C
C                      f = omega**2 / (rho + omega**2) ,
C                    rho = p q / (p + q) ,
C
C                and the substitution t**2 = f s**2 turns it back into
C                an integral from 0 to 1. Hence, for every exponent
C                quadruple:
C
C                  1) T       ->  f T       (before the roots)
C                  2) prefact ->  sqrt(f) prefact
C                  3) roots   ->  f roots   (the roots are t**2)
C
C                Steps 1) and 2) are done with DOROOT = .FALSE., before
C                the roots and weights are calculated; step 3) with
C                DOROOT = .TRUE., after. Everything else (the VRR
C                coefficients, the 2D integrals) stays as it is.
C
C                  Input:
C
C                    MIJ(KL)      =  current # of ij (kl) primitive
C                                    index pairs
C                    MIJKL        =  MIJ * MKL
C                    NGQP         =  # of roots per quadruple
C                    P,Q          =  current MIJ (MKL) exponent sums
C                    PQPINV       =  current MIJKL values of 1/(p+q)
C                    DOROOT       =  see above
C                    TVAL         =  current MIJKL T values
C                    SCALEPQ      =  current MIJKL prefactors
C                    RTS          =  current MIJKL*NGQP roots
C
C                  Output:
C
C                    TVAL,SCALEPQ =  attenuated (DOROOT = .FALSE.)
C                    RTS          =  attenuated (DOROOT = .TRUE.)
C
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INCLUDE     'erd__rsomega.inc'

         LOGICAL     DOROOT

         INTEGER     I,IJ,KL,M,N
         INTEGER     MIJ,MKL,MIJKL,NGQP

         DOUBLE PRECISION  F,RHO,W2

         DOUBLE PRECISION  P       (1:MIJ)
         DOUBLE PRECISION  Q       (1:MKL)
         DOUBLE PRECISION  PQPINV  (1:MIJKL)
         DOUBLE PRECISION  TVAL    (1:MIJKL)
         DOUBLE PRECISION  SCALEPQ (1:MIJKL)
         DOUBLE PRECISION  RTS     (1:MIJKL*NGQP)
C
C
C------------------------------------------------------------------------
C
C
         W2 = RSOMEG * RSOMEG

         M = 0
         N = 0
         DO 100 IJ = 1,MIJ
            DO 110 KL = 1,MKL
               M = M + 1
               RHO = P (IJ) * Q (KL) * PQPINV (M)
               F = W2 / (RHO + W2)
               IF (DOROOT) THEN
                   DO 120 I = 1,NGQP
                      N = N + 1
                      RTS (N) = F * RTS (N)
  120              CONTINUE
               ELSE
                   TVAL (M) = F * TVAL (M)
                   SCALEPQ (M) = DSQRT (F) * SCALEPQ (M)
               END IF
  110       CONTINUE
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
