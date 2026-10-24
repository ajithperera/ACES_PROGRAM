










      SUBROUTINE EVEC_SHIFT(EIGEN, VEC, SCR, NOPT)
C 
C This routine moves the 3N-6(or 5) translational and rotational
C modes to the end (ie. This eigenvalues and vectors that have
C zero eigenvalues). Ajith Perera, 11/2012 (the date is when
C the comment was added, the routine was written 3 or 4 years 
C prior to that).  

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION EIGEN(NOPT, NOPT), VEC(NOPT, NOPT), SCR(NOPT, NOPT)
C
      Write(6,"(a)") "Unsorted Eigen values and vectors"
      call output(eigen, 1, nopt, 1, nopt, nopt, nopt, 1)
      call output(vec, 1, nopt, 1, nopt, nopt, nopt, 1)

C Move each zero mode to the end by adjacent swaps. Scanning from
C the right keeps the order of both the zero and non-zero modes, and
C needs no work space (the old version kept the vectors in rows
C 1-NTR of SCR and the eigenvalues in row 8, which overruns SCR for
C NOPT < 8, and it used an uninitialized ISHIFT when no zero modes
C were present, e.g. constrained runs where the TR modes are not
C projected). SCR is no longer used.
C
      NEND=NOPT
      DO K=NOPT,1,-1
         IF (ABS(EIGEN(K,K)).LT.1.0D-09) THEN
            DO J=K,NEND-1
               TMP=EIGEN(J,J)
               EIGEN(J,J)=EIGEN(J+1,J+1)
               EIGEN(J+1,J+1)=TMP
               CALL DSWAP(NOPT, VEC(1,J), 1, VEC(1,J+1), 1)
            ENDDO
            NEND=NEND-1
         ENDIF
      ENDDO
C
      Write(6,"(a)") "Sorted Eigen values and vectors"
      call output(eigen, 1, nopt, 1, nopt, nopt, nopt, 1)
      call output(vec, 1, nopt, 1, nopt, nopt, nopt, 1)
      
      RETURN
      END


