










      SUBROUTINE SALC2ATM(SPNC_TENSOR, SMTRN, SCR, NTPERT)
C
C THE COLUMNS OF SPNC_TENSOR ARE THE SYMMETRY-ADAPTED PERTURBATIONS
C VCCEH SOLVED FOR, THE ROWS THE ATOMIC PROPERTY INTEGRALS, I.E.
C SPNC_TENSOR = S*J(SALC). VCCEH (JJSPNCUP) FORMS J = S*J(SALC)*S^T,
C SO MULTIPLY BY S^T FROM THE RIGHT TO GET THE SAME ATOM-ATOM TENSOR.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL JFC, JPSO, JSD
C


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



C
      DIMENSION SPNC_TENSOR(NTPERT, NTPERT), SMTRN(NTPERT, NTPERT),
     &          SCR(NTPERT, NTPERT)
      COMMON /SPN_SPN_TYPE/ JFC, JPSO, JSD
C
      DATA AZERO, ONE /0.0D0, 1.0D0/
C
      IF (JFC) THEN
         CALL GETREC(20, 'JOBARC', 'SMTRNFC ', NTPERT*NTPERT*IINTFP,
     &               SMTRN)
      ELSE IF (JPSO) THEN
         CALL GETREC(20, 'JOBARC', 'SMTRNPSO', NTPERT*NTPERT*IINTFP,
     &               SMTRN)
      ELSE IF (JSD) THEN
         CALL GETREC(20, 'JOBARC', 'SMTRNSD ', NTPERT*NTPERT*IINTFP,
     &               SMTRN)
      ENDIF
C
      CALL XGEMM('N', 'T', NTPERT, NTPERT, NTPERT, ONE, SPNC_TENSOR,
     &           NTPERT, SMTRN, NTPERT, AZERO, SCR, NTPERT)
      CALL DCOPY(NTPERT*NTPERT, SCR, 1, SPNC_TENSOR, 1)
C
      RETURN
      END
