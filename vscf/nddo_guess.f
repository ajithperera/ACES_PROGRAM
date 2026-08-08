











      SUBROUTINE NDDO_GUESS(DENS,LDIM1,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION DENS((IUHF+1)*LDIM1)


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



      CALL GETREC(20,'JOBARC','NDDODENA',LDIM1*IINTFP,DENS)
      IF (IUHF.GT.0) THEN
      CALL GETREC(20,'JOBARC','NDDODENB',LDIM1*IINTFP,DENS(LDIM1+1))
      END IF
      RETURN
      END

