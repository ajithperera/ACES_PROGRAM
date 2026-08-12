










      SUBROUTINE CHECK_NATORB(NBAS,NBASP)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION SCR1(NBASP*NBAS),SCR2(NBAS*NBAS)

      CHARACTER*1 TYPE(2)

      DATA TYPE/"A","B"/



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




      FCTR = 1.0D0
      CALL GETREC(20,"JOBARC","AONTORBA",NBAS*NBAS*IINTFP,SCR1)
      Call checksum("CHECK_NATORB",SCR1,NBAS*NBAS)
      CALL GETREC(20,"JOBARC","AONTORBB",NBAS*NBAS*IINTFP,SCR2)
      Call checksum("CHECK_NATORB",SCR2,NBAS*NBAS)

      CALL DAXPY(NBAS*NBAS,-1.0D0,SCR1,1,SCR2,1)
      Call checksum("NA-NB       ",SCR2,NBAS*NBAS)
      WRITE(6,*)


      RETURN
      END

