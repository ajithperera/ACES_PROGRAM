










      SUBROUTINE MKFOCK(WORK, MAXCOR, IUHF)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WORK(MAXCOR)
      INTEGER POP, VRT, DIRPRD
      CHARACTER*1 SP(2)
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
      COMMON /SYM/ POP(8,2),VRT(8,2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      DATA SP /'A', 'B'/
C      
      CALL GETREC(20, 'JOBARC', 'NBASTOT ', 1, NBAS)
C
      IFOCK = 1
      ISCR1 = IFOCK + NBAS*NBAS
      ISCR2 = ISCR1 + NBAS*NBAS
      IEND  = ISCR2 + 2*NBAS*NBAS
C      
      DO ISPIN = 1, 1 + MIN(1,IUHF)
         CALL GETREC(20,'JOBARC','FOCK'//SP(ISPIN),NBAS*NBAS, 
     &               WORK(IFOCK))
         CALL AO2MO2(WORK(IFOCK), WORK(IFOCK), WORK(ISCR1), 
     &               WORK(ISCR2), NBAS, NBAS, ISPIN)

         CALL MKFOCK_BLCK(WORK(IFOCK), WORK(ISCR1), WORK(ISCR2),
     &                    NOCCO(ISPIN), NVRTO(ISPIN), NBAS, ISPIN)
      ENDDO
      
      RETURN
      END
 
