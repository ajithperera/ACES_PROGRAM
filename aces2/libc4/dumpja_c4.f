










      SUBROUTINE DUMPJA_C4(IENTRY)
      IMPLICIT INTEGER (A-Z)
      LOGICAL YESNO,ISOPN
      CHARACTER*1 IENTRY
      CHARACTER*8 MARKER
      CHARACTER*80 FNAME1,FNAME2
      INTEGER LUFIL
      PARAMETER (LUFIL=75)
      PARAMETER (IBUFLN=128)
      INTEGER IBUF(IBUFLN)
      PARAMETER (LENGTH=1000)


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



      COMMON /JOBARC/ MARKER(LENGTH),LOC(LENGTH),SIZE(LENGTH),NRECS,
     &                IRECWD,IRECLN
      COMMON/RECSIZE/ISCALE
      CALL GFNAME('JOBARC  ',FNAME1,ILENGTH1)
      CALL GFNAME('JAINDX  ',FNAME2,ILENGTH2)
      IF(IENTRY.EQ.'O')THEN
       CLOSE (LUFIL)
       OPEN(UNIT=LUFIL,FILE=FNAME2(1:ILENGTH2),FORM='UNFORMATTED',
     &        STATUS='UNKNOWN')
       REWIND(LUFIL)
       WRITE(LUFIL)MARKER,LOC,SIZE,NRECS
       CLOSE(UNIT=LUFIL,STATUS='KEEP')
      ELSE
       INQUIRE(FILE=FNAME2(1:ILENGTH2),EXIST=YESNO)
       IF(YESNO)THEN
        OPEN(UNIT=LUFIL,FILE=FNAME2(1:ILENGTH2),
     &       FORM='UNFORMATTED',STATUS='OLD')
        REWIND(LUFIL)
        READ(LUFIL)MARKER,LOC,SIZE,NRECS
        CLOSE(UNIT=LUFIL,STATUS='KEEP')
       ENDIF
       INQUIRE(FILE=FNAME1(1:ILENGTH1),EXIST=YESNO,OPENED=ISOPN)
       IF(.NOT.YESNO)THEN
        IRECWD=IBUFLN
        IRECLN=IBUFLN*IINTLN/ISCALE
        OPEN(UNIT=LUFIL,FILE=FNAME1(1:ILENGTH1),
     &       FORM='UNFORMATTED',STATUS='NEW',
     &       ACCESS='DIRECT',RECL=IRECLN)
        NRECS=0
        WRITE(LUFIL,REC=1)IBUF
        DO 5 I=1,1000
         MARKER(I)='OPENSLOT'
5       CONTINUE
       ENDIF
       IF(.NOT.ISOPN)THEN
        IRECWD=IBUFLN
        IRECLN=IBUFLN*IINTLN/ISCALE
        OPEN(UNIT=LUFIL,FILE=FNAME1(1:ILENGTH1),
     &       FORM='UNFORMATTED',STATUS='OLD',
     &       ACCESS='DIRECT',RECL=IRECLN)
       ENDIF
      ENDIF
      RETURN
      END
