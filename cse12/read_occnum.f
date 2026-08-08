










      SUBROUTINE READ_OCCNUM(DOCC,MAXCOR,ACT_MIN_A,ACT_MAX_A,ACT_MIN_B,
     +                       ACT_MAX_B,NBASIS,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      DIMENSION DOCC(MAXCOR), NSUM(2), NOC_ORB(2)
      CHARACTER*80 FNAME, Blank
      INTEGER ACT_MIN_A,ACT_MIN_B,ACT_MAX_A,ACT_MAX_B


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
      Write(6,"(a,a)") " The occupation numbers are read",
     &                " from OCCNUMS file"
      Write(6,*)

      IUNIT = 5
      OPEN(UNIT=IUNIT,FILE="OCCNUMS",FORM="FORMATTED")

      READ(IUNIT, "(80a)") Blank
      READ(IUNIT, "(80a)") Blank
C
      READ(IUNIT,10,END=19) NOC_ORB(1)
      READ(IUNIT,10,END=19) NOC_ORB(2)

      READ(IUNIT, "(80a)") Blank

      READ(IUNIT,10,END=19) NSUM(1)
      READ(IUNIT,10,END=19) NSUM(2)
C
      READ(IUNIT, "(80a)") Blank

      DO ISPIN =1, (IUHF+1)
         ISTART = 0
         DO IBF = 1, NOC_ORB(ISPIN)
            ISTART = ISTART + 1
            READ(IUNIT,*,END=19) DOCC_NUM
            DOCC(ISTART) = DOCC_NUM 
         ENDDO

         DO IBF = 1, (NSUM(ISPIN)-NOC_ORB(ISPIN))
            ISTART = ISTART + 1
            DOCC(ISTART) = 0.0D0
         ENDDO
         READ(IUNIT, "(80a)") Blank
         IF (ISPIN .EQ. 1) THEN
            CALL PUTREC(20,"JOABRC","ORB_OCCA",NBASIS,DOCC)
      Write(6,"(a)")" The Alpha and Beta occupation numbers"
      write(6,*)
      Write(6,"(6F10.5)")(DOCC(I), I=1, NBASIS)
         ELSE 
            CALL PUTREC(20,"JOABRC","ORB_OCCB",NBASIS,DOCC)
      write(6,*)
      Write(6,"(6F10.5)")(DOCC(I), I=1, NBASIS)
      write(6,*)
         ENDIF
      ENDDO 

C Read the active info blocks.

      READ(IUNIT,*,END=19) ACT_MIN_A
      READ(IUNIT,*,END=19) ACT_MAX_A
      READ(IUNIT, "(80a)") Blank
      READ(IUNIT,*,END=19) ACT_MIN_B
      READ(IUNIT,*,END=19) ACT_MAX_B

 10   FORMAT(16I5)
 19   CLOSE(IUNIT)
C
C Copy the occupation in the case of RHF
C
      IF (IUHF .EQ. 0) 
     &    CALL PUTREC(20,"JOABRC","ORB_OCCB",NBASIS,DOCC)

      RETURN
      END
                
