













































































































































































































      SUBROUTINE DRIVE_RPA(ICORE,MAXCOR,IUHF)
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      LOGICAL RCCD_EOM
C
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /EOM_METHOD/RCCD_EOM
C
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end


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




C If any of the r/dr-CCD EOM flags are set and RPA is also run, then
C write the RPA T2 vectors to lists and set a flag to inform the r/dr-CCD
C codes that T2 is already available.

      RCCD_EOM = ((IFLAGS(87) .EQ. 13) .OR.
     +            (IFLAGS(87) .EQ. 14) .OR.
     +            (IFLAGS(87) .EQ. 15) .OR.
     +            (IFLAGS(87) .EQ. 16) .OR.
     +            (IFLAGS(87) .EQ. 17) .OR.
     +            (IFLAGS(87) .EQ. 18) .OR.
     +            (IFLAGS(87) .EQ. 19) .OR.
     +            (IFLAGS(87) .EQ. 20))

      IF (IUHF.EQ.0)THEN
C
          NSIZ=IRPDPD(1,ISYTYP(1,25))
C
          I000=1
          I010=I000+2*IINTFP*NSIZ
          MXCOR=MAXCOR-I010+1
          CALL DO_RHF_RPA(ICORE(I010),MXCOR,ICORE(I000),NSIZ,IUHF)
C
      ELSEIF (IUHF.EQ.1) THEN

          NSIZ=IRPDPD(1,ISYTYP(1,19))+IRPDPD(1,ISYTYP(1,20))

          I000=1
          I010=I000+2*IINTFP*NSIZ
          MXCOR=MAXCOR-I010+1
          CALL DO_UHF_RPA(ICORE(I010),MXCOR,ICORE(I000),NSIZ,IUHF)

      ENDIF

      IF (RCCD_EOM) CALL PROCESS_RPAT2(ICORE,MAXCOR,IUHF)
C
CSS      I000 = 1
CSS      call modf_amps(ICORE(I000),MAXCOR,IUHF,0,.FALSE.,'T',"SAVE",
CSS     &               "AMPS")
C
C Change the refrence to RHF to compute the T1 energy contribution 
C corresponding to direct RPA.
C
      iflags(11) = 1
      Iuhf = 0
CSSS      call Do_t1crs_2RPA(ICORE(I000), MAXCOR, IUHF)
C
      RETURN
      END
