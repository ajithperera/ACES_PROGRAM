













































































































































































































      SUBROUTINE PSI4DBG_SETMET()

      IMPLICIT INTEGER(A-Z)

      LOGICAL NONHF
      LOGICAL CIS,EOM
      LOGICAL PCCD,CCD,LCCD,MBPT2








c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/






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




      COMMON/METH/PCCD,CCD,LCCD,MBPT2
      COMMON/EXCITE/EOM,CIS
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE 
C
      PCCD    = .FALSE.
      CCD     = .FALSE.
      LCCD    = .FALSE.
      EOM     = .FALSE.
      CIS     = .FALSE.
      NONHF   = .FALSE.
      RLE     = .FALSE.
      MBPT2   = .FALSE.

      IF(IFLAGS(38)  .EQ.1) NONHF = .TRUE.
      IF(IFLAGS(87) .EQ.3) EOM   = .TRUE.
      IF(IFLAGS(87) .EQ.5) CIS   = .TRUE.

      IF(IFLAGS(21) .Ne.3) RLE  = .TRUE.

      IF (IFLAGS(2)  .EQ.52) PCCD   = .TRUE.
      IF (IFLAGS(2)  .EQ.53) LCCD   = .TRUE.
      IF (IFLAGS(2)  .EQ.54)  CCD   = .TRUE.
      IF (IFLAGS(2)  .EQ.1) MBPT2   = .TRUE.

C Here LCCD and CCD stands for oo-LCCD and oo-CCD respectively.

      IF (PCCD) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right Pair coupled cluster",
     +                     " doubles (pCCD) equations are solved."
          Write(6,*)
      ELSE IF (LCCD) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right linear coupled cluster",
     +                     " doubles (LCCD) equations are solved."
          Write(6,*)
      ELSE IF (CCD) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right coupled cluster",
     +                     " doubles (CCD) equations are solved."
          Write(6,*) 
      ELSE IF (MBPT2) THEN
          Write(6,*) 
          Write(6,"(a,a)") " The left and right MBPT(2) equations",
     +                     " are solved."
          Write(6,*) 
      ELSE
          Write(6,*) 
          Write(6,"(2a)") " Incorrect choice of method. Choose pCCD",
     +                    " OO-CCD or OO-LCCD instead"
          Call Errex
      ENDIF 

      RETURN
      END
