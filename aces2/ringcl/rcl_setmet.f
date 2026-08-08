













































































































































































































      SUBROUTINE RCL_SETMET(HBAR_ONLY)

      IMPLICIT INTEGER(A-Z)

      LOGICAL NONHF
      LOGICAL CIS,EOM,RPA,EOM_SDRCCD,EOM_SRCCD,EOM_SFDRCCD,EOM_SFRCCD
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD,RLE,EOM_SDXRCCD,EOM_SFDXRCCD

      LOGICAL TRIPIT,TRIPNI,TRIPNI1,EOM_SRXCCD,EOM_SFRXCCD,
     &        HBAR_ONLY
     







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




      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/EXCITE/EOM,CIS,RPA,EOM_SRCCD,EOM_SDRCCD,EOM_SFRCCD,
     &              EOM_SFDRCCD,EOM_SDXRCCD,EOM_SFDXRCCD
      COMMON/NHFREF/NONHF
      COMMON/EXTRAPO/RLE 
C
      MBPT2       = .FALSE.
      MBPT3       = .FALSE.
      M4DQ        = .FALSE.
      M4SDQ       = .FALSE.
      M4SDTQ      = .FALSE.
      CCD         = .FALSE.
      CC2         = .FALSE.
      QCISD       = .FALSE.
      CCSD        = .FALSE.
      UCC         = .FALSE.
      NONHF       = .FALSE.
      TRIPIT      = .FALSE.
      TRIPNI      = .FALSE.
      TRIPNI1     = .FALSE.
      RCCD        = .FALSE. 
      DRCCD       = .FALSE.
      EOM         = .FALSE.
      CIS         = .FALSE.
      NONHF       = .FALSE.
      RLE         = .FALSE.
      RPA         = .FALSE. 
      EOM_SFDRCCD = .FALSE.
      EOM_SFRCCD  = .FALSE.
      EOM_SDRCCD  = .FALSE.
      EOM_SRCCD   = .FALSE.
      EOM_SRXCCD  = .FALSE.
      EOM_SFRXCCD = .FALSE.

      IF(IFLAGS(38)  .EQ.1)  NONHF        = .TRUE.
      IF(IFLAGS(87) .EQ.3)  EOM          = .TRUE.
      IF(IFLAGS(87) .EQ.5 .OR. 
     +   IFLAGS(87) .EQ.1)  CIS          = .TRUE.
      IF(IFLAGS(87) .EQ.2)  RPA          = .TRUE.
      IF(IFLAGS(87) .EQ.17) EOM_SFRCCD   = .TRUE.
      IF(IFLAGS(87) .EQ.18) EOM_SFDRCCD  = .TRUE.
      IF(IFLAGS(87) .EQ.13) EOM_SDRCCD   = .TRUE.
      IF(IFLAGS(87) .EQ.15) EOM_SRCCD    = .TRUE.
      IF(IFLAGS(87) .EQ.19) EOM_SDXRCCD  = .TRUE.
      IF(IFLAGS(87) .EQ.20) EOM_SFDXRCCD = .TRUE.

      IF(IFLAGS(21) .Ne.3) RLE  = .TRUE.

      IF (IFLAGS(2)  .EQ.48) RCCD   = .TRUE.
      IF (IFLAGS(2)  .EQ.49) DRCCD  = .TRUE.

      IF (RPA .AND. RCCD)  EOM_SFRCCD  = .TRUE.
      IF (RPA .AND. DRCCD) EOM_SFDRCCD = .TRUE.

      IF (HBAR_ONLY) RETURN
 
      IF (RCCD) THEN
          Write(6,"(a,a)") " Ring coupled cluster doubles (rCCD)",
     +                     " lambda equations are solved."
          Write(6,*)
      ELSE IF (DRCCD) THEN
          Write(6,"(a,a)") " Direct Ring coupled cluster doubles",
     +                     " (drCCSD) lambda equations are solved."
          Write(6,*)
      ENDIF 
      RETURN
      END
