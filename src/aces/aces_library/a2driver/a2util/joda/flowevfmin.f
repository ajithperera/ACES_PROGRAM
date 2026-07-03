










      SUBROUTINE FLOWEVFMIN(SCRATCH, GRDHES, HESMOD, DIAGHES, LMBDAN,
     &                      MORSE, NX, NOPT)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DOUBLE PRECISION LMBDAN
      LOGICAL MORSE
C

































































































































































































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




c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end
C
      DIMENSION SCRATCH(NX*NX), GRDHES(NOPT), HESMOD(NOPT, NOPT),
     &          DIAGHES(NOPT, NOPT)
C
      DO 10 I = 1, NOPT
         DO 20 J = 1, NOPT

            IF (LAMBDN .EQ. 0.0D0 .AND. 
     &          DABS(HESMOD(I,I)) .LT. 1.0D0-8)  THEN
C
                 SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
C
            ELSE 
C
               IF (DABS((HESMOD(I,I) - LMBDAN)) .LT. 1.0D0-8) THEN

                  SCRATCH(J+NOPT) = SCRATCH(J+NOPT) - 0.0D0
C
               ELSE
                  SCRATCH(J+NOPT) = SCRATCH(J+NOPT)-
     &                             GRDHES(I)*DIAGHES(J,I)/
     &                             (HESMOD(I,I) - LMBDAN)
C
              ENDIF
            ENDIF
C
 20      CONTINUE
 10   CONTINUE 
        Write(6,*) "The unscaled RFO step"
        Write(6, "(3F10.5)") (SCRATCH(J+NOPT), J=1, NOPT)
C
C Notice that there is no Morse adjustments can be made for
C pure Cartesian OLPtimization (no connections), Also, no
C Morse scaling can be done for redundent internals.
C Ajith Perera,08/2008, 2011.

      IF (MORSE) THEN
         IF (iFlags2(5) .ge. 3) THEN
CSSS           CALL DOMORSEXYZ_RIC(SCRATCH, NOPT, NX, LUOUT)
         ELSE IF (iFlags2(5) .eq. 1) THEN
           CALL DOMORSEZMT(SCRATCH, NOPT, NX, LUOUT)
         ENDIF
      ENDIF
C
      RETURN
      END
