










      SUBROUTINE GEN_RICS(CARTCOORD,TOTREDNCO, TOTNOFBND, TOTNOFANG,
     &                    TOTNOFDIH, IREDUNCO, REDUNCO)
C
C Generate the redundent internal coordinates given the bond, angle and
C dihedral angle assignments.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
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
C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end


C 
      PARAMETER(THRESHOLD = 5.0D00, EPSILON = 1.0D-10)
C
      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH,
     &        FRAGSCR
      LOGICAL TEST_4CHANGE
C
      DIMENSION CARTCOORD(3*MXATMS), BNDLENTHS(MXATMS, MXATMS), 
     &          IREDUNCO(4, MAXREDUNCO), 
     &          REDUNCO(MAXREDUNCO),        
     &          VECBA(3), VECBC(3), 
     &          VECBAD(3), VECCB(3), VECAB(3), VECCD(3), VECABC(3),
     &          VECBCD(3)
C      
      DATA  IZERO, MONE, ZERO /0, -1, 0.0D0/
C
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0 
      PI     = (ATAN(DFLOAT(1))*DFLOAT(4))
C
      NW_TOTREDNCO = IZERO
      BHORTOANG    = 0.529177249D0
C
      Write(6,"(a)"), "The redundent internal coordinate assignments"
      Write(6,*)
      Do i = 1, TOTREDNCO
         Write(6,111) (iredunco(j, i), j=1, 4)
      Enddo
 111  Format(5X, 4(I3, 1X))
C
C Assign bond coordiante
C
      Write(6,*)
  
      DO IBNDS = 1, TOTNOFBND
         IF (IREDUNCO(2, IBNDS) .NE. 0) THEN
            ICON1  = IREDUNCO(1, IBNDS)
            ICON2  = IREDUNCO(2, IBNDS)
            DISTAB = DIST(CARTCOORD(3*ICON1 - 2), 
     &               CARTCOORD(3*ICON2 - 2))
            REDUNCO(IBNDS) = DISTAB
C     
            WRITE(6,"(a,F10.5)") "The bond distance =",DISTAB*0.529177249d0
         ENDIF
      ENDDO
C
C Assign bond angle Coordinates.
C
      Write(6,*)
      DO IANGS = (TOTNOFBND + 1), (TOTNOFANG + TOTNOFBND)
     
         IF (IREDUNCO(4, IANGS) .EQ. MONE) THEN
C
C This is for linear arrangements (two angle coordinates are
C required). 
C
            REDUNCO(IANGS) = 180.0D0*DINVPI
            WRITE(6,"(a,F10.5)") "The Bond Angle =", ANGL/DINVPI
C
         ELSE
C
C This is for the non-linear arrangements.
C
            ICON1  = IREDUNCO(1, IANGS)
            ICON2  = IREDUNCO(2, IANGS)  
            ICON3  = IREDUNCO(3, IANGS)
            CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2),
     &               VECBC, 1)
            CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON1 - 2),
     &               VECBA, 1)
            ANGL    = ANGLE(VECBC, VECBA, 3)*DINVPI
            REDUNCO(IANGS) = ANGL
            WRITE(6,"(a,F10.5)") "The Bond Angle =", ANGL/DINVPI
         ENDIF
C
      ENDDO
C
C Assign bond dihedral angle Coordinates.
C
      Write(6,*)
      DO IDIHS = (TOTNOFANG + TOTNOFBND + 1),  TOTREDNCO
         ICON1  = IREDUNCO(1, IDIHS)
         ICON2  = IREDUNCO(2, IDIHS)
         ICON3  = IREDUNCO(3, IDIHS)
         ICON4  = IREDUNCO(4, IDIHS)
         CALL VEC(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2),
     &         VECAB, 1)
         CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2),
     &         VECBC, 1)
         CALL VEC(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2),
     &         VECCB, 1)
         CALL VEC(CARTCOORD(3*ICON4 - 2), CARTCOORD(3*ICON3 - 2),
     &         VECCD, 1)
C
         DISTAB = DIST(CARTCOORD(3*ICON1 - 2), 
     &            CARTCOORD(3*ICON2 - 2))
         DISTBC = DIST(CARTCOORD(3*ICON3 - 2), 
     &            CARTCOORD(3*ICON2 - 2))
         DISTCD = DIST(CARTCOORD(3*ICON3 - 2), 
     &            CARTCOORD(3*ICON4 - 2))
C
C First evaluate the dihedral angle. This is calculated as the
C angle between the two unit vectors that are perpendicular to
C the ABC and BCD planes for A-B-C-D pattern.
C
         CALL CROSS(VECAB, VECBC, VECABC, 1)
C
C We need the vector CD, obtain that from vec DC (note the misnomar)
C
         CALL DSCAL(3, -1.0D0, VECCD, 1)
         CALL CROSS(VECBC, VECCD, VECBCD, 1)
C
         DANG = (ANGLE(VECABC, VECBCD, 3))*DINVPI
C
C--- Obtainig the sign of the dihedral angle:
C
         CALL DSCAL(3, -1.0D0, VECCD, 1)
         CALL CROSS(VECCD, VECCB, VECBCD, 0)
C
         SENSE = DDOT(3, VECAB, 1, VECBCD, 1)

         IF (SENSE .GT. 0.0D0) DANG = -DANG
         REDUNCO(IDIHS) = DANG
         
            WRITE(6,"(a,F10.5)") "The Dihedral Angle =", DANG/DINVPI
C
C
      ENDDO 
C
      RETURN
      END

