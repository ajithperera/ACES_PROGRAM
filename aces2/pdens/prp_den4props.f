










C
      SUBROUTINE PRP_DEN4PROPS(ICORE, MAXCOR, NAO, NATOM, 
     &                         NCENTR, IPERT, IRREPX, IUHF)
C



































































































































































































      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C
      CHARACTER*8 RECNAMEA, RECNAMEB
      LOGICAL TRIPLET_PERT, ANTI
C
      DIMENSION LENVV(2), LENOO(2), LENVO(2), INDEX(MXATMS)
      DIMENSION DMOM(3), DSO(3), ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /FILES/ LUOUT, MOINTS
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
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




C
      IONE =  1
      ONE  =  1.00D+00
      ONEM = -1.00D+00
      TWO  =  2.00D+0
      ANTI = .FALSE.
      If (Iflags(18) .EQ. 8) ANTI = .TRUE.
    
C
      CALL IZERO (LENOO, 2)
      CALL IZERO (LENVV, 2)
      CALL IZERO (LENVO, 2)
      NMO = NOCCO(1) + NVRTO(1)
      FACT = 2.00D+00 - DFLOAT(IUHF)
C
      CALL GETREC (20, 'JOBARC', 'NBASTOT', IONE, NAO)
      CALL GETREC (20, 'JOBARC', 'NATOMS  ', IONE, NATOM)
      CALL GETREC (20, 'JOBARC', 'NREALATM', IONE, NCENTR)
C
      DO 10 ISPIN = 1, (IUHF + 1)
         LENOO(ISPIN) = IRPDPD(IRREPX, 20+ISPIN)
         LENVV(ISPIN) = IRPDPD(IRREPX, 18+ISPIN)
         LENVO(ISPIN) = IRPDPD(IRREPX,  8+ISPIN)
 10   CONTINUE
C
C Allocate memory for ALPHA density
C
      I000 = IONE
      I010 = I000 + LENOO(1)*IINTFP
      I020 = I010 + LENVV(1)*IINTFP
      I030 = I020 + LENVO(1)*IINTFP
      I040 = I030 + LENVO(1)*IINTFP
      I050 = I040 + NAO*NAO*IINTFP
      I060 = I050 + NAO*NAO*IINTFP
      I070 = I060 + NAO*NAO*IINTFP
      I080 = I070 + NAO*NAO*IINTFP
C
C Allocate memory for BETA density (UHF only)
C
      IF (IUHF .NE. 0) THEN 
         I120 = I080
         I130 = I120 + LENOO(2)*IINTFP
         I140 = I130 + LENVV(2)*IINTFP
         I150 = I140 + LENVO(2)*IINTFP
         I160 = I150 + LENVO(2)*IINTFP
         I170 = I160 + NAO*NAO*IINTFP
         I180 = I170 + NAO*NAO*IINTFP
         I190 = I180 + NAO*NAO*IINTFP
         I200 = I190 + NAO*NAO*IINTFP
      ENDIF
c
C      DO  ICOUNT = 1, MXATM
C         INDEX(ICOUNT) = ICOUNT
C      ENDDO
C    
C Read the OO, VV, OV blocks of the perturbed density GAMLAM file
C
      CALL READ_PDENS(ICORE(I000), ICORE(I010), ICORE(I020),
     &               ICORE(I030), LENOO(1), LENVV(1), LENVO(1), 1)
CSSS      CALL DCOPY(LENVO(1)*IINTFP,ICORE(I020),1,ICORE(I030),1) 
      CALL EXPDEN (ICORE(I000), ICORE(I010), ICORE(I020),
     &             ICORE(I030), ICORE(I040), NMO, 1, 1, .FALSE.)
      IF (ANTI) CALL SCALEBLOCK(ICORE(I040),NOCCO(1),NOCCO(2),
     &                          NVRTO(1),NVRTO(2), NMO,.TRUE.,
     &                          .FALSE.,.FALSE.,.FALSE.,1,-1.0D0)
CSSS      CALL SYMMET2 (ICORE(I040), NMO)
C-----------ZERO-BLOCK TEST
CSSS      CALL ZEROBLOCK(ICORE(I040),NOCCO(1),NOCCO(2),NVRTO(1),NVRTO(2),
CSSS     &               NMO,.FALSE.,.FALSE.,.FALSE.,.FALSE.,1)
C--------------------
      CALL MO2AO3 (ICORE(I040), ICORE(I050), ICORE(I060), ICORE(I070),
     &             NAO, NMO, 1)
      WRITE(RECNAMEA, "(A,I2)") "PTDENA", IPERT
      CALL PUTREC(20, "JOBARC", RECNAMEA, NMO*NMO*IINTFP,
     &            ICORE(I050))
C
C Put ALPHA + BETA density into ICORE(I050). 
C
      If (IUHF .NE. 0) THEN
         CALL READ_PDENS(ICORE(I120), ICORE(I130), ICORE(I140),
     &                  ICORE(I150), LENOO(2), LENVV(2), LENVO(2), 2)
CSSS         CALL DCOPY(LENVO(2)*IINTFP,ICORE(I140),1,ICORE(I150),1) 
         cALL EXPDEN (ICORE(I120), ICORE(I130), ICORE(I140),
     &                ICORE(I150), ICORE(I160), NMO, 1, 2,
     &               .FALSE.)
      IF (ANTI) CALL SCALEBLOCK(ICORE(I040),NOCCO(1),NOCCO(2),
     &                          NVRTO(1),NVRTO(2), NMO,.TRUE.,
     &                          .FALSE.,.FALSE.,.FALSE.,2,-1.0D0)
C-----------ZERO-BLOCK TEST
CSSS      CALL ZEROBLOCK(ICORE(I160),NOCCO(1),NOCCO(2),NVRTO(1),NVRTO(2),
CSSS     &               NMO,.FALSE.,.FALSE.,.FALSE.,.FALSE.,2)
C------------------
CSSS         CALL SYMMET2 (ICORE(I160), NMO)
         CALL MO2AO3 (ICORE(I160), ICORE(I170), ICORE(I180), 
     &                ICORE(I190), NAO, NMO, 2)
         WRITE(RECNAMEB, "(A,I2)") "PTDENB", IPERT
         CALL PUTREC(20, "JOBARC", RECNAMEB, NMO*NMO*IINTFP,
     &               ICORE(I170))
C
C Build the total density
C
CSSS         CALL SAXPY (NAO*NAO, ONE, ICORE(I170), 1, ICORE(I050), 1)
      ENDIF
C
      TRIPLET_PERT = .FALSE.
      IF (IFLAGS(18) .EQ.  9 .OR. IFLAGS(18) .EQ. 10) 
     &    TRIPLET_PERT = .TRUE.
C
      IF (IUHF .NE. 0) THEN
          IF (TRIPLET_PERT) THEN
             CALL SAXPY (NAO*NAO, ONEM, ICORE(I170), 1, ICORE(I050), 1)
          ELSE
             CALL SAXPY (NAO*NAO, ONE,  ICORE(I170), 1, ICORE(I050), 1)
          END IF
      ELSE
             CALL DSCAL(NAO*NAO, TWO, ICORE(I050), 1)
      ENDIF
C
C Copy the total density to begining to the heap.
C      
      CALL DCOPY(NAO*NAO, ICORE(I050), 1, ICORE(I000), 1)
C
      RETURN
      END
