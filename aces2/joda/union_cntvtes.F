#include "flags.h"
      SUBROUTINE UNION_CNTVTES(IBNDTO, IATOMICNMBER, NRATMS)
C
C For LST/QST transition state searches the redundant internal
C coordinates are built from the union of the bonds of the current
C structure and of the reactant and product structures (C. Peng, P.
C Y. Ayala, H. B. Schlegel, M. J. Frisch, J. Comput. Chem. 17, 49
C (1996)), so that the bonds being formed and broken are part of the
C coordinate set. IBNDTO holds the connectivities of the current
C structure on entry; the bonds of the reactant and product
C (RXSTRUCT, PRSTRUCT) are added. Nothing is done unless OPT_CONTROL
C is LST or QST and both structures are present. 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
#include "mxatms.par"
#include "machsp.com"
#include "jodaflags.com"
C
      INTEGER FRAGSCR
      CHARACTER*8 RECNAM(2)
      DIMENSION IBNDTO(NRATMS, NRATMS), IATOMICNMBER(NRATMS),
     &          XYZ(3*MXATMS), BNDLENTHS(MXATMS*MXATMS),
     &          SMOFCOVRADI(MXATMS*MXATMS), NCONPRCNTR(MXATMS),
     &          IBND2(MXATMS*MXATMS), MARK_FRAGMENTS(MXATMS*MXATMS),
     &          LENGTH_FRAGMENTS(MXATMS), FRAGSCR(MXATMS*MXATMS),
     &          BLNGT_IJFRGS(MXATMS)
      DATA RECNAM /'RXSTRUCT', 'PRSTRUCT'/
C
      IF (IFLAGS2(h_IFLAGS2_opt_control) .NE. 1 .AND.
     &    IFLAGS2(h_IFLAGS2_opt_control) .NE. 2) RETURN
      CALL GETREC(0, 'JOBARC', RECNAM(1), LEN1, XYZ)
      CALL GETREC(0, 'JOBARC', RECNAM(2), LEN2, XYZ)
      IF (LEN1 .LE. 0 .OR. LEN2 .LE. 0) RETURN
C
      NADD = 0
      DO ISTR = 1, 2
         CALL GETREC(20, 'JOBARC', RECNAM(ISTR), 3*NRATMS*IINTFP, XYZ)
         CALL IZERO(IBND2, NRATMS*NRATMS)
         CALL ASIGN_CNTVTES(XYZ, BNDLENTHS, IATOMICNMBER,
     &                      SMOFCOVRADI, NCONPRCNTR, IBND2,
     &                      MAXCNTVS, NRATMS, MARK_FRAGMENTS,
     &                      LENGTH_FRAGMENTS, FRAGSCR, BLNGT_IJFRGS)
         DO J = 2, NRATMS
            DO I = 1, J - 1
               IJ = I + (J-1)*NRATMS
               JI = J + (I-1)*NRATMS
               IF ((IBND2(IJ) .GT. 0 .OR. IBND2(JI) .GT. 0) .AND.
     &             IBNDTO(I,J) .LE. 0 .AND. IBNDTO(J,I) .LE. 0) THEN
                  IBNDTO(I,J) = 1
                  IBNDTO(J,I) = 1
                  NADD = NADD + 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
      IF (NADD .GT. 0) WRITE(6,"(a,I3,a)") "  LST/QST: ", NADD,
     &   " bonds of the reactant/product added to the coordinates."
C
      RETURN
      END
