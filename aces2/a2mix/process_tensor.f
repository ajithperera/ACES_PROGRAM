










      SUBROUTINE PROCESS_TENSOR(CC_ARRAY, CC_TENSOR, GRID_POINTS,
     &                          NTPERT, IPOINT, IJPAIR, MAXCENT,
     &                          MAXGRD, MAXPAIRS, NATOMS, 
     &                          NPOINT_PAIR)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION CC_ARRAY(MAXCENT,MAXGRD*MAXPAIRS), 
     &          NPOINT_PAIR(MAXPAIRS), GRID_POINTS(3,MAXPAIRS*MAXGRD),
     &          CC_TENSOR(MAXPAIRS,MAXGRD*MAXPAIRS)
C
      CHARACTER*12 FILENAME
      LOGICAL FILE_EXIST

C
C Build the coupling tensor from the coupling array.
C
      IJPAIR = 0
      IKEEP  = 0
      DO IATOM = 1, NATOMS
         DO JATOM = 1, NATOMS
C
            IF (IATOM .LT. JATOM)  THEN
                IJPAIR = IJPAIR + 1

                IMV   = (IJPAIR - 1)*NPOINT_PAIR(IJPAIR)
                IFX   =  IJPAIR*NPOINT_PAIR(IJPAIR)
C
                DO NPOINT = 1, NPOINT_PAIR(IJPAIR)
                     IMV = IMV + 1
                   IKEEP = IKEEP + 1
                   CC_TENSOR(IJPAIR, IKEEP) = CC_ARRAY(JATOM,IMV)
     &                                      + CC_ARRAY(IATOM,IFX)
                ENDDO
C
            ENDIF
         ENDDO
      ENDDO
C
C
      RETURN
      END
