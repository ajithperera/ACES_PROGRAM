










      SUBROUTINE SETGD
C
C     Set DCORGD, indicating which GD vectors to calculate
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)

c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
C
C Used from common blocks:
C  ABAINF: GDALL
C  CBITRO: COMPAR
C  NUCLEI: DCORD, DCORGD, NTRACO, ITRACO
C
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      LOGICAL         COMPAR, SKIP, NOROT,
     *                HESTRO, GDTRO, RDTRO, TROGRD, TROHES, TRODIP
      COMMON /CBITRO/ THRESH, IPRINT, COMPAR, SKIP, NOROT,
     *                HESTRO, GDTRO, RDTRO, TROGRD, TROHES, TRODIP
      CHARACTER NAMEX*6
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON /NUCLEC/ NAMEX(MXCOOR)
C
C     DCORGD is initialized to all .true. in HESINP subroutine/
C
      IF (GDALL .AND. .NOT. COMPAR) THEN
         DO 100 J = 1,NTRACO
            ICOOR = ITRACO(J)
            INUC  = (ICOOR-1)/3 + 1
            ICOOR = ICOOR - (INUC-1)*3
            DCORGD(INUC,ICOOR) = .FALSE.
  100    CONTINUE
      ELSE
         DO 300 I = 1, MXCENT
            DO 200 J = 1, 3
               DCORGD(I,J) = DCORD(I,J)
  200       CONTINUE
  300    CONTINUE
      END IF
      RETURN
      END
