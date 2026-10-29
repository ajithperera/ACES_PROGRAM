










      SUBROUTINE BLKCHK(ICORE, MAXCOR, IUHF, IRREPX, IPERT, NPERTX,
     &                  NTPERT)
C
C DEBUG: BLOCK-BY-BLOCK CHECK OF THE RESPONSE DENSITY (UHF, C1 ONLY).
C BOTH VCCEH'S X(I,J) = <M(L)_I|T(1)_J> (LISTS 374/373) AND PDENS'
C TR(D(1)_J H_I) ARE LINEAR IN T(1)_J, SO EACH CAN BE SPLIT OVER THE
C BLOCKS OF T(1)_J: T1(A), T1(B), T2(AB), T2(BB), T2(AA) (LOADVEC1
C ORDER). FOR EACH BLOCK THE DENSITY IS REBUILT WITH ONLY THAT BLOCK
C OF T(1)_J NONZERO AND CONTRACTED WITH THE ATOMIC FC INTEGRALS.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPN_SPN, POLAR
      CHARACTER*4 BLKNAM(6), DNAM(5)
      COMMON /DBLKMASK/ IDMASK
      DATA DNAM /'all ','OO  ','VV  ','VO  ','OV  '/
C


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
      DIMENSION ICORE(MAXCOR), NOFF(6), DOTS(5), SCRT(3600), POLT(9)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA BLKNAM /'T1A ','T1B ','T2AB','T2BB','T2AA','ALL '/
C
      IF (IUHF .EQ. 0 .OR. NTPERT .GT. 60 .OR. NTPERT .LT. 1) RETURN
      POLAR   = .FALSE.
      SPN_SPN = .FALSE.
C
      N1A = IRPDPD(IRREPX, 9)
      N1B = IRPDPD(IRREPX,10)
      NAB = IDSYMSZ(IRREPX, ISYTYP(1,46), ISYTYP(2,46))
      NBB = IDSYMSZ(IRREPX, ISYTYP(1,45), ISYTYP(2,45))
      NAA = IDSYMSZ(IRREPX, ISYTYP(1,44), ISYTYP(2,44))
      NOFF(1) = 0
      NOFF(2) = NOFF(1) + N1A
      NOFF(3) = NOFF(2) + N1B
      NOFF(4) = NOFF(3) + NAB
      NOFF(5) = NOFF(4) + NBB
      NOFF(6) = NOFF(5) + NAA
      NSIZE   = NOFF(6)
C
      IV = 1
      IL = IV + NSIZE*IINTFP
      IF (IL + NSIZE*IINTFP .GT. MAXCOR) RETURN
C
C VCCEH SIDE
C
      CALL GETLST(ICORE(IV), IPERT, 1, 1, IRREPX, 373)
      DO I = 1, NPERTX
         CALL GETLST(ICORE(IL), I, 1, 1, IRREPX, 374)
         CALL BLKDOT(ICORE(IL), ICORE(IV), NOFF, DOTS)
         WRITE(6,"(a,2I3,5F16.10)") " BLKCHK vcceh I,J ", I, IPERT,
     &                              (DOTS(K), K = 1, 5)
      ENDDO
C
C PDENS SIDE, ONE BLOCK OF T(1)_J AT A TIME
C
      DO K = 1, 6
         CALL GETLST(ICORE(IV), IPERT, 1, 1, IRREPX, 373)
         IF (K .LE. 5) CALL BLKMASK(ICORE(IV), NOFF, K)
         CALL PUTLST(ICORE(IV), 1, 1, 1, 1, 490)
         CALL PUTLST(ICORE(IV+NOFF(2)*IINTFP), 1, 1, 1, 2, 490)
         CALL PUTALL(ICORE(IV+NOFF(3)*IINTFP), NAB, IRREPX, 446)
         CALL PUTALL(ICORE(IV+NOFF(4)*IINTFP), NBB, IRREPX, 445)
         CALL PUTALL(ICORE(IV+NOFF(5)*IINTFP), NAA, IRREPX, 444)
C
         WRITE(6,"(a,a4)") " BLKCHK density from block ", BLKNAM(K)
         CALL GEN_PERTRB_RSPNS_DENS(ICORE, MAXCOR, IRREPX, IUHF)
C
C DENSITY-BLOCK SPLIT (OO, VV, VO, OV) OF THIS T(1) BLOCK'S CONTRIBUTION
C
         DO IDM = 2, 5
            IDMASK = IDM - 1
            CALL PRP_DEN4PROPS(ICORE, MAXCOR, NAO, NATOMS, NCENTER, 99,
     &                         IRREPX, IUHF)
            CALL ZERO(SCRT, NTPERT*NTPERT)
            CALL ZERO(POLT, 9)
            ILEFT = 1 + NAO*NAO*IINTFP
            CALL COMP_2NDORD_PROPS(ICORE(1), ICORE(ILEFT), NTPERT,
     &                             MAXCOR-ILEFT, NAO, NATOMS, 1, IUHF,
     &                             SPN_SPN, POLAR, SCRT, POLT)
            IF (SPN_SPN) WRITE(6,"(a,a4,a,a4,9F14.8)") " BLKCHK dblk ",
     &         BLKNAM(K), " ", DNAM(IDM), (SCRT(I), I = 1, NTPERT)
         ENDDO
         IDMASK = 0
         CALL PRP_DEN4PROPS(ICORE, MAXCOR, NAO, NATOMS, NCENTER, 99,
     &                      IRREPX, IUHF)
C
         CALL ZERO(SCRT, NTPERT*NTPERT)
         CALL ZERO(POLT, 9)
         ILEFT = 1 + NAO*NAO*IINTFP
         CALL COMP_2NDORD_PROPS(ICORE(1), ICORE(ILEFT), NTPERT,
     &                          MAXCOR-ILEFT, NAO, NATOMS, 1, IUHF,
     &                          SPN_SPN, POLAR, SCRT, POLT)
         IF (POLAR) THEN
            WRITE(6,"(a,a4,I3,3F18.10)") " BLKCHK pdens ", BLKNAM(K),
     &           IPERT, (POLT(I), I = 1, 3)
         ENDIF
         IF (SPN_SPN) THEN
            WRITE(6,"(a,a4,I3,9F14.8)") " BLKCHK pdens ", BLKNAM(K),
     &           IPERT, (SCRT(I), I = 1, NTPERT)
         ENDIF
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE BLKDOT(A, B, NOFF, DOTS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(*), B(*), NOFF(6), DOTS(5)
      DO K = 1, 5
         N = NOFF(K+1) - NOFF(K)
         DOTS(K) = 0.0D0
         IF (N .GT. 0) DOTS(K) = SDOT(N, A(NOFF(K)+1), 1,
     &                                B(NOFF(K)+1), 1)
      ENDDO
      RETURN
      END
C
      SUBROUTINE BLKMASK(V, NOFF, KEEP)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION V(*), NOFF(6)
      DO K = 1, 5
         IF (K .NE. KEEP) THEN
            DO I = NOFF(K)+1, NOFF(K+1)
               V(I) = 0.0D0
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END
