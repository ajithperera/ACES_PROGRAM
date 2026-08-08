











C CREATES LISTS USED IN EOM CALCULATIONS











C
      SUBROUTINE MAKE_LISTS4_PERTRB_DENS(ICORE, MAXCOR, IUHF)
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL EOMCC, NODAVID, CIS, ESPROP
C
      DIMENSION ICORE(MAXCOR),MAXOO(2),MAXVV(2),MAXVO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
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
      CALL ACES_IO_REMOVE(54,'DERGAM')
      INEWFIL= 0
      IZILCH = 0
      EOMCC  = .FALSE.
      NODAVID= .FALSE.
          CIS= .FALSE.
       ESPROP= .FALSE.        
C
C CALCULATE MAXIMUM SIZES OF VV, OO AND VO VECTORS
C
      MAXVV(1)=0
      MAXOO(1)=0
      MAXVO(1)=0

      DO IRREP = 1, NIRREP
         MAXVV(1)=MAX(MAXVV(1),IRPDPD(IRREP,19))
         MAXOO(1)=MAX(MAXOO(1),IRPDPD(IRREP,21))
         MAXVO(1)=MAX(MAXVO(1),IRPDPD(IRREP,9))
      END DO
C
      IF (IUHF .EQ. 1) THEN
         MAXVV(2)=0
         MAXOO(2)=0
         MAXVO(2)=0
         DO IRREP = 1, NIRREP
            MAXVV(2)=MAX(MAXVV(2),IRPDPD(IRREP,20))
            MAXOO(2)=MAX(MAXOO(2),IRPDPD(IRREP,22))
            MAXVO(2)=MAX(MAXVO(2),IRPDPD(IRREP,10))
         END DO
      ENDIF
C
C CREATE SINGLES VECTOR LISTS AND DIPOLE MOMENT LISTS
C
      DO ISPIN=1,1+IUHF
C
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,INEWFIL,0)
        INEWFIL=0
CSSS        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,493,INEWFIL,0)
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,INEWFIL,0)
CSSS        CALL UPDMOI(1,MAXVO(ISPIN),9,447+ISPIN,INEWFIL,0)
        call aces_list_memset(ispin,  490,0)
CSSS        call aces_list_memset(ispin,  493,0)
        call aces_list_memset(ispin+2,490,0)
CSSS        call aces_list_memset(9,      447+ISPIN,0)
C
C
      ENDDO
C
C
C CREATE AREA FOR Q(AB) AND Q(IJ) THREE-BODY INTERMEDIATES
C
C
C CALCULATE MAXIMUM SIZE OF DAVIDSON LISTS
C
C
C NOW MAKE DENOMINATOR AND T2 LISTS   
C       
      DO ISPIN=3,3-2*IUHF,-1
C
C DENOMINATOR AND T2 LISTS
C
         TTYPEL=ISYTYP(1,43+ISPIN)
         TTYPER=ISYTYP(2,43+ISPIN)
         CALL INIPCK2(1,TTYPEL,TTYPER,460+ISPIN,IZILCH,IZILCH,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,447+ISPIN,IZILCH,IZILCH,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,443+ISPIN,IZILCH,IZILCH,1)
C
C FOR AO-BASED ALGORITHMS, NEED AO T2 LISTS
C
         IF (IFLAGS(93).EQ.2) THEN
            CALL INIPCK2(1,TTYPEL,TTYPER,280+ISPIN,IZILCH,IZILCH,1)
            TTYPEL=15
            CALL INIPCK2(1,TTYPEL,TTYPER,213+ISPIN,IZILCH,IZILCH,1)
            CALL INIPCK2(1,TTYPEL,TTYPER,463+ISPIN,IZILCH,IZILCH,1)
         ENDIF
C
      ENDDO
C
C MAKE RESORTED R2 AND L2 LISTS
C

      CALL INIPCK2(1,9,10,426,IZILCH,0,1)
      CALL INIPCK2(1,11,11,428,IZILCH,0,1)
      CALL INIPCK2(1,9,9,434,IZILCH,0,1)
      CALL INIPCK2(1,9,10,437,IZILCH,0,1)
      CALL INIPCK2(1,11,12,439,IZILCH,0,1)
      CALL INIPCK2(1,9,9,440,IZILCH,0,1)
      CALL INIPCK2(1,10,10,441,IZILCH,0,1)
      CALL INIPCK2(1,9,10,442,IZILCH,0,1)
      CALL INIPCK2(1,11,12,443,IZILCH,0,1)
      CALL INIPCK2(1,9,9,454,IZILCH,0,1)
      CALL INIPCK2(1,9,10,457,IZILCH,0,1)
      CALL INIPCK2(1,11,12,459,IZILCH,0,1)

      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,426,0)
         call aces_list_memset(iGrp,428,0)
         call aces_list_memset(iGrp,434,0)
         call aces_list_memset(iGrp,437,0)
         call aces_list_memset(iGrp,439,0)
         call aces_list_memset(iGrp,440,0)
         call aces_list_memset(iGrp,441,0)
         call aces_list_memset(iGrp,442,0)
         call aces_list_memset(iGrp,443,0)
         call aces_list_memset(iGrp,454,0)
         call aces_list_memset(iGrp,457,0)
         call aces_list_memset(iGrp,459,0)
      enddo
C
      IF (IUHF.NE.0) THEN
         CALL INIPCK2(1,9,9,424,IZILCH,0,1)
         CALL INIPCK2(1,10,10,425,IZILCH,0,1)
         CALL INIPCK2(1,10,9,427,IZILCH,0,1)
         CALL INIPCK2(1,12,12,429,IZILCH,0,1)
         CALL INIPCK2(1,10,10,435,IZILCH,0,1)
         CALL INIPCK2(1,10,9,436,IZILCH,0,1)
         CALL INIPCK2(1,12,11,438,IZILCH,0,1)
         CALL INIPCK2(1,10,10,455,IZILCH,0,1)
         CALL INIPCK2(1,10,9,456,IZILCH,0,1)
         CALL INIPCK2(1,12,11,458,IZILCH,0,1)

      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,424,0)
         call aces_list_memset(iGrp,425,0)
         call aces_list_memset(iGrp,427,0)
         call aces_list_memset(iGrp,429,0)
         call aces_list_memset(iGrp,435,0)
         call aces_list_memset(iGrp,436,0)
         call aces_list_memset(iGrp,438,0)
         call aces_list_memset(iGrp,455,0)
         call aces_list_memset(iGrp,456,0)
         call aces_list_memset(iGrp,458,0)
      end do
      END IF
C
      IF (IFLAGS(91).GT.1) THEN
         CALL INIPCK2(1,1,3,114,IZILCH,0,1)
         CALL INIPCK2(1,2,4,115,IZILCH,0,1)
         CALL INIPCK2(1,13,  14,  116,IZILCH,0,1)
      END IF
C
      RETURN
      END
