      SUBROUTINE DENS11(C,DOO,DVV,IRREPC,ISPIN,IUHF)
C
C CALCULATES CONTRIBUTIONS TO THE TDA ONE-PARTICLE DENSITY MATRICES 
C
C                       +        
C    D(ij) = - SUM C(ai) C(aj)
C               a
C
C                           +
C    D(ab) = SUM C(am) C(bm)
C             m
C
C THE OCCUPIED-VIRTUAL DENSITY MATRIX ELEMENTS DO NOT INVOLVE 
C CONTRACTIONS AND ARE NOT CALCULATED HERE
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION C(*),DOO(*),DVV(*)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
      IRREPX=1
      LENGTHC=IRPDPD(IRREPC,9)
C
      IF(IUHF.EQ.0)THEN
       CALL SSCAL(LENGTHC,1.0D0/DSQRT(2.0D0),C,1)
      ENDIF
C 
C FORM OCCUPIED-OCCUPIED CONTRIBUTION
C
      DO 10 IRREPJ=1,NIRREP
       IRREPI=DIRPRD(IRREPJ,IRREPX)
       NUMJ=POP(IRREPJ,ISPIN)
       NUMI=POP(IRREPI,ISPIN)
       IRREPA=DIRPRD(IRREPI,IRREPC)
       NUMA=VRT(IRREPA,ISPIN)
       IOFFCL=ISYMOFF(IRREPI,IRREPC,8+ISPIN)
       IOFFCR=ISYMOFF(IRREPJ,IRREPC,8+ISPIN)
       IOFFD =ISYMOFF(IRREPJ,IRREPX,20+ISPIN)
C The contraction dimension (NUMA, the number of virtuals in IRREPA) can
C legitimately be 0 for a low-symmetry irrep with no virtual orbitals --
C XGEMM's own quick-return guard (tools/xgemm.F) only covers K=0 when
C BETA=1 (meant for the "leave C unchanged" case); this call uses
C BETA=ZILCH=0, so an empty contraction should ZERO this block of DOO
C instead, which XGEMM's guard doesn't do. Found via case 057 (EE_SYM=
C 3-2-2-2, small BASIS=3-21G): NUMA=0 with NUMI=NUMJ=1 reached MKL's own
C DGEMM with an invalid LDA (parameter 8), which MKL's stricter checking
C caught as a hard error where an older/looser BLAS apparently did not.
C Same pre-existing gap in vcceh's own copy of this routine, unfixed
C there -- not currently reachable via any passing pyaces case.
       IF (NUMA.GT.0) THEN
        CALL XGEMM('T','N',NUMI,NUMJ,NUMA,ONEM,C(IOFFCL),NUMA,
     &              C(IOFFCR),NUMA,ZILCH,DOO(IOFFD),NUMI)
       ELSE IF (NUMI*NUMJ.GT.0) THEN
        CALL ZERO(DOO(IOFFD),NUMI*NUMJ)
       ENDIF
10    CONTINUE
      CALL PUTLST(DOO,1,1,1,ISPIN,160)

#ifdef _DEBUG_LVLM
      call output(doo, 1, 5, 1, 5, 5, 5, 1)
      sum = ddot(numj*numi, doo, 1, doo, 1)
      Write(6,*) "Sum =", 4*Sum
#endif
C
C FORM VIRTUAL-VIRTUAL CONTRIBUTION
C
      DO 20 IRREPB=1,NIRREP
       IRREPA=DIRPRD(IRREPB,IRREPX)
       NUMB=VRT(IRREPB,ISPIN)
       NUMA=VRT(IRREPA,ISPIN)
       IRREPM=DIRPRD(IRREPB,IRREPC)
       NUMM=POP(IRREPM,ISPIN)
       IOFFCL=ISYMOFF(IRREPM,IRREPC,8+ISPIN)
       IOFFCR=ISYMOFF(IRREPM,IRREPC,8+ISPIN)
       IOFFD =ISYMOFF(IRREPB,IRREPX,18+ISPIN)
C Same gap as the OO block above -- NUMM (occupied count for IRREPM, the
C contraction dimension here) can legitimately be 0.
       IF (NUMM.GT.0) THEN
        CALL XGEMM('N','T',NUMA,NUMB,NUMM,ONE,C(IOFFCR),NUMA,
     &             C(IOFFCL),NUMB,ZILCH,DVV(IOFFD),NUMA)
       ELSE IF (NUMA*NUMB.GT.0) THEN
        CALL ZERO(DVV(IOFFD),NUMA*NUMB)
       ENDIF
20    CONTINUE
      CALL PUTLST(DVV,1,1,1,2+ISPIN,160)
#ifdef _DEBUG_LVLM
      call output(dvv, 1, 9, 1, 9, 9, 9, 1)
      sum = ddot(numa*numb, dvv, 1, dvv, 1)
      Write(6,*) "Sum =", 4*Sum
#endif
C
      CALL PUTLST(C,1,1,1,ISPIN,490)
      CALL PUTLST(C,1,1,1,2+ISPIN,490)
C
      RETURN
      END
