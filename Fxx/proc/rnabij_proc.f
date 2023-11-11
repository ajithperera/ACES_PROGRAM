






































































































































































































      SUBROUTINE RNABIJ_PROC(ICORE,MAXCOR,IUHF,TYPE)
C
C THIS SUBROUTINE DRIVES THE FORMATION OF THE RESORTED <Ab|Ij>
C  AND <AB||IJ> INTEGRALS OR T(Ab;Ij) AND T(AB;IJ) AMPLITUDES AND 
C  WRITES THEM TO DISK.
C
C    INPUT - TYPE (CHARACTER*1): 'W' FOR INTEGRAL RESORTS.
C                                'T' FOR T2 VECTOR RESORTS.
C                                'L' FOR L2 VECTOR RESORTS.
C                                'R' FOR T2 INCREMENT RESORTS.
C                                'C' FOR T2 INCREMENT RESORTS FROM
C                                    A DIFFERENT SET OF LISTS
C
C    FOR RHF, THE AI-bj, Aj-bI AND AI-BJ PACKED LISTS ARE WRITTEN.
C    FOR UHF(ROHF), THE bj-AI, bI-Aj AND ai-bj PACKED LISTS ARE ALSO
C        WRITTEN.
C
C   The 'C' option is only used for iterative MBPT, such as in ROHF
C   with standard orbitals.
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 SPCASE(2),SPCAS2,RERDTP
      CHARACTER*1 TYPE
      LOGICAL DRCCD
      DIMENSION ICORE(MAXCOR),LIST(6)
      DIMENSION ISIZOO(8), ISIZVV(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SHIFT/ ISHIFT 

C Parameters for parametrized CC

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      DATA SPCASE /'AABB','BBAA'/

      NNM1O2(I)=(I*(I-1))/2
      NOCA=NOCCO(1)
      NOCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      IF(TYPE.EQ.'W')THEN
       LIST(1)=18 + ISHIFT 
       LIST(2)=17 + ISHIFT 
       LIST(3)=21 + ISHIFT 
       LIST(4)=22 + ISHIFT 
       LIST(5)=19 + ISHIFT 
       LIST(6)=20 + ISHIFT 
       REFLST=16 + ISHIFT 
       RERDTP='AIBJ'
      ELSEIF(TYPE.EQ.'R')THEN
       LIST(1)=37
       LIST(2)=36
       LIST(3)=39
       LIST(4)=38
       LIST(5)=34
       LIST(6)=35
       REFLST=63
       RERDTP='AJBI'
      ELSEIF(TYPE.EQ.'T')THEN
       LIST(1)=37
       LIST(2)=36
       LIST(3)=39
       LIST(4)=38
       LIST(5)=34
       LIST(6)=35
       REFLST=46
       RERDTP='AJBI'
      ELSEIF(TYPE.EQ.'C')THEN
       LIST(1)=37
       LIST(2)=36
       LIST(3)=39
       LIST(4)=38
       LIST(5)=34
       LIST(6)=35
       REFLST=96
       RERDTP='AJBI'
      ELSEIF(TYPE.EQ.'L') THEN
       LIST(1)=137
       LIST(2)=136
       LIST(3)=139
       LIST(4)=138
       LIST(5)=134
       LIST(6)=135
       REFLST=146
       RERDTP='AJBI'
      ENDIF
      ISZTOT=ISYMSZ(13,14)
      NSCRSZ=NOCA*NOCB+NVRTA*NVRTB+NVRTA*NOCA+NVRTB*NOCB
      ISZTAR=ISYMSZ(9,10)
      ISZTAR2=ISYMSZ(11,12)
      I000=1
      I010=I000+IINTFP*MAX(ISZTOT,ISZTAR,ISZTAR2)
      I020=I010+IINTFP*MAX(ISZTOT,ISZTAR,ISZTAR2)
      I030=I020+NSCRSZ
      IF(I030.GT.MAXCOR)CALL INSMEM('RNABIJ',I030,MAXCOR)
C
C WRITE AI-bj AND Aj-bI ORDERINGS AND bj-AI AND bI-Aj IF UHF.
C
      DO 20 I=1,1+IUHF
       CALL GETALL(ICORE(I000),ISZTOT,1,REFLST)
       CALL SST002(ICORE(I000),ICORE(I010),ISZTOT,ISZTAR,ICORE(I020),
     &             SPCASE(I))
       CALL PUTALL(ICORE(I010),ISZTAR,1,LIST(I))
       CALL SSTRNG(ICORE(I010),ICORE(I000),ISZTAR,ISZTAR2,ICORE(I020),
     &             SPCASE(I))
       CALL PUTALL(ICORE(I000),ISZTAR2,1,LIST(2+I))
20    CONTINUE
C
C WRITE AI-BJ INTEGRALS (ai-bj FOR UHF IF ISPIN=2)
C
      SPCAS2='AAAA'
      DO 30 ISPIN=1,1+IUHF
       ISZTOT=ISYMSZ(ISPIN,ISPIN+2)
       ISZTAR=ISYMSZ(8+ISPIN,8+ISPIN)
       NSCRSZ=NNM1O2(NVRTO(ISPIN))+NNM1O2(NOCCO(ISPIN))+
     &        NVRTO(ISPIN)*NOCCO(ISPIN)
       I000=1
       I010=I000+IINTFP*MAX(ISZTOT,ISZTAR)
       I020=I010+IINTFP*ISZTAR
       I030=I020+NSCRSZ
       IF(I030.GT.MAXCOR)CALL INSMEM('RNABIJ',I030,MAXCOR)
       CALL GETALL(ICORE(I000),ISZTOT,1,REFLST-3+ISPIN)
       CALL SST003(ICORE(I000),ICORE(I010),ISZTOT,ISZTAR,ICORE(I020),
     &             SPCAS2,RERDTP)
       CALL PUTALL(ICORE(I010),ISZTAR,1,LIST(4+ISPIN))
       SPCAS2='BBBB'
30    CONTINUE

      DRCCD = (IFLAGS(2) .EQ. 49)
      IF ((ISPAR .AND. COULOMB) .OR. DRCCD) THEN

      Write(*,*)
      Write(*,"(a)") "The Coulomb only AI,BJ orderd <AB|IJ> is formed"
      Write(*,*)
      SPCAS2 ='AAAA'
      IRREPX = 1
      DO 40 ISPIN=1+IUHF, 1, -1
       ISZTOT=ISYMSZ(18+ISPIN,20+ISPIN)
       ISZTAR=ISYMSZ(8+ISPIN,8+ISPIN)
       NSCRSZ=NVRTO(ISPIN)*NVRTO(ISPIN)+NOCCO(ISPIN)*NOCCO(ISPIN)+
     &        NVRTO(ISPIN)*NOCCO(ISPIN)
       I000=1
       I010=I000+IINTFP*MAX(ISZTOT,ISZTAR)
       I020=I010+IINTFP*ISZTAR
       I030=I020+NSCRSZ
       IF(I030.GT.MAXCOR)CALL INSMEM('RNABIJ',I030,MAXCOR)
CSSS       CALL GETALL(ICORE(I000),ISZTOT,1,113+ISPIN)
       CALL GETALL(ICORE(I000),ISZTOT,1,213+ISPIN)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISZTOT,VRT(1,ISPIN),
     &             VRT(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),ICORE(I020),
     &             IRREPX,"1324")
       CALL PUTALL(ICORE(I010),ISZTAR,1,118+ISPIN)
       SPCAS2='BBBB'


40    CONTINUE
      ENDIF 
      RETURN
      END
