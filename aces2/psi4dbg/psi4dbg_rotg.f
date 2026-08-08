










      SUBROUTINE Psi4dbg_ROTG(Grd,Grd_oo,Grd_vv,Grd_ov,Grd_vo,
     +                     Grad_stat,Icore,Maxcor,Lenoo,Lenvv,Lenvo,
     +                     Nocci,Nvrti,Nbas,Iuhf)
C
C THIS ROUTINE ROTATES THE MOLECULAR ORBITALS BY THE T1-LIKE
C  AMPLITUDES ON LISTS 3,90 (AND 4,90 FOR UHF) CALCULATIONS.  THIS
C  ROUTINE GENERATES A NEW SET OF ORTHONORMALIZED MOLECULAR ORBITALS.
C
C Extended to work with dropmo. Ajith Perera, 07/2005.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DOUBLE PRECISION GRD(NBAS,NBAS)
      DOUBLE PRECISION GRD_OO(LENOO)
      DOUBLE PRECISION GRD_VV(LENVV)
      DOUBLE PRECISION GRD_VO(LENVO)
      DOUBLE PRECISION GRD_OV(LENVO)
      DOUBLE PRECISION GRAD_STAT(6)
      DOUBLE PRECISION DSUM,DDOT
      CHARACTER*1 ISP(2)
      DIMENSION ICORE(MAXCOR)
      DIMENSION NOCC_EXPND(8,2),NVRT_EXPND(8,2),
     &          NOCCO_EXPND(2),NVRTO_EXPND(2),
     &          NBF4IRREP(8)
      COMMON /SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(4)
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM /-1.0/
      DATA ISP /'A','B'/
C
      Write(6,"(a)") "--------------------pCCD_rotg--------------------"
      CALL GETREC(20, "JOBARC", 'NUMDROPA', 1, NDROPMO)
      NBAS2=NBAS*NBAS
      IF (NDROPMO.GT.0) THEN
         NBAS_DROP = NOCCO(1)+NVRTO(1)
         CALL FILL_DROP_MOS(NIRREP, IUHF, NOCC_EXPND, NVRT_EXPND,
     &                      NBAS_FULL, NOCCO_EXPND, NVRTO_EXPND)
         NBAS  = NBAS_FULL
         NBAS2 = NBAS_FULL*NBAS_FULL
      END IF
C
C NOCC_EXPND, POP = The # of occ per irrep in the full and reduced spaces.
C NVRT_EXPND, VRT = The # of vrt per irrep in the full and reduced spaces.
C NOCCO_EXPND, NOCCO = The # of occ in the full and reduced spaces.
C NVRTO_EXPND, NVRTO = The # of vrt in the full and reduced spaces.
C NBAS = The number of basis functions in the full space.
C
      DO 10 ISPIN=1,1+IUHF
       NOCC=NOCCO(ISPIN)
       NVRT=NVRTO(ISPIN)
       IF (NDROPMO.GT.0) THEN
          NOCC=NOCCO_EXPND(ISPIN)
          NVRT=NVRTO_EXPND(ISPIN)
       END IF
       I000=1
       I010=I000+NBAS2*IINTFP
       I020=I010+NBAS2*IINTFP
       I030=I020+NOCC*NVRT*IINTFP
       I040=I030+NOCC*NVRT*IINTFP
       IOFF=IINTFP*(NBAS*NOCC)

C GET MO COEFFICIENTS AND T1 VECTOR
C (Always read the MO vectors in the full space.)
C
       IF (NDROPMO.GT.0) THEN
          CALL GETREC(20,'JOBARC','SCFEVC'//ISP(ISPIN)//'0',
     &                NBAS2*IINTFP,ICORE(I000))
       ELSE
          CALL GETREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),
     &                NBAS2*IINTFP,ICORE(I000))
       END IF
       CALL SCOPY (NBAS2,ICORE(I000),1,ICORE(I010),1)

       MEMLEFT = MAXCOR-I040 
       CALL PSI4DBG_NR_SEARCH(GRD,GRD_OO,GRD_VV,GRD_VO,GRD_OV,GRAD_STAT,
     &                        LENOO,LENVV,LENVO,NBAS,NOCC,NVRT,
     &                        ICORE(I040),MEMLEFT)

       write(6,"(a)") "The U"
       call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
C DO ORBITAL ROTATIONS 
C
       CALL XGEMM ('N','N',NBAS,NBAS,NBAS,ONE,ICORE(I000),NBAS,GRD,
     &             NBAS,ZILCH,ICORE(I010),NBAS)
C
       CALL PUTREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),NBAS2*IINTFP,
     &             ICORE(I010))
C For norm conserving orbital rotations, the extra orthonormalization is
C not necessary. But to be safe it is good to monitor it and have an error
C exit.
       CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS2*IINTFP,ICORE(I000))
       CALL GETREC(20,'JOBARC','SCFEVEC'//ISP(ISPIN),NBAS2*IINTFP,
     &             ICORE(I020))
       call xgemm("T", "N",nbas,nbas,nbas,one,icore(i020),nbas,
     &             icore(i000),nbas,zilch,icore(i010),nbas)

       call xgemm("N", "N",nbas,nbas,nbas,one,icore(i010),nbas,
     &             icore(i020),nbas,zilch,icore(i000),nbas)
       call output(icore(i000), 1, nbas, 1, nbas, nbas, nbas, 1)
       Dsum = Ddot(Nbas*Nbas,icore(i000),1,icore(i000),1)
       If ((Int(Dsum) - Dble(Nbas)) .GT. 1.0D-12) Then
           Write(6,"(2a)") " The norm condition (C^(t)SC=1) is",
     &                     " violated and can not proceed."
           Write(6,*)
           Call Errex
       Endif

       IF(IUHF.EQ.0)THEN
        CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS2*IINTFP,ICORE(I020))
       ENDIF

10    CONTINUE

C
      RETURN
      END
