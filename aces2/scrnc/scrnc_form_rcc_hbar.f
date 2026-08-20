













































































































































































































      SUBROUTINE SCRNC_FORM_RCC_HBAR(ICORE,MAXCOR,IUHF,IMULT)
C
C THIS ROUTINE FORMS THE SINGLE-SINGLE BLOCK OF EITHER THE REGULAR
C OR THE CCSD EFFECTIVE HAMILTONIAN
C
C       H(ai,bj) = F(ab) d(ij) - F(ij) d(ab) - W(ai,bj)
C
C WHERE F AND W ARE THE ONE AND TWO PARTICLE PARTS OF EXP(-T) H EXP(T).
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      LOGICAL NONHF 
      LOGICAL HEFF,TDA,EOM_SFRCCD,EOM_SFDRCC
      DOUBLE PRECISION ONE,ONEM,TWO,HALF
      DIMENSION ICORE(MAXCOR),IOFFO(2),IOFFV(2)
C
      COMMON/FLAGS/IFLAGS(100)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C 
      DATA ONE,ONEM /1.0D0,-1.0D0/
      DATA TWO,HALF /2.0D0,0.50D0/
      Write(6,"(1x,a)") "-----Entering scrnc_form_rcc_hbar------"

      HEFF  = .FALSE.
      TDA   = .FALSE.
      CIS   = .FALSE.
      RPA   = .TRUE.
      NONHF =  (IFLAGS(38)  .GT. 0)

C Occ-Occ block

      I000=1
      I010=I000+NFMI(1)*IINTFP
      IF(IUHF.NE.0)THEN
       I020=I010+NFMI(2)*IINTFP
      ELSE
       I020=I010
      ENDIF

C Vrt-Vrt  block

      I030=I020+NFEA(1)*IINTFP
      IF(IUHF.NE.0)THEN
       I040=I030+NFEA(2)*IINTFP
      ELSE
       I040=I030
      ENDIF

      IF (CIS) THEN
C This simply do CIS and is usefull for debugging.
          HEFF = .FALSE.
          TDA  = .TRUE. 

      ELSEIF(RPA) THEN

C This is drRPA or RPA, but using Hbar from drCCD or rCCD with HF diagonals.
C These methods are designated as EOM(Sf)-RCCD or  EOM(Sf)-DRCCD

          HEFF = .FALSE.
          TDA  = .FALSE. 
      ENDIF 

      IF (HEFF) THEN

C If Hbar(i,j) and Hbar(a,b) must be used; load them here. 
         
          CALL GETLST(ICORE(I000),1, 1, 1, 1, 91)
          CALL GETLST(ICORE(I020),1, 1, 1, 1, 92)
          IF (IUHF .EQ. 1) THEN
            CALL GETLST(ICORE(I010),1, 1, 1, 2, 91)
            CALL GETLST(ICORE(I030),1, 1, 1, 2, 92)
          ENDIF
      ELSE 

C For NON-HF we need to F(a,b) and F(i,j) elements. For EOM(S)-RCC
C or DRCC methods these have already been added when we cnstruct
C the Hbar(a,b) and Hbar(i,j). But for EOM(SF)-RCC or -DRCC we need to
C add them here.

          IF (NONHF) THEN
             IF (RPA .OR. EOM_SFRCCD .OR. EOM_SFDRCCD) THEN
                CALL GETLST(ICORE(I000),1, 1, 1, 3, 91)
                CALL GETLST(ICORE(I020),1, 1, 1, 3, 92)
                IF (IUHF .EQ. 1) THEN
                  CALL GETLST(ICORE(I010),1, 1, 1, 4, 91)
                  CALL GETLST(ICORE(I030),1, 1, 1, 4, 92)
                ENDIF
             ELSE
                CALL IZERO(ICORE(I000),I040-1)
             ENDIF 
          ELSE
             CALL IZERO(ICORE(I000),I040-1)
C 2026-08-20 (Ajith's DCCSD-precedent port, see SCRNC_DEBUG_LOG.md):
C this Hbar(i,j)/Hbar(a,b) one-particle block was PURE ZERO here for
C RHF RCCD/DRCCD -- i.e. the naive bare-Fock-only Hbar the JCP2020
C (d)rCCD paper explicitly says is wrong. Fill it with the Lagrangian-
C derived quadratic Fmi/Fae correction (Scrnc_rcc_lagrangian_fix,
C called once from Scrnc_driver, stores into list 91/92 sub-index 9).
C TEMPORARILY DISABLED to isolate Part-1-only (two-particle) effect,
C see rcc_lagrangian_fix.F. See SCRNC_DEBUG_LOG.md.
C            IF (IUHF .EQ. 0) THEN
C               CALL GETLST(ICORE(I000),1,1,1,9,91)
C               CALL GETLST(ICORE(I020),1,1,1,9,92)
C            ENDIF
          ENDIF
       ENDIF

C The Diagonal elements. First alpha occupied and followed by virtuals 

       NBAS=NOCCO(1)+NVRTO(1)
       CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*NBAS,ICORE(I040))

       IOFFE =I040
       IOFFTO=I000
       IOFFTV=I020

       DO 10 IRREP=1,NIRREP
          NOCC=POP(IRREP,1)
          CALL SAXPY(NOCC,ONE,ICORE(IOFFE),1,ICORE(IOFFTO),NOCC+1)
          IOFFE =IOFFE +NOCC*IINTFP
          IOFFTO=IOFFTO+NOCC*NOCC*IINTFP
10     CONTINUE

       DO 11 IRREP=1,NIRREP
          NVRT=VRT(IRREP,1)
          CALL SAXPY(NVRT,ONE,ICORE(IOFFE),1,ICORE(IOFFTV),NVRT+1)
          IOFFE =IOFFE +NVRT*IINTFP
          IOFFTV=IOFFTV+NVRT*NVRT*IINTFP
11     CONTINUE
C
       IF(IUHF.EQ.1)THEN
         CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*NBAS,ICORE(I040))
         IOFFE =I040
         IOFFTO=I010
         IOFFTV=I030

        DO 12 IRREP=1,NIRREP
           NOCC=POP(IRREP,2)
           CALL SAXPY(NOCC,ONE,ICORE(IOFFE),1,ICORE(IOFFTO),NOCC+1)
           IOFFE =IOFFE+NOCC*IINTFP
           IOFFTO=IOFFTO+NOCC*NOCC*IINTFP
12      CONTINUE

        DO 13 IRREP=1,NIRREP
           NVRT=VRT(IRREP,2)
           CALL SAXPY(NVRT,ONE,ICORE(IOFFE),1,ICORE(IOFFTV),NVRT+1)
           IOFFE =IOFFE+NVRT*IINTFP
           IOFFTV=IOFFTV+NVRT*NVRT*IINTFP
13        CONTINUE
       ENDIF
C
      IOFFO(1)=I000
      IOFFO(2)=I010
      IOFFV(1)=I020
      IOFFV(2)=I030

C Read in Hbar(ib,aj) elements; ordered as (bi,aj)

      DO ISPIN=1,1+IUHF
         IF(IUHF.EQ.0)THEN
C 2026-08-19: ISIZE must match whichever list (56 singlet, 54 triplet)
C this Imult pass actually reads/writes -- was hardcoded to list 56's
C type regardless of IMULT. Numerically harmless for H2/STO-3G (both
C lists happen to share the same size there) but wrong in general, and
C risky now that IMULT=2 also WRITES via this same ISIZE (see the
C LISTDUMP fix below) -- an undersized/oversized PUTALL to list 54
C would silently corrupt it. See SCRNC_DEBUG_LOG.md.
           LISTHBAR = 56
           IF (IMULT .EQ. 2) LISTHBAR = 54
           ISIZE=ISYMSZ(ISYTYP(1,LISTHBAR),ISYTYP(1,LISTHBAR))
           I050=I040+ISIZE*IINTFP
           I060=I050+ISIZE*IINTFP
           IEND=I060
           IF (IEND .GT. MAXCOR) CALL INSMEM("RCL_HEFFXC",
     &                                       IEND,MAXCOR)
           IF (TDA) THEN
              CALL GETALL(ICORE(I040),ISIZE,1,23)
              CALL VMINUS(ICORE(I040),ISIZE)
              CALL GETALL(ICORE(I050),ISIZE,1,18)
              CALL SAXPY (ISIZE,ONE,ICORE(I050),1,ICORE(I040),1)
           ELSE 
              IF (IMULT .EQ. 1) THEN
                 CALL GETALL(ICORE(I040),ISIZE,1,56)
              ELSE IF (IMULT .EQ. 2) THEN
                 CALL GETALL(ICORE(I040),ISIZE,1,54)
              ENDIF
C 2026-08-20 TEST (Ajith): raw W(ai,bj) straight off list 56, BEFORE
C the (ai,bj)->(ab,ij) reorder and SCRNC_ADDONEH Fock-diagonal add
C below -- same (ai,bj) index convention as the final "drCCD Hbar
C matrix" dump in scrnc_form_rcc_hbar_inv.F, so row/col 9,10 here are
C directly comparable to that dump's row/col 9,10. See SCRNC_DEBUG_LOG.md.
              IF (IMULT .EQ. 1) THEN
                MATDIMX = NINT(SQRT(DBLE(ISIZE)))
                Write(6,"(a)") " TEST: raw W(ai,bj) pre-Fock-diag,"//
     +                         " list 56 (singlet)"
                call output(ICORE(I040),1,MATDIMX,1,MATDIMX,
     +                      MATDIMX,MATDIMX,1)
              ENDIF
           ENDIF
         ELSE
           IF (TDA) THEN
               LISTAA=22+ISPIN
           ELSE 
               LISTAA=53+ISPIN
           ENDIF 
           ISIZE=ISYMSZ(ISYTYP(1,LISTAA),ISYTYP(1,LISTAA))
           I050=I040+ISIZE*IINTFP
           I060=I050+ISIZE*IINTFP
           IEND=I060 
           IF (IEND .GT. MAXCOR) CALL INSMEM("RCL_HEFFXC",
     &                                        IEND,MAXCOR)
           CALL GETALL(ICORE(I040),ISIZE,1,LISTAA)
           IF (TDA) CALL VMINUS(ICORE(I040),ISIZE)
          enDIF

C Reorder; (ai,bj) -> (ab,ij)

          CALL SSTGEN(ICORE(I040),ICORE(I050),ISIZE,
     &                VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &                POP(1,ISPIN),ICORE(I060),1,'1324')

C Add the one particle contributions

          DISSZ1=IRPDPD(1,18+ISPIN)
          NUMDS1=IRPDPD(1,20+ISPIN)
          CALL SCRNC_ADDONEH(ICORE(I050),ICORE(I040),ISIZE,
     &                       ICORE(IOFFV(ISPIN)),
     &                       ICORE(IOFFO(ISPIN)),
     &                       DISSZ1,NUMDS1,ISPIN)

C Reorder (ab,ij) -> (ai,bj)

          CALL SSTGEN(ICORE(I050),ICORE(I040),ISIZE,
     &                VRT(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &                POP(1,ISPIN),ICORE(I060),1,'1324')
     
          IF (TDA) THEN
            LISTDUMP = 42
          ELSE IF (IMULT .EQ. 2) THEN
C 2026-08-19: the completed (Fock-diagonal-added) triplet Hbar must be
C written back to its OWN list (54), matching the raw-W read at line
C ~173 above -- previously this always wrote to list 56 regardless of
C IMULT, clobbering the singlet Hbar and leaving list 54 stuck at its
C raw, Fock-diagonal-incomplete state forever. See SCRNC_DEBUG_LOG.md.
            LISTDUMP = 54
          ELSE
            LISTDUMP = 56
          ENDIF
          IF (IUHF .NE. 0) LISTDUMP = LISTAA
          CALL PUTALL(ICORE(I040),ISIZE,1,LISTDUMP)

      ENDDO 

      RETURN
      END
