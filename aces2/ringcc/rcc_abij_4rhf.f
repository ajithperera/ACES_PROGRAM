










      SUBROUTINE RCC_ABIJ_4RHF(ICORE,MAXCOR,IUHF)

      IMPLICIT INTEGER (A-Z)

      DIMENSION ICORE(MAXCOR)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),IJUNK(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)

C#ifdef _NOSKIP
      IRREPX = 1
      NSIZE_AJBI  = IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))
      NSIZE_ABIJ  = IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))

      I000  = 1
      I010  = I000 + NSIZE_AJBI * IINTFP
      I020  = I010 + NSIZE_ABIJ * IINTFP
      I030  = I020 + NSIZE_ABIJ
      IEND  = I030
      IF (IEND .GE. MAXCOR) CALL INSMEM("rcc_initin",
     +                                        IEND,MAXCOR)
C Read <Aj|bi> integrals 

      CALL GETALL(ICORE(I000),NSIZE_AJBI,IRREPX,21)

C <Aj|bi> -> <Ab|ji>

      CALL SSTGEN(ICORE(I000),ICORE(I010),NSIZE_ABIJ,VRT(1,1),
     +            POP(1,2),VRT(1,2),POP(1,1),ICORE(I020),IRREPX,
     +           "1324")

C These are exchange integrals of <Ab|Ij> (RHF only) 

C Form the triplet part of the spin-adapted <Ab|ji> by -<Ab|ji>
      
      CALL DSCAL(NSIZE_ABIJ,-1.0D0,ICORE(I010),1)
      CALL PUTALL(ICORE(I010),NSIZE_ABIJ,IRREPX,198)
      call checksum("List-198:",ICORE(I010),NSIZE_ABIJ)

c Form the singlet part of the spin-adapted <Ab|ij> by 2<Ab|Ij>-<Ab|jI>

      CALL GETALL(ICORE(I000),NSIZE_ABIJ,IRREPX,16)
      CALL DSCAL(NSIZE_ABIJ,2.0D0,ICORE(I000),1)
      CALL DAXPY(NSIZE_ABIJ,1.0D0,ICORE(I010),1,ICORE(I000),1)
      CALL PUTALL(ICORE(I000),NSIZE_ABIJ,IRREPX,197)
      call checksum("List-197:",ICORE(I000),NSIZE_ABIJ)
C#endif

        nsize = idsymsz(1,isytyp(1,198),isytyp(2,198)) 
        i000=1
        i010=i000 + nsize
        i020=i010 + nsize 
        call getall(icore(i000),nsize,1,198)
        call checksum("198",icore(i000),nsize)
        call getall(icore(i010),nsize,1,197)
        call checksum("197",icore(i010),nsize)
      
      RETURN
      END
