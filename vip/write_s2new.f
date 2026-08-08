










       SUBROUTINE WRITE_S2NEW(ICORE,MAXCOR,IUHF,ISPIN)



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
C S(m=MIXSPIN,aj=ISPIN)
C
      DO MIXSPIN = 1, IUHF+1 
         LISTS2EX = LS2OUT(ISPIN, MIXSPIN + 1 - IUHF)

         IF ((IUHF.NE.0) .AND. (MIXSPIN.EQ.ISPIN)) THEN
            CALL ASSYMALL(ICORE(I020),LENS,NUMSZS,IRPDPD(1,18+ISPIN),
     &                    VRT(1,ISPIN),ICORE(I030),MAXCOR-I030+1)
            CALL GETEXP2_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     &                       LISTS2EX,VRT(1,ISPIN),
     &                       IRPDPD(1,ISYTYP(1,LISTS2EX)))
       ELSE
          CALL GETALLS2_HHP(ICORE(I010),LENS,POP(1,SOUTSPIN),NS(1),
     &                      1,LISTS2EX)
       ENDIF
C
       CALL SAXPY (LENS,ONE,ICORE(I020),1,ICORE(I010),1)
C
       IF ((IUHF.NE.0).AND. (SOUTSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ_HHP(ICORE(I010),LENS,NUMSZS,IRPDPD(1,20+ISPIN),
     &                    LISTS2EX, VRT(1, ISPIN),
     &                    IRPDPD(1, ISYTYP(1, LISTS2EX)),ICORE(I020),
     &                     MAXCOR-I020+1)
       ELSE
          CALL PUTALLS2_HHP(ICORE(I010),LENS,VRT(1,ISPIN),NS(1),1,
     $                      LISTS2EX)
       ENDIF

      ENDDO
C
      RETURN
      END

