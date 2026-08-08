

























































































































































































































      SUBROUTINE DDMPLGSS(SCR, MAXCOR, IUHF, NLIST, NAME)
C
C  This routine figures out if it is safe to dump a L
C vector to use as a later guess.  If so, it loads R into core and calls 
c DMPRGSS.

      IMPLICIT NONE
C
      INTEGER MAXCOR, IUHF, LEngth, NLIST
      DOUBLE PRECISION SCR(MAXCOR), S
      CHARACTER *8 NAME
C
      INTEGER IFLAGS
      COMMON /FLAGS/ IFLAGS(100)
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
      INTEGER POP, VRT, NT, NFMI, NFEA, i
      COMMON /SYM/ POP(8,2), VRT(8,2), NT(2), NFMI(2), NFEA(2)
      INTEGER IRPDPD, ISYTYP, ID

      INTEGER IRREPR, IRREPL, DISSZT, NUMDST, DISSZF, NUMDSF
      INTEGER IOFFF, IOFFT, NSIZE, I060, IRREPX, ITMP1, ITMP2, IDIS

      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
C
      INTEGER ISYMSZ
      CHARACTER *4 FPGRP, CPGRP
C
      INTEGER I000, I010, I020, I030, I040, I050, I070
C
C Dump the L vector.
C
      I000 = 1
      I010 = I000 + NT(1)
      I020 = I010 + IUHF*NT(2)
      I030 = I020 + IUHF*ISYMSZ(2,4)
      I040 = I030 + ISYMSZ(1,3)
      I050 = I040 + ISYMSZ(13,14)

      IF (I050 .GT. MAXCOR) THEN
         WRITE(6,1000)
 1000    FORMAT(T3, '@DDMPLGSS-I, not enough memory to dump L',
     &              ' as initial guess')
      ELSE

         I000 = 1
         CALL GETLST(SCR(I000), 1, 1, 1, 1, 90+NLIST)
         IF (IUHF .NE. 0) THEN
            CALL GETLST(SCR(I010), 1, 1, 1, 2, 90+NLIST)
            CALL GETALL(SCR(I020),
     &                  ISYMSZ(ISYTYP(1,45+NLIST), ISYTYP(2,45+NLIST)),
     &                  1, 45+NLIST)
         ENDIF

         CALL GETALL(SCR(I030),
     &               ISYMSZ(ISYTYP(1,44+NLIST), ISYTYP(2,44+NLIST)),
     &               1, 44+NLIST)
C 
         CALL GETALL(SCR(I040),
     &               ISYMSZ(ISYTYP(1,46+NLIST), ISYTYP(2,46+NLIST)),
     &               1, 46+NLIST)

         CALL DMPLGSS(SCR(I000), I050-I000, NAME)
      ENDIF
C
      RETURN
      END
