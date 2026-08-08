










      SUBROUTINE ABTOAA(ICORE, MAXCOR, IUHF, LISTT2)
C 
C Generate AA list from a AB list.
C 
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)     
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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
C
      IF (IUHF .EQ. 0) THEN
   
         DO IRREP=1, NIRREP
            LISTAB   = LISTT2 + 2 
            LISTAA   = LISTT2
            NUMABCOL = IRPDPD(IRREP,ISYTYP(2,LISTAB))
            NUMABROW = IRPDPD(IRREP,ISYTYP(1,LISTAB))
C
            NUMAACOL = IRPDPD(IRREP,ISYTYP(2,LISTAA))
            NUMAAROW = IRPDPD(IRREP,ISYTYP(1,LISTAA))
C
                I000 = 1
                I010 = I000 + IINTFP*NUMABCOL*NUMABROW
                I020 = I010 + IINTFP*NUMAACOL*NUMAAROW 
                IF (I020 .GE. MAXCOR) CALL INSMEM("@ABTOAA",
     &                                I020, MAXCOR)
C
            CALL GETLST(ICORE(I010), 1, NUMABCOL, 1, IRREP, LISTAB)
            CALL ASSYM(IRREP, POP, NUMABROW, NUMABROW, ICORE(I000), 
     &                 ICORE(I010))
            CALL SQSYM(IRREP, VRT, NUMAAROW, NUMABROW, NUMAACOL, 
     &                 ICORE(I010), ICORE(I000))
            CALL PUTLST(ICORE(I010), 1, NUMAACOL, 1, IRREP, LISTAA)
          ENDDO
      ENDIF
C
      RETURN
      END
