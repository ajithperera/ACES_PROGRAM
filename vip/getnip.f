










      SUBROUTINE GETNIP(NIP,SIRREP,ISPIN,IUHF)
C
      IMPLICIT INTEGER (A-Z)
      logical lefthand,  singonly, DROPCORE
      DIMENSION NS(8)
C
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      COMMON/IPCALC/LEFTHAND,SINGONLY,DROPCORE
C
C FIRST FIND THE DIMENSION OF THE MATRICES
C
      CALL IZERO(NS, 8)
      NS(SIRREP) = 1
      IF (IUHF.EQ.0) THEN
         CALL GETLEN_HHP(LENAB,POP(1,1),POP(1,1),VRT(1,1),NS(1))
         LENAA = 0
      ELSE
         MSPIN = 3 - ISPIN
         CALL GETLEN_HHP2(LENAA,IRPDPD(1,ISPIN+2),VRT(1,ISPIN),
     +                     NS(1))
         CALL GETLEN_HHP2(LENAB,IRPDPD(1,14),VRT(1,MSPIN),NS(1))
      ENDIF
      NIP = POP(SIRREP,ISPIN) + LENAA + LENAB
C
      RETURN
      END
