










      SUBROUTINE PCCD_LGENINT(ICORE,MAXCOR,IUHF)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TCPU,TSYS,ONE
      DIMENSION ICORE(MAXCOR)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      DATA ONE /1.0D0/

      CALL PCCD_FORMV1(ICORE,MAXCOR,IUHF) 
      CALL PCCD_FORMG1(ICORE,MAXCOR,IUHF,ONE)
      CALL PCCD_FORMG2(ICORE,MAXCOR,IUHF,ONE)

      RETURN
      END
