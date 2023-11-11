










      SUBROUTINE INIT_NOREDNT_LIST(IRREPX, SYTYPL, SYTYPR, LIST)
      IMPLICIT NONE

      integer irrepx, sytypl, sytypr, list

c#include "lists.com"
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

CSSS      DO IRREPR = 1, NIRREP
CSSS         IRREPL = DIRPRD(IRREPR, IRREPX)
CSSS         NCOLS  = IRPDPD(IRREPR, SYTYPR)
CSSS         NROWS  = IRPDPD(IRREPR, SYTYPL)
CSSS         ILEFT  = IRREPR
CSSS         IRIGHT = LIST
CSSS         MOIOSZ(ILEFT, IRIGHT) = NROWS
CSSS         MOIODS(ILEFT, IRIGHT) = NCOLS
CSSS      END DO

      ISYTYP(1, LIST) = SYTYPL
      ISYTYP(2, LIST) = SYTYPR

      RETURN
      END
