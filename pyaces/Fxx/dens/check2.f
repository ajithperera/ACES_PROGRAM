      SUBROUTINE CHECK2(ICORE,A,LENGTH,LISTF,ISPIN,FACT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD
      DIMENSION ICORE(1),A(length)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/ADD/SUM
      return
      E=0.0D+0
      CALL GETLST(ICORE,1,1,1,Ispin,LISTf)
      E=E+SDOT(length,ICORE,1,a,1)
1000  CONTINUE
      write(6,*) 'listf',listf
      write(6,*) 'energy contribution',e
      sum=sum+FACT*e
      write(6,*) ' cumulated e',sum
      return
      end
