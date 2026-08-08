










      SUBROUTINE CALCEXCP(KSPIN,IUHF,SCR,MAXCOR,H0COLUMN,LIST)
C
C AN EXCITATION PATTERN IS FORMED IN THE SAME FORMAT AS THE S-VECTORS ARE
C STORED IN ROUTINE SOLVE_IPEOM SEE LOADS). THE EXCITATION PATTERN CONTAINS A
C ONE (1.0) FOR THOSE OPERATORS (CONFIGURATIONS) THAT ARE OF INTEREST,
C A ZERO (0.0)  OTHERWHISE. THE PATTERN IS PUT ON COLUMN H0COLUMN OF 
C LISTH0
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR, OCCI, OCCBI, OCC
      DIMENSION SCR(MAXCOR), IOFFVRT(8,2), IOFFPOP(8,2), NDUMS(8)
      LOGICAL LEFTHAND, EXCICORE, SINGONLY, DROPOPEN, STSPLIT

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      COMMON/SINFO/NS(8), SIRREP
      COMMON/IPCALC/LEFTHAND, SINGONLY, DROPCORE
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
      I000 = 1
      I010 = I000 + NSIZEC
C
      DO I=I000,I000+NSIZEC -1
           SCR(I) = 1.0D0
      ENDDO 

      IF (SINGONLY) THEN
         DO  I=I000+POP(SIRREP,KSPIN), I000+NSIZEC-1
               SCR(I) = 0.0D0
         ENDDO 
      ENDIF

      CALL PUTLST(SCR(I000), H0COLUMN,1,1,1,LIST)
C
      RETURN
      END
