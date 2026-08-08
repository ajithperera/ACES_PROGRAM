      SUBROUTINE PCCD_GAMDRV(ICORE,MAXCOR,IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ICORE(MAXCOR)
      LOGICAL pCCD,CCD,LCCD

      COMMON /CALC/ pCCD,CCD,LCCD
C 
      CALL PCCD_FORMV1(ICORE,MAXCOR,IUHF)
      CALL PCCD_FORMX(ICORE,MAXCOR,IUHF)

      write(6,*)
      write(6,"(a)") "Checksums of G(ab,cd)"
      CALL PCCD_GAMMA2(ICORE,MAXCOR,IUHF,0)

      write(6,*)
      write(6,"(a)") "Checksums of G(ij,kl)"
      CALL PCCD_GAMMA3(ICORE,MAXCOR,IUHF)
      write(6,*)
      write(6,"(a)") "Checksums of G(ia,jb)"
      CALL PCCD_GAMMA4(ICORE,MAXCOR,IUHF)

C There is no G(ij,ka),G(ak,ij), G(ab,ci) or G(ci,ab) for pCCD. 

      RETURN
      END
