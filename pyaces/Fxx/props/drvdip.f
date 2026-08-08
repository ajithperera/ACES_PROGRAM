      SUBROUTINE DRVDIP(SCR,MAXCOR,NORBS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C



































































































































































































C
      DOUBLE PRECISION DIP(3), xprop(19)

      DOUBLE PRECISION BUF(600)
      INTEGER IBUF(600)
C
      CHARACTER LABELS(3)
      CHARACTER*80 FNAME 
      DIMENSION SCR(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /FILES/ LUOUT,MOINTS
      DATA LABELS /'X','Y','Z'/
      I010=1+NORBS*NORBS
      CALL SEEKLB_PROPS('     X  ',IERR,0)
      IF(IERR.NE.0)RETURN
      CALL COMPPR(DIP(1),SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB_PROPS('     Y  ',IERR,1)
      CALL COMPPR(DIP(2),SCR,SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB_PROPS('     Z  ',IERR,1)
      CALL COMPPR(DIP(3),SCR, SCR(I010),NORBS,.TRUE.)
      x = SDOT(norbs*norbs, scr, 1, scr, 1)
      y = SDOT(norbs*norbs, scr(i010), 1, scr(i010), 1)
      Z = SDOT(norbs*norbs, scr(i010), 1, scr, 1)
      write(6,"(a,F20.10)") ' @props/drvdip, <Z> ',  Z
      write(6,*) ' @drvdip, Trace Z-integrals ', y
      call output(SCR(I010), 1, norbs, 1, norbs, norbs, norbs, 1)
C      write(6,*) ' @drvdip, Trace dens ', x
C      call output(SCR, 1, norbs, 1, norbs, norbs, norbs, 1)
      WRITE(LUOUT,100)
      WRITE(LUOUT,101)(LABELS(I),DIP(I),I=1,3)
100   FORMAT(T3,' Components of electric dipole moment ')
101   FORMAT((T6,3(1X,A3,' = ',F15.10)))
c
c write properties to jobarc
c
      if (iflags2(163) .ne. 0) then
      call getrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      do i = 1, 3
         xprop(i) = dip(i)
      enddo
      call putrec(20,'JOBARC','MULTMOM ',
     $              19*iintfp, xprop)
      endif

      RETURN
      END
