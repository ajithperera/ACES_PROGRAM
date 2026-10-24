











      SUBROUTINE QSTLST_OR_EVEC(QSTLST_TANGENT, GRDMOD, HESMOD, DIAGHES,
     &                          HES, SCRATCH, IMODE, QSTLST_CLIMB)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL QSTLST_CLIMB
      PARAMETER (THRESHOLD = 5.0D-2, WGHT= 1.0D-3)

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end



      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      DATA HALF, ONE /0.50D0, 1.0D0 /

      DIMENSION HESMOD(NOPT, NOPT), GRDMOD(NOPT), QSTLST_TANGENT(NOPT),
     &          SCRATCH(NX*NX), DIAGHES(NOPT, NOPT),
     &          HES(NOPT, NOPT)

      CALL ZERO(SCRATCH, NOPT)

      Write(6,*) 
      Write(6,"(a)") "The QST/LST tangent in QSTLST_OR_EVEC"
      Write(6,"(3(1x,F12.6))") (QSTLST_TANGENT(I), I= 1, NOPT)
      CALL XGEMM('N','N', NOPT, 1, NOPT, 1.0D0, HES, NOPT,
     &            QSTLST_TANGENT, NOPT, 0.0D0, SCRATCH, NOPT)

      EIGVALUE = DDOT(NOPT, QSTLST_TANGENT, 1, SCRATCH, 1)
      GRDVALUE = DDOT(NOPT, QSTLST_TANGENT, 1, GRDMOD,  1)

      EST_STEP = (EIGVALUE + DSQRT(EIGVALUE**2 + 4.0D0*(GRDVALUE)
     &            **2))*HALF
      EST_STEP =  -GRDVALUE/(EST_STEP - EIGVALUE)

      Write(6,*)
      Write(6,"(a)")"Relevent parameter that is used to decide a climb"
      Write(6,"(a,2(1x,F12.6))") "The EIGVALUE,GRDVALUE : ", EIGVALUE, 
     &                            GRDVALUE
      Write(6,"(a,F12.6)") "The estimated step size     :",EST_STEP
      Write(6,*)
      Write(6,"(a)") "The HT vector"
      Write(6,"(3(1x,F12.6))") (Scratch(i), i=1, NOPT)
C
C Let's check that the estimated step size is greater than 0.05 au.
C If so, then follow the QST or LST tangent (I am interpreting Schlegel's
C "estimated displacement along the tangent vector is greater than the 0.05".
C At this point I am not sure whether he meant the largest absolute
C displacement.
C
C Peng and Schlegel: always climb in steps 1 and 2 (NCYCLE 0 and 1),
C in steps 3 and 4 only if the estimated displacement along the
C tangent exceeds 0.05 au, never afterwards. A zero tangent (not
C defined at this point) is never climbed. 10/2026.
C
      TNRM = DDOT(NOPT, QSTLST_TANGENT, 1, QSTLST_TANGENT, 1)
      QSTLST_CLIMB = (TNRM .GT. 0.5D0) .AND.
     &               ((DABS(EST_STEP).GE.THRESHOLD.AND.NCYCLE.LE.3).OR.
     &                 NCYCLE.LE.1)

      Write(6,*)
      Write(6,"(a,l)") "The climbing Phase if true", QSTLST_CLIMB
C
C A climbing step is taken along the tangent T only (eq. 11 of Peng and
C Schlegel with v = T and b = T'HT). T is put in as mode 1 (IMODE) with
C curvature b; GRDHES(1) becomes T'g in TKSTEP, Lambda(p) follows from
C b and T'g in GETLAMBDA and FLOWEVFTS takes the step along mode 1
C only. (This replaces MODFY_HESSIAN, which made T an approximate
C eigenvector of a modified Hessian.) 10/2026.
C
      IF (QSTLST_CLIMB) THEN
         CALL DCOPY(NOPT, QSTLST_TANGENT, 1, DIAGHES(1,1), 1)
         HESMOD(1,1) = EIGVALUE
         IMODE = 1
      END IF
C
C This is just to be safe in subsequent steps.
C
      CALL ZERO(SCRATCH,NOPT*NOPT)

      RETURN
      END

