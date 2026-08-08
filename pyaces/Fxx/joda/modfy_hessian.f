
      SUBROUTINE MODFY_HESSIAN(DIAGHES, HESMOD, HES, QSTLST_TANGENT,
     &                         SCRATCH, EIGVALUE, WEIGHT, NOPT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)









































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/





      DIMENSION HESMOD(NOPT, NOPT), QSTLST_TANGENT(NOPT),
     &          SCRATCH(NOPT), DIAGHES(NOPT, NOPT),
     &          HES(NOPT, NOPT)
      DATA ONE /1.0D0/, TWO /2.0D0/
C
C Note that at the moment SCRATCH array contains the HxT. Let's label
C that product as A vector. As we can see the A vector is an arbitrary
C direction in the space of the eigenvectors of H.
C Impetus for this formula is from Ayala & Schlegel, JCP, 375, 1997.
C Eqn. 7a. 
 
         SCALED_WEIGHT = TWO*EIGVALUE*WEIGHT
        
         CALL XGEMM('T','N', NOPT, NOPT, 1, 1.0D0, QSTLST_TANGENT,
     &               1, SCRATCH, 1, 0.0D0, DIAGHES, NOPT)
         CALL TRANSP(DIAGHES, HESMOD, NOPT, NOPT)
         CALL DAXPY(NOPT*NOPT, 1.0D0, HESMOD, 1, DIAGHES, 1) 
         CALL DSCAL(NOPT*NOPT, -WEIGHT, DIAGHES, 1)
         CALL DAXPY(NOPT*NOPT, 1.0D0, DIAGHES, 1, HES, 1)
         CALL XGEMM('T','N', NOPT, NOPT, 1, SCALED_WEIGHT, 
     &               QSTLST_TANGENT, 1,
     &               QSTLST_TANGENT, 1, 1.0D0, HES, NOPT)

      Write(6,*)
      Write(6,"(a)") "The modified. Hessian"
      Call output(HES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
C
C Diagonalize the modified Hessian and print the eigenvalues and
C vectors. These new vectors and the values dictate the climbing
C phase of the search.
C
         CALL DCOPY(NOPT*NOPT, HES, 1, HESMOD, 1)
         CALL EIG(HESMOD, DIAGHES, NOPT, NOPT, 1)

      Write(6,"(a)") "The Eigen vectors of modified. Hessian"
      Call output(DIAGHES, 1, NOPT, 1, NOPT, NOPT, NOPT, 1)
      Write(6,"(a)") "The Eigen values of modified Hessian"
      Write(6, "(3(1x,F12.6)))") (hesmod(i,i), i=1,nopt)
C
C
C This block of code is obsoleted. 

      RETURN
      END

