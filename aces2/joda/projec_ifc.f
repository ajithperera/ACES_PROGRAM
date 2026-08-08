













































































































































































































      SubrOutine projec_IFC(HESS, PHESS, Tmp, PMAT, GMATRX_N, 
     &                      GMATRX_M, GRD, GRDTMP, NXM6)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Logical Constrain_opt

C Built the projected Hessian and Gradients for redundent internal
C optimizations (see JCP, 117, 9160, 2002). 
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



C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)







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




C
      Dimension HESS(NXM6, NXM6), PMAT(NXM6, NXM6), PHESS(NXM6, NXM6),
     &          GMATRX_M(NXM6, NXM6), GMATRX_N(NXM6, NXM6),
     &          GRD(NXM6), Scr(Maxredunco), Tmp(NXM6, NXm6),
     &          GRDTM(NXM6)
 
      LENGMAT=NXM6*NXM6
      CALL GETREC(20,'JOBARC','GI-MATRX',LENGMAT*IINTFP,GMATRX_M)
      CALL GETREC(20,'JOBARC','G-MATRX ',LENGMAT*IINTFP,GMATRX_N)

      Call XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, GMATRX_M, 
     &            NXM6, GMATRX_N, NXM6, 0.0D0, PMAT, NXM6)
C
      Write(6,*)
      Write(6,"(a)") "The Hessian projector:P"
      CALL OUTPUT(PMAT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
C
C Constrained redundent Internal optimizations; Built the C matrix
C as described in JCC, 17, 49, 1996. This applies to constrained 
C redundent internal optimiztions. Use Phes to store the C matrix,
C 
      Call Getrec(0, "JOBARC", "CONSTRNS", Ilength, Ijunk)
      Constrain_opt = .FALSE.
C
      Write(6,"(a,1x,I4,1x,I4)") "Constrained opt?", Ilength,
     &                            Iflags2(169)
C
      If (Ilength .Gt. 0) Then
          Constrain_opt = .TRUE.
          If (Iflags2(169) .NE. 3) Then
             Write(6,"(3x,2a)") "Inconsistency: There are constrains in" 
     &             ," input, but opt_control=constrained is not set."
              Call Errex
           Endif
      Endif
C 
      If (Constrain_opt) Then

         Call Build_cmatrix(NXM6, Phess)
C
C Built the CPC; Not that Gmatrx_N{I} are used as scratch matrices.
C     
         CALL XGEMM("T", "N", NXM6, NXM6, NXM6, 1.0D0, Phess,
     &               NXM6, Pmat, NXM6, 0.0D0, GMATRX_N, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, GMATRX_N,
     &               NXM6, Phess, NXM6, 0.0D0, GMATRX_M, NXM6)
      Write(6,*)
      Write(6,"(a)") "The CPC"
      CALL OUTPUT(GMATRX_M, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
C
C Invert the CPC 
C
CSSS         Call Dgefa(GMATRX_M, NXM6, NXM6, Ipvt, Info)
C
C#ifdef _NOSKIP
CSSS      CALL MINV(GMATRX_M, NXM6, NXM6, Tmp, Det, 1.0D-8, 0, 1)
      
      CALL EIG(GMATRX_M, Tmp, NXM6, NXM6, 1)
C
      Write(6,*)
      Write(6,"(a)") "The eigen vectors of the CPC"
      CALL OUTPUT(Tmp, 1, nxm6, 1, nxM6, nxm6, nxm6, 1)
      Write(6,"(a)") "The eigenvalues of the CPC "
      Write(6, "(4F10.5)") (Gmatrx_M(I,I), I=1, Nxm6)

      DO I = 1, NXM6
         IF (GMATRX_M(I, I) .LE. 1.0D-10) THEN
             GMATRX_M(I, I) = 0.0D0
         ELSE
             GMATRX_M(I, I) = 1.0D0/GMATRX_M(I, I)
         ENDIF
      ENDDO
C
C Built the generalized inverse of CPC matrix,
C
      CALL XGEMM('N', 'N', NXM6, NXM6, NXM6, 1.0D0, Tmp,
     &            NXM6, GMATRX_M, NXM6, 0.0D0, GMATRX_N, NXM6)
      CALL XGEMM('N', 'T', NXM6, NXM6, NXM6, 1.0D0, GMATRX_N,
     &            NXM6, Tmp, NXM6, 0.0D0, GMATRX_M, NXM6)
C#endif 
C
      Write(6,*) "The inverse of CPC"
      CALL OUTPUT(GMATRX_M, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
C
C Create CP and PC
C
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Phess,
     &               NXM6, Pmat, NXM6, 0.0D0, GMATRX_N, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Pmat,
     &               NXM6, Phess, NXM6, 0.0D0, Tmp, NXM6)
C
C Built PC(CPC)^(-1)CP
C 
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Tmp,
     &               NXM6, GMATRX_M, NXM6, 0.0D0, Phess, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, Phess,
     &               NXM6, GMATRX_N, NXM6, 0.0D0, GMATRX_M, NXM6)

      Write(6,*)
      Write(6,"(a)") "The PC(CPC)^(-1)CP"
      CALL OUTPUT(GMATRX_M, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
      CALL OUTPUT(PMAT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
C
C Built P = P - PC(CPC)^(-1)CP

       Call Daxpy(NXM6*NXM6, -1.0D0, GMATRX_M, 1, Pmat, 1)

      Write(6,*)
      Write(6,"(a)") "The constrained Hessian projector:P"
      CALL OUTPUT(PMAT, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
C
C Built PHP (projected Hessian)
C
     
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)

         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)

C Also project the gradients ie. form PG. Then copy the projected gradient
C to the FI array so that the rest of the code can proceed.
C
      Write(6,"(a)") "The gradients, g"
      Write(6,"(6(1x,F15.7))") (Grd(I), I=1,NXM6)
         CALL XGEMM('N', 'N', NXM6, 1, NXM6, 1.0D0, PMAT,
     &               NXM6, GRD, NXM6, 0.0D0, GRDTMP, NXM6)
         CALL DCOPY(NXM6, GRDTMP, 1, GRD, 1)
C
C
      Write(6,*)
      Write(6,"(a)") "The projected Hessian, (PHP)"
      CALL OUTPUT(HESS, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
      Write(6,"(a)") "The projected gradients, Hg"
      Write(6,"(6(1x,F15.7))") (Grd(I), I=1,NXM6)
      
      Else 
C
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PMAT,
     &            NXM6, HESS, NXM6, 0.0D0, PHESS, NXM6)
         CALL XGEMM("N", "N", NXM6, NXM6, NXM6, 1.0D0, PHESS,
     &            NXM6, PMAT, NXM6, 0.0D0, HESS, NXM6)
C
C Also project the gradients. Then copy the projected gradient
C to the FI array so that the rest of the code can proceed.
C
         CALL XGEMM('N', 'N', NXM6, 1, NXM6, 1.0D0, PMAT,
     &               NXM6, GRD, NXM6, 0.0D0, GRDTMP, NXM6)
         CALL DCOPY(NXM6, GRDTMP, 1, GRD, 1)

C
          Write(6,"(a)")"The projected hessian (PHP)"
          CALL OUTPUT(HESS, 1, NXM6, 1, NXM6, NXM6, NXM6, 1)
          Write(6,*) 
          Write(6,"(a)") "The projected gradients Pg"
          Write(6,"(6(1x,F15.7))") (Grd(I), I=1,NXM6)
      Endif 
           
C------The lines below are for Debugging only--------------

C------------------------------------------------------
      RETURN
      END
     
