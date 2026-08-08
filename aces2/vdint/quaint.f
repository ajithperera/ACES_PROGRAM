










       SUBROUTINE QUAINT(NHKTA,NHKTB,KHKTA,KHKTB,LDIAG,ISTEPA,ISTEPB,
     *                  ONECEN,CORPX,CORPY,CORPZ,EXPPI,DIFDIP,WORK1)
C
C  CALCULATES QUADRUPOL INTEGRALS: <mu|x_i x_j|nu>
C
CEND
C
C JG 4/91
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)

c These parameters are gathered from vmol and vdint and are used by ecp
c as well. It just so happens that the vmol parameters do not exist in
c vdint and vice versa. LET'S TRY TO KEEP IT THAT WAY!

c VMOL PARAMETERS ------------------------------------------------------

C     MAXPRIM - Maximum number of primitives for a given shell.
      INTEGER    MAXPRIM
      PARAMETER (MAXPRIM=72)

C     MAXFNC  - Maximum number of contracted functions for a given shell.
C               (vmol/readin requires this to be the same as MAXPRIM)
      INTEGER    MAXFNC
      PARAMETER (MAXFNC=MAXPRIM)

C     NHT     - Maximum angular momentum
      INTEGER    NHT
      PARAMETER (NHT=7)

C     MAXATM  - Maximum number of atoms
      INTEGER    MAXATM
      PARAMETER (MAXATM=100)

C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers.
      INTEGER    MXTNPR
      PARAMETER (MXTNPR=MAXPRIM*MAXPRIM)

C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers.
      INTEGER    MXTNCC
      PARAMETER (MXTNCC=180000)

C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers.
      INTEGER    MXTNSH
      PARAMETER (MXTNSH=200)

C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions).
c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

c VDINT PARAMETERS -----------------------------------------------------

C     MXPRIM - Maximum number of primitives for all symmetry
C              inequivalent centers.
      INTEGER    MXPRIM
      PARAMETER (MXPRIM=MXTNPR)

C     MXSHEL - Maximum number of shells for all symmetry inequivalent centers.
      INTEGER    MXSHEL
      PARAMETER (MXSHEL=MXTNSH)

C     MXCORB - Maximum number of contracted basis functions.
      INTEGER    MXCORB
      PARAMETER (MXCORB=MXCBF)

C     MXORBT - Length of the upper or lower triangle length of MXCORB.
      INTEGER    MXORBT
      PARAMETER (MXORBT=MXCORB*(MXCORB+1)/2)

C     MXAOVC - Maximum number of subshells per center.
      INTEGER    MXAOVC,    MXAOSQ
      PARAMETER (MXAOVC=32, MXAOSQ=MXAOVC*MXAOVC)

c     MXCONT - ???
      INTEGER    MXCONT
      PARAMETER (MXCONT=MXAOVC)


c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN
     *           *MXCONT*MXCONT)
      LOGICAL LDIAG, ONECEN, DIFDIP
      DIMENSION WORK1(1)

c Dimensions of the WORK2 and WORK3 arrays of the VDINT AME. For
c machines that have small stack space, these values must be small.
c Marshall Cory, who has access to CRAY machines with large memory, uses
c 900000 and 500000 is recommended for the others. Small values for
c these arrays will crash the program for basis sets that have large
c primitive spaces. For example the ANO basis sets. A. Perera, 03/2005. 

      PARAMETER (LWORK2 = 300 0000)
      PARAMETER (LWORK3 = 900 0000)
      PARAMETER (LWKRLM = LWORK3 - 8020)
      COMMON /CWORK3/ WK3LOW, SHGTF, AHGTF((MXCENT+1)*(2*MXQN+1)**3)
CSSS                  RLMCOF(LWKRLM), WK3HGH

      COMMON/ODCS/ODC10X(mxqn*mxqn*(2*mxqn+4)),
     &            ODC20X(mxqn*mxqn*(2*mxqn+4)),
     &            ODC30X(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC40X(mxqn*mxqn*(2*mxqn+4)),
     &            ODC00X(mxqn*(mxqn+5)*(2*mxqn+4)),
     &            ODC10Y(mxqn*mxqn*(2*mxqn+4)),
     &            ODC20Y(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC30Y(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC40Y(mxqn*mxqn*(2*mxqn+4)),
     &            ODC00Y(mxqn*(mxqn+5)*(2*mxqn+4)),
     &            ODC10Z(mxqn*mxqn*(2*mxqn+4)),
     &            ODC20Z(mxqn*mxqn*(2*mxqn+4)),
     &            ODC30Z(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC40Z(mxqn*mxqn*(2*mxqn+4)),
     &            ODC00Z(mxqn*(mxqn+5)*(2*mxqn+4)) 



      COMMON /POINTER/IS0000, IS000X, IS000Y, IS000Z,
     *                IS00XX, IS00XY, IS00XZ, IS00YY,
     *                IS00YZ, IS00ZZ, IT0000, IT000X,
     *                IT000Y, IT000Z, IT00XX, IT00XY,
     *                IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     *                ID0000, ID000X, ID000Y, ID000Z,
     *                ID00XX, ID00XY, ID00XZ, ID00YX,
     *                ID00YY, ID00YZ, ID00ZX, ID00ZY,
     *                ID00ZZ,
     *                IA0000, IA0X00, IA0Y00, IA0Z00,
     *                IAXX00, IAXY00, IAXZ00, IAYY00,
     *                IAYZ00, IAZZ00, IA000X, IA000Y,
     *                IA000Z, IA00XX, IA00XY, IA00XZ,
     *                IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     *                IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     *                IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON/POINT2/IQ00XX,IQ00XY,IQ00XZ,IQ00YY,IQ00YZ,IQ00ZZ
      COMMON/SCRPOIN/ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,
     *               ISCR9,ISCR10,ISCR11,ISCR12,ISCR13,ISCR14,ISCR15,
     *               ISCR16,ISCR17,ISCR18,ISCR19,ISCR20,
     *               ISCR21,ISCR22,ISCR23,ISCR24,ISCR(24)
      COMMON /LMNS/ LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN),
     *              LVALUB(MXAQN), MVALUB(MXAQN), NVALUB(MXAQN)
      COMMON /GENCON/ NRCA, NRCB, CONTA(MXCONT), CONTB(MXCONT)
      DATA THRSH /1.D-20/,TWO/2.D0/,HALF/0.5D0/,HALFM/-0.5D0/
C
C  FIRST FORM THE QUADRUPOLE INTEGRALS OVER THE PRIMITIVES
C
      IF(.NOT.DIFDIP.OR.ONECEN) THEN
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
       DO 100 I = 1,KHKTA*KHKTB
       ICOMPA=(I-1)/KHKTB+1
       ICOMPB=I-(ICOMPA-1)*KHKTB
       LVALA = LVALUA(ICOMPA)
       MVALA = MVALUA(ICOMPA)
       NVALA = NVALUA(ICOMPA)
       ISTRAT = 1 + ISTEPA*LVALA
       ISTRAU = 1 + ISTEPA*MVALA
       ISTRAV = 1 + ISTEPA*NVALA
       LVALB = LVALUB(ICOMPB)
       MVALB = MVALUB(ICOMPB)
       NVALB = NVALUB(ICOMPB)
       ISTRET = ISTRAT + ISTEPB*LVALB
       ISTREU = ISTRAU + ISTEPB*MVALB
       ISTREV = ISTRAV + ISTEPB*NVALB
C
       SX0 = SHGTF*ODC00X(ISTRET)
       SY0 = SHGTF*ODC00Y(ISTREU)
       SZ0 = SHGTF*ODC00Z(ISTREV)
       DX0 = SHGTF*ODC00X(ISTRET + 1) + CORPX*SX0
       DY0 = SHGTF*ODC00Y(ISTREU + 1) + CORPY*SY0
       DZ0 = SHGTF*ODC00Z(ISTREV + 1) + CORPZ*SZ0
       QX0 = TWO*SHGTF*ODC00X(ISTRET + 2) + TWO*CORPX*DX0 
     &                                  -(CORPX**2+HALFM*EXPPI)*SX0
       QY0 = TWO*SHGTF*ODC00Y(ISTREU + 2) + TWO*CORPY*DY0 
     &                                  -(CORPY**2+HALFM*EXPPI)*SY0
       QZ0 = TWO*SHGTF*ODC00Z(ISTREV + 2) + TWO*CORPZ*DZ0
     &                                  -(CORPZ**2+HALFM*EXPPI)*SZ0
C
       WORK1(ISCR1+I) =  QX0*SY0*SZ0
       WORK1(ISCR2+I) =  DX0*DY0*SZ0
       WORK1(ISCR3+I) =  DX0*SY0*DZ0
       WORK1(ISCR4+I) =  SX0*QY0*SZ0
       WORK1(ISCR5+I) =  SX0*DY0*DZ0
       WORK1(ISCR6+I) =  SX0*SY0*QZ0
C
100   CONTINUE
C
C  NOW FORM THE CONTRACTIONS
C
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP= NRCA*(NRCA+1)/2
      IOFF=0
      DO 200 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 200 IRCB=1,MAXB
        IOFF=IOFF+1
        CONT=CONTA(IRCA)*CONTB(IRCB) 
        IF(ABS(CONT).GT.THRSH) THEN
         CALL SAXPY(KHKTA*KHKTB*6,CONT,WORK1(ISCR1+1),1,
     *              WORK1(IQ00XX+IOFF),ISKIP)
        ENDIF
200   CONTINUE
C
      ELSE
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
       DO 1100 I = 1,KHKTA*KHKTB
       ICOMPA=(I-1)/KHKTB+1
       ICOMPB=I-(ICOMPA-1)*KHKTB
       LVALA = LVALUA(ICOMPA)
       MVALA = MVALUA(ICOMPA)
       NVALA = NVALUA(ICOMPA)
       ISTRAT = 1 + ISTEPA*LVALA
       ISTRAU = 1 + ISTEPA*MVALA
       ISTRAV = 1 + ISTEPA*NVALA
       LVALB = LVALUB(ICOMPB)
       MVALB = MVALUB(ICOMPB)
       NVALB = NVALUB(ICOMPB)
       ISTRET = ISTRAT + ISTEPB*LVALB
       ISTREU = ISTRAU + ISTEPB*MVALB
       ISTREV = ISTRAV + ISTEPB*NVALB
C
       SX0 = SHGTF*ODC00X(ISTRET)
       SY0 = SHGTF*ODC00Y(ISTREU)
       SZ0 = SHGTF*ODC00Z(ISTREV)
       DX0 = SHGTF*ODC00X(ISTRET + 1) + CORPX*SX0
       DY0 = SHGTF*ODC00Y(ISTREU + 1) + CORPY*SY0
       DZ0 = SHGTF*ODC00Z(ISTREV + 1) + CORPZ*SZ0
       QX0 = SHGTF*ODC00X(ISTRET + 2)
       QY0 = SHGTF*ODC00Y(ISTREU + 2)
       QZ0 = SHGTF*ODC00Z(ISTREV + 2)
C
       WORK1(ISCR1+I) =  QX0*SY0*SZ0
       WORK1(ISCR2+I) =  DX0*DY0*SZ0
       WORK1(ISCR3+I) =  DX0*SY0*DZ0
       WORK1(ISCR4+I) =  SX0*QY0*SZ0
       WORK1(ISCR5+I) =  SX0*SY0*DZ0
       WORK1(ISCR6+I) =  SX0*SY0*QZ0
c       SX1 = SHGTF*ODC10X(ISTRET)
c       SY1 = SHGTF*ODC10Y(ISTREU)
c       SZ1 = SHGTF*ODC10Z(ISTREV)
c       DX1 = SHGTF*ODC10X(ISTRET + 1) + CORPX*SX1
c       DY1 = SHGTF*ODC10Y(ISTREU + 1) + CORPY*SY1
c       DZ1 = SHGTF*ODC10Z(ISTREV + 1) + CORPZ*SZ1
c       QX1 = SHGTF*ODC10X(ISTRET + 2)
c       QY1 = SHGTF*ODC10Y(ISTREU + 2)
c       QZ1 = SHGTF*ODC10Z(ISTREV + 2)
c       WORK1(ISCR7+I) = QX1*SY0*SZ0
c       WORK1(ISCR8+I) = DX1*DY0*SZ0
c       WORK1(ISCR9+I) = DX1*SY0*DZ0
c       WORK1(ISCR10+I) = SX1*QY0*SZ0
c       WORK1(ISCR11+I) = SX1*DY0*DZ0
c       WORK1(ISCR12+I) = SX1*SY0*QZ0
c       WORK1(ISCR13+I) = QX0*SY1*SZ0
c       WORK1(ISCR14+I) = DX0*DY1*SZ0
c       WORK1(ISCR15+I) = DX0*SY1*DZ0
c       WORK1(ISCR16+I) = SX0*QY1*SZ0
c       WORK1(ISCR17+I) = SX0*DY1*DZ0
c       WORK1(ISCR18+I) = SX0*SY1*QZ0
c       WORK1(ISCR19+I) = QX0*SY0*SZ1
c       WORK1(ISCR20+I) = DX0*DY0*SZ1
c       WORK1(ISCR21+I) = DX0*SY0*DZ1
c       WORK1(ISCR22+I) = SX0*QY0*SZ1
c       WORK1(ISCR23+I) = SX0*DY0*DZ1
c       WORK1(ISCR24+I) = SX0*SY0*QZ1
1100   CONTINUE
C
C  NOW FORM THE CONTRACTIONS
C
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP= NRCA*(NRCA+1)/2
      IOFF=0
      DO 1200 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 1200 IRCB=1,MAXB
        IOFF=IOFF+1
        CONT=CONTA(IRCA)*CONTB(IRCB) 
        IF(ABS(CONT).GT.THRSH) THEN
         CALL SAXPY(KHKTA*KHKTB*6,CONT,WORK1(ISCR1+1),1,
     *              WORK1(IQ00XX+IOFF),ISKIP)
c         CALL SAXPY(KHKTA*KHKTB*18,CONT,WORK1(ISCR7+1),1,
c     *              WORK1(IQ00XX+IOFF),ISKIP)
        ENDIF
1200   CONTINUE
      END IF
      RETURN
      END
