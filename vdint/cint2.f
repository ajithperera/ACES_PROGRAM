










      SUBROUTINE CINT2(NHKTA,NHKTB,KHKTA,KHKTB,ISTEPA,ISTEPB,
     *                 ISTEPU,ISTEPV,NAHGTF,NATOMC,SECOND,LDIAG,
     *                 WORK1)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (DP5 = 0.5 D00)
      LOGICAL SECOND
      LOGICAL LDIAG
      DIMENSION WORK1(1)

c Dimensions of the WORK2 and WORK3 arrays of the VDINT AME. For
c machines that have small stack space, these values must be small.
c Marshall Cory, who has access to CRAY machines with large memory, uses
c 900000 and 500000 is recommended for the others. Small values for
c these arrays will crash the program for basis sets that have large
c primitive spaces. For example the ANO basis sets. A. Perera, 03/2005. 

      PARAMETER (LWORK2 = 300 0000)
      PARAMETER (LWORK3 = 900 0000)

c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)

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

      PARAMETER (LWKRLM = LWORK3 - 8020)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN
     *           *MXCONT*MXCONT)
      COMMON /CWORK3/ WK3LOW, SHGTF, AHGTF((MXCENT+1)*(2*MXQN+1)**3)
CSSS                      RLMCOF(LWKRLM),WK3HGH

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



      COMMON /ADER/ ADER0 (MXAQNS)
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
      COMMON/SCRPOIN/ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,
     *               ISCR9,ISCR10,ISCR11,ISCR12,ISCR13,ISCR14,ISCR15,
     *               ISCR16,ISCR17,ISCR18,ISCR19,ISCR20,
     *               JSCR1,JSCR2,JSCR3,JSCR4,JSCR5,JSCR6,JSCR7,JSCR8,
     *               JSCR9,JSCR10,JSCR11,JSCR12,JSCR13,JSCR14,JSCR15,
     *               JSCR16,JSCR17,JSCR18,JSCR19,JSCR20,JSCR21,JSCR22,
     *               JSCR23,JSCR24,JSCR25,JSCR26,JSCR27,JSCR28 

      COMMON /CENTC/ SIGNCX(MXCENT), SIGNCY(MXCENT), SIGNCZ(MXCENT),
     *               NCENTC(MXCENT), JSYMC(MXCENT),  JCENTC(MXCENT),
     *               ICXVEC(MXCENT), ICYVEC(MXCENT), ICZVEC(MXCENT)
      COMMON /LMNS/ LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN),
     *              LVALUB(MXAQN), MVALUB(MXAQN), NVALUB(MXAQN)
      COMMON/GENCON/NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA THRSH /1.D-20/
C
      IF(.NOT.SECOND) THEN
C
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
CSSS       Write(6,"(3(1x,I4))") LVALA, MVALA, NVALA
CSSS       Write(6,"(3(1x,I4))") LVALB, MVALB, NVALB
C
C    **********************************************************
C    ***** CALCULATE OVERLAP AND KINETIC ENERGY INTEGRALS *****
C    **********************************************************
C
       DERX0 = SHGTF*ODC00X(ISTRET)
       DERY0 = SHGTF*ODC00Y(ISTREU)
       DERZ0 = SHGTF*ODC00Z(ISTREV)
       DERX1 = SHGTF*ODC10X(ISTRET)
       DERY1 = SHGTF*ODC10Y(ISTREU)
       DERZ1 = SHGTF*ODC10Z(ISTREV)
       DERX2 = SHGTF*ODC20X(ISTRET)
       DERY2 = SHGTF*ODC20Y(ISTREU)
       DERZ2 = SHGTF*ODC20Z(ISTREV)
       DERX3 = SHGTF*ODC30X(ISTRET)
       DERY3 = SHGTF*ODC30Y(ISTREU)
       DERZ3 = SHGTF*ODC30Z(ISTREV)
C
       WORK1(ISCR1+I) =  DERX0*DERY0*DERZ0
       WORK1(ISCR2+I) =  DERX1*DERY0*DERZ0
       WORK1(ISCR3+I)=  DERX0*DERY1*DERZ0
       WORK1(ISCR4+I) =  DERX0*DERY0*DERZ1

       WORK1(ISCR5+I)  =  - (DERX2*DERY0*DERZ0
     *      + DERX0*DERY2*DERZ0 + DERX0*DERY0*DERZ2)*DP5
       WORK1(ISCR6+I)  =  - (DERX3*DERY0*DERZ0
     *      + DERX1*DERY2*DERZ0 + DERX1*DERY0*DERZ2)*DP5
       WORK1(ISCR7+I)  =  - (DERX2*DERY1*DERZ0
     *      + DERX0*DERY3*DERZ0 + DERX0*DERY1*DERZ2)*DP5
       WORK1(ISCR8+I)  =  - (DERX2*DERY0*DERZ1
     *      + DERX0*DERY2*DERZ1 + DERX0*DERY0*DERZ3)*DP5
100   CONTINUE
CSSS      Write(6,*) " Primitive Overlap derivatie in CINT2"
CSSS      Write(6,"(4(1x,f20.10))") (WORK1(IS000Z+IOFF), ioff=1,KHKTA*KHKTB)
CSSS      Write(6,*)
C
C  NOW FORM THE CONTRACTIONS
C
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      IOFF=0
      DO 150 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 150 IRCB=1,MAXB
        IOFF=IOFF+1
        CONT=CONTA(IRCA)*CONTB(IRCB)
        IF(ABS(CONT).GT.THRSH) THEN
         CALL SAXPY(KHKTA*KHKTB*8,CONT,WORK1(ISCR1+1),1,
     *              WORK1(IS0000+IOFF),ISKIP)
        ENDIF
150   CONTINUE
C
CSSS      Write(6,*) "Contracted Overlap derivatie in CINT2"
CSSS      Write(6,"(4(1x,f20.10))") (WORK1(IS000Z+IOFF), ioff=1,KHKTA*KHKTB)
CSSS      Write(6,*)
      ELSE
C
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
C    **********************************************************
C    ***** CALCULATE OVERLAP AND KINETIC ENERGY INTEGRALS *****
C    **********************************************************
C
       DERX0 = SHGTF*ODC00X(ISTRET)
       DERY0 = SHGTF*ODC00Y(ISTREU)
       DERZ0 = SHGTF*ODC00Z(ISTREV)
       DERX1 = SHGTF*ODC10X(ISTRET)
       DERY1 = SHGTF*ODC10Y(ISTREU)
       DERZ1 = SHGTF*ODC10Z(ISTREV)
       DERX2 = SHGTF*ODC20X(ISTRET)
       DERY2 = SHGTF*ODC20Y(ISTREU)
       DERZ2 = SHGTF*ODC20Z(ISTREV)
       DERX3 = SHGTF*ODC30X(ISTRET)
       DERY3 = SHGTF*ODC30Y(ISTREU)
       DERZ3 = SHGTF*ODC30Z(ISTREV)
       DERX4 = SHGTF*ODC40X(ISTRET)
       DERY4 = SHGTF*ODC40Y(ISTREU)
       DERZ4 = SHGTF*ODC40Z(ISTREV)
C
       WORK1(ISCR1+I) =  DERX0*DERY0*DERZ0
       WORK1(ISCR2+I) =  DERX1*DERY0*DERZ0
       WORK1(ISCR3+I) =  DERX0*DERY1*DERZ0
       WORK1(ISCR4+I) =  DERX0*DERY0*DERZ1
       WORK1(ISCR5+I)  =  - (DERX2*DERY0*DERZ0
     *      + DERX0*DERY2*DERZ0 + DERX0*DERY0*DERZ2)*DP5
       WORK1(ISCR6+I)  =  - (DERX3*DERY0*DERZ0
     *      + DERX1*DERY2*DERZ0 + DERX1*DERY0*DERZ2)*DP5
       WORK1(ISCR7+I)  =  - (DERX2*DERY1*DERZ0
     *      + DERX0*DERY3*DERZ0 + DERX0*DERY1*DERZ2)*DP5
       WORK1(ISCR8+I)  =  - (DERX2*DERY0*DERZ1
     *      + DERX0*DERY2*DERZ1 + DERX0*DERY0*DERZ3)*DP5
       WORK1(ISCR9+I) = + DERX2*DERY0*DERZ0
       WORK1(ISCR10+I) = + DERX1*DERY1*DERZ0
       WORK1(ISCR11+I) = + DERX1*DERY0*DERZ1
       WORK1(ISCR12+I) = + DERX0*DERY2*DERZ0
       WORK1(ISCR13+I) = + DERX0*DERY1*DERZ1
       WORK1(ISCR14+I) = + DERX0*DERY0*DERZ2
       WORK1(ISCR15+I) = - (DERX4*DERY0*DERZ0
     *         + DERX2*DERY2*DERZ0 + DERX2*DERY0*DERZ2)*DP5
       WORK1(ISCR16+I) = - (DERX3*DERY1*DERZ0
     *         + DERX1*DERY3*DERZ0 + DERX1*DERY1*DERZ2)*DP5
       WORK1(ISCR17+I) = - (DERX3*DERY0*DERZ1
     *         + DERX1*DERY2*DERZ1 + DERX1*DERY0*DERZ3)*DP5
       WORK1(ISCR18+I) = - (DERX2*DERY2*DERZ0
     *         + DERX0*DERY4*DERZ0 + DERX0*DERY2*DERZ2)*DP5
       WORK1(ISCR19+I) = - (DERX2*DERY1*DERZ1
     *         + DERX0*DERY3*DERZ1 + DERX0*DERY1*DERZ3)*DP5
       WORK1(ISCR20+I) = - (DERX2*DERY0*DERZ2
     *         + DERX0*DERY2*DERZ2 + DERX0*DERY0*DERZ4)*DP5
C
1100   CONTINUE
C
C  NOW FORM THE CONTRACTIONS
C
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      IOFF=0
      DO 1150 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 1150 IRCB=1,MAXB
        IOFF=IOFF+1
        CONT=CONTA(IRCA)*CONTB(IRCB)
        IF(ABS(CONT).GT.THRSH) THEN
         CALL SAXPY(KHKTA*KHKTB*8,CONT,WORK1(ISCR1+1),1,
     *              WORK1(IS0000+IOFF),ISKIP)
         CALL SAXPY(KHKTA*KHKTB*12,CONT,WORK1(ISCR9+1),1,
     *              WORK1(IS00XX+IOFF),ISKIP)
        ENDIF
1150   CONTINUE

CSSS      Write(6,*) "Overlap derivative"
CSSS      Write(6,"(4(1x,f20.10))") (WORK1(IS00XX+IOFF), ioff=0,nrca*nrcb)
CSSS
      ENDIF
C
C     **************************************************
C     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C     **************************************************
C
         MAXADD=1
         IF(SECOND) MAXADD=2
         INT20 = - NATOMC
         DO 250 I=1,KHKTA*KHKTB
          ICOMPA=(I-1)/KHKTB+1
          ICOMPB=I-(ICOMPA-1)*KHKTB
          LVALA=LVALUA(ICOMPA)
          MVALA=MVALUA(ICOMPA)
          NVALA=NVALUA(ICOMPA)
          ISTRAT=1+ISTEPA*LVALA
          ISTRAU=1+ISTEPA*MVALA
          ISTRAV=1+ISTEPA*NVALA
          LVALB=LVALUB(ICOMPB)
          MVALB=MVALUB(ICOMPB) 
          NVALB=NVALUB(ICOMPB)
          ISTRET=ISTRAT+ISTEPB*LVALB
          ISTREU=ISTRAU+ISTEPB*MVALB
          ISTREV=ISTRAV+ISTEPB*NVALB
          MAXT = LVALA + LVALB + MAXADD
          MAXU = MVALA + MVALB + MAXADD
          MAXV = NVALA + NVALB + MAXADD
          INT20 = INT20 + NATOMC
          IADRAV = 1
          DO 200 IV = 0,MAXV
           IADREV = ISTREV + IV
           EV = ODC00Z(IADREV)
           FV = ODC10Z(IADREV)
           GV = ODC20Z(IADREV)
           IADRAU = IADRAV
           DO 300 IU = 0,MAXU
            IADREU = ISTREU + IU
            EU = ODC00Y(IADREU)
            FU = ODC10Y(IADREU)
            GU = ODC20Y(IADREU)
            EE = EU*EV
            FE = FU*EV
            GE = GU*EV
            EF = EU*FV
            FF = FU*FV
            EG = EU*GV
            DO 400 IT = 0,MAXT
             IADRET = ISTRET + IT
             ET = ODC00X(IADRET)
             FT = ODC10X(IADRET)
             EEE = ET*EE
             FEE = FT*EE
             EFE = ET*FE
             EEF = ET*EF
             IADR00 = IADRAU + IT
             IADR0T = IADR00 + 1
             IADR0U = IADR00 + ISTEPU
             IADR0V = IADR00 + ISTEPV
             IF (SECOND) THEN
              GT = ODC20X(IADRET)
              FFE = FT*FE
              FEF = FT*EF
              EFF = ET*FF
              GEE = GT*EE
              EGE = ET*GE
              EEG = ET*EG
              IADRTT = IADR0T + 1
              IADRTU = IADR0T + ISTEPU
              IADRTV = IADR0T + ISTEPV
              IADRUU = IADR0U + ISTEPU
              IADRUV = IADR0U + ISTEPV
              IADRVV = IADR0V + ISTEPV
             END IF
             IADD = - NAHGTF
C
C                 ***** LOOP OVER NUCLEI *****
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                  DO 500 IATOM = 1,NATOMC
C
C                    PICK UP HGTF INTEGRALS
C
                     IADD = IADD + NAHGTF
                     AH00 = AHGTF(IADR00 + IADD)
                     AH0T = AHGTF(IADR0T + IADD)
                     AH0U = AHGTF(IADR0U + IADD)
                     AH0V = AHGTF(IADR0V + IADD)
C
C                    MULTIPLY BY EXPANSION COEFFICIENTS
C                    AND ADD TO APPROPRIATE CGTF INTEGRAL
C
                     INT2 = INT20 + IATOM
C
C                    UNDIFFERENTIATED INTEGRAL:
C
                     WORK1(JSCR1+INT2)=WORK1(JSCR1+INT2)+EEE*AH00
c                     WORK1(IA0000+INT2) = WORK1(IA0000+INT2) + EEE*AH00
C
C                    A DIFFERENTIATED INTEGRALS:
C
                     WORK1(JSCR2+INT2)=WORK1(JSCR2+INT2)+FEE*AH00
                     WORK1(JSCR3+INT2)=WORK1(JSCR3+INT2)+EFE*AH00
                     WORK1(JSCR4+INT2)=WORK1(JSCR4+INT2)+EEF*AH00
 
c                     WORK1(IA0X00+INT2) = WORK1(IA0X00+INT2) + FEE*AH00
c                     WORK1(IA0Y00+INT2) = WORK1(IA0Y00+INT2) + EFE*AH00
c                     WORK1(IA0Z00+INT2) = WORK1(IA0Z00+INT2) + EEF*AH00
C
C                    C DIFFERENTIATED INTEGRALS:
C
                     WORK1(JSCR5+INT2)=WORK1(JSCR5+INT2)-EEE*AH0T
                     WORK1(JSCR6+INT2)=WORK1(JSCR6+INT2)-EEE*AH0U
                     WORK1(JSCR7+INT2)=WORK1(JSCR7+INT2)-EEE*AH0V
c                     WORK1(IA000X+INT2) = WORK1(IA000X+INT2) - EEE*AH0T
c                     WORK1(IA000Y+INT2) = WORK1(IA000Y+INT2) - EEE*AH0U
c                     WORK1(IA000Z+INT2) = WORK1(IA000Z+INT2) - EEE*AH0V
C
C                    SECOND DERIVATIVES
C
                     IF (SECOND) THEN
                        AHTT = AHGTF(IADRTT + IADD)
                        AHTU = AHGTF(IADRTU + IADD)
                        AHTV = AHGTF(IADRTV + IADD)
                        AHUU = AHGTF(IADRUU + IADD)
                        AHUV = AHGTF(IADRUV + IADD)
                        AHVV = AHGTF(IADRVV + IADD)
C
C                       A-A DIFFERENTIATED INTEGRALS:
C
c                        WORK1(IAXX00+INT2) = WORK1(IAXX00+INT2)+GEE*AH00
c                        WORK1(IAXY00+INT2) = WORK1(IAXY00+INT2)+FFE*AH00
c                        WORK1(IAXZ00+INT2) = WORK1(IAXZ00+INT2)+FEF*AH00
c                        WORK1(IAYY00+INT2) = WORK1(IAYY00+INT2)+EGE*AH00
c                        WORK1(IAYZ00+INT2) = WORK1(IAYZ00+INT2)+EFF*AH00
c                        WORK1(IAZZ00+INT2) = WORK1(IAZZ00+INT2)+EEG*AH00
                      WORK1(JSCR8+INT2) = WORK1(JSCR8+INT2)+GEE*AH00
                      WORK1(JSCR9+INT2) = WORK1(JSCR9+INT2)+FFE*AH00
                      WORK1(JSCR10+INT2) = WORK1(JSCR10+INT2)+FEF*AH00
                      WORK1(JSCR11+INT2) = WORK1(JSCR11+INT2)+EGE*AH00
                      WORK1(JSCR12+INT2) = WORK1(JSCR12+INT2)+EFF*AH00
                      WORK1(JSCR13+INT2) = WORK1(JSCR13+INT2)+EEG*AH00
C
C                       A-C DIFFERENTIATED INTEGRALS:
C
c                        WORK1(IA0X0X+INT2) = WORK1(IA0X0X+INT2)-FEE*AH0T
c                        WORK1(IA0X0Y+INT2) = WORK1(IA0X0Y+INT2)-FEE*AH0U
c                        WORK1(IA0X0Z+INT2) = WORK1(IA0X0Z+INT2)-FEE*AH0V
c                        WORK1(IA0Y0X+INT2) = WORK1(IA0Y0X+INT2)-EFE*AH0T
c                        WORK1(IA0Y0Y+INT2) = WORK1(IA0Y0Y+INT2)-EFE*AH0U 
c                        WORK1(IA0Y0Z+INT2) = WORK1(IA0Y0Z+INT2)-EFE*AH0V 
c                        WORK1(IA0Z0X+INT2) = WORK1(IA0Z0X+INT2)-EEF*AH0T 
c                        WORK1(IA0Z0Y+INT2) = WORK1(IA0Z0Y+INT2)-EEF*AH0U 
c                        WORK1(IA0Z0Z+INT2) = WORK1(IA0Z0Z+INT2)-EEF*AH0V 
                      WORK1(JSCR14+INT2) = WORK1(JSCR14+INT2)-FEE*AH0T
                      WORK1(JSCR15+INT2) = WORK1(JSCR15+INT2)-FEE*AH0U 
                      WORK1(JSCR16+INT2) = WORK1(JSCR16+INT2)-FEE*AH0V
                      WORK1(JSCR17+INT2) = WORK1(JSCR17+INT2)-EFE*AH0T
                      WORK1(JSCR18+INT2) = WORK1(JSCR18+INT2)-EFE*AH0U
                      WORK1(JSCR19+INT2) = WORK1(JSCR19+INT2)-EFE*AH0V
                      WORK1(JSCR20+INT2) = WORK1(JSCR20+INT2)-EEF*AH0T
                      WORK1(JSCR21+INT2) = WORK1(JSCR21+INT2)-EEF*AH0U
                      WORK1(JSCR22+INT2) = WORK1(JSCR22+INT2)-EEF*AH0V
C
C                       C-C DIFFERENTIATED INTEGRALS:
C
c                        WORK1(IA00XX+INT2) = WORK1(IA00XX+INT2)+EEE*AHTT
c                        WORK1(IA00XY+INT2) = WORK1(IA00XY+INT2)+EEE*AHTU
c                        WORK1(IA00XZ+INT2) = WORK1(IA00XZ+INT2)+EEE*AHTV
c                        WORK1(IA00YY+INT2) = WORK1(IA00YY+INT2)+EEE*AHUU
c                        WORK1(IA00YZ+INT2) = WORK1(IA00YZ+INT2)+EEE*AHUV
c                        WORK1(IA00ZZ+INT2) = WORK1(IA00ZZ+INT2)+EEE*AHVV
                      WORK1(JSCR23+INT2) = WORK1(JSCR23+INT2)+EEE*AHTT
                      WORK1(JSCR24+INT2) = WORK1(JSCR24+INT2)+EEE*AHTU
                      WORK1(JSCR25+INT2) = WORK1(JSCR25+INT2)+EEE*AHTV
                      WORK1(JSCR26+INT2) = WORK1(JSCR26+INT2)+EEE*AHUU
                      WORK1(JSCR27+INT2) = WORK1(JSCR27+INT2)+EEE*AHUV
                      WORK1(JSCR28+INT2) = WORK1(JSCR28+INT2)+EEE*AHVV
                    END IF
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
  250 CONTINUE
C
C  NOW FORM THE CONTRACTIONS
C
      INUMB=7
      IF(SECOND) INUMB=28
      ISKIP=NRCA*NRCB*NATOMC
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2*NATOMC
      IOFF=0
      DO 700 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 700 IRCB=1,MAXB
        CONT=CONTA(IRCA)*CONTB(IRCB)
        IF(ABS(CONT).GT.THRSH) THEN
         DO 710 IATOM=1,NATOMC
         CALL SAXPY(KHKTA*KHKTB*INUMB,CONT,WORK1(JSCR1+IATOM),NATOMC,
     *              WORK1(IA0000+IOFF+IATOM),ISKIP)
710     CONTINUE
        ENDIF
       IOFF=IOFF+NATOMC
700   CONTINUE
      RETURN
      END
