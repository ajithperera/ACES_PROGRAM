










      SUBROUTINE SDINT(NHKTA,NHKTB,KHKTA,KHKTB,ICENTA,LDIAG,
     &                 ISTEPA,ISTEPB,ISTEPU,ISTEPV,NAHGTF,
     &                 NATOMC,CHARGE,WORK1)
C
C   THIS ROUTINE CALCULATES INTEGRALS REQUIRED FOR THE
C   SPIN-DIPOLAR CONTRIBUTION TO NMR COUPLING CONSTANTS
C
CEND
C
C   JG 4/93
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LDIAG
      DIMENSION WORK1(1),CHARGE(1)

c Dimensions of the WORK2 and WORK3 arrays of the VDINT AME. For
c machines that have small stack space, these values must be small.
c Marshall Cory, who has access to CRAY machines with large memory, uses
c 900000 and 500000 is recommended for the others. Small values for
c these arrays will crash the program for basis sets that have large
c primitive spaces. For example the ANO basis sets. A. Perera, 03/2005. 

      PARAMETER (LWORK2 = 300 0000)
      PARAMETER (LWORK3 = 900 0000)

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
     &           *MXCONT*MXCONT)
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



      COMMON /ADER/ ADER0 (MXAQNS)
      COMMON /POINTER/IS0000, IS000X, IS000Y, IS000Z,
     &                IS00XX, IS00XY, IS00XZ, IS00YY,
     &                IS00YZ, IS00ZZ, IT0000, IT000X,
     &                IT000Y, IT000Z, IT00XX, IT00XY,
     &                IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     &                ID0000, ID000X, ID000Y, ID000Z,
     &                ID00XX, ID00XY, ID00XZ, ID00YX,
     &                ID00YY, ID00YZ, ID00ZX, ID00ZY,
     &                ID00ZZ,
     &                IA0000, IA0X00, IA0Y00, IA0Z00,
     &                IAXX00, IAXY00, IAXZ00, IAYY00,
     &                IAYZ00, IAZZ00, IA000X, IA000Y,
     &                IA000Z, IA00XX, IA00XY, IA00XZ,
     &                IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     &                IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     &                IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON/SCRPOIN/ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,
     &               ISCR8,ISCR9,ISCR10,ISCR11,ISCR12,ISCR13,ISCR14,
     &               ISCR15,ISCR16,ISCR17,ISCR18,ISCR19,ISCR20,
     &               JSCR1,JSCR2,JSCR3,JSCR4,JSCR5,JSCR6,
     &               JSCR7,JSCR8,JSCR9,JSCR10,JSCR11,JSCR12,JSCR13,
     &               JSCR14,JSCR15,JSCR16,JSCR17,JSCR18,JSCR19,
     &               JSCR20,JSCR21,JSCR22,JSCR23,JSCR24,JSCR25,JSCR26,
     &               JSCR27,JSCR28 
      COMMON /CENTC/ SIGNCX(MXCENT), SIGNCY(MXCENT), SIGNCZ(MXCENT),
     &               NCENTC(MXCENT), JSYMC(MXCENT),  JCENTC(MXCENT),
     &               ICXVEC(MXCENT), ICYVEC(MXCENT), ICZVEC(MXCENT)
      COMMON /LMNS/ LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN),
     &              LVALUB(MXAQN), MVALUB(MXAQN), NVALUB(MXAQN)
C
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
C
      DATA ONE,ONEM,TWO,THREE/1.D0,-1.D0,2.D0,3.D0/
      DATA THRSH /1.D-20/
C
      INT20=-NATOMC
      DO 100 I=1,KHKTA*KHKTB
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
       MAXT = LVALA + LVALB
       MAXU = MVALA + MVALB
       MAXV = NVALA + NVALB
       INT20 = INT20 + NATOMC
       IADRAV = 1
       DO 200 IV = 0,MAXV
        EV = ODC00Z(ISTREV + IV)
        IADRAU = IADRAV
        DO 300 IU = 0,MAXU
         EE = ODC00Y(ISTREU + IU)*EV
         DO 400 IT = 0,MAXT
          EEE = ODC00X(ISTRET + IT)*EE
          IADR00 = IADRAU + IT
          IADR0T = IADR00 + 1
          IADR0U = IADR00 + ISTEPU
          IADR0V = IADR00 + ISTEPV
          IADRTT = IADR0T + 1
          IADRTU = IADR0T + ISTEPU
          IADRTV = IADR0T + ISTEPV
          IADRUU = IADR0U + ISTEPU
          IADRUV = IADR0U + ISTEPV
          IADRVV = IADR0V + ISTEPV
          IADD = - NAHGTF
C
C   LOOP OVER NUCLEI 
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
          DO 500 IATOM = 1,NATOMC
           IADD=IADD+NAHGTF
           FACT=ONEM/CHARGE(JCENTC(IATOM))
           INT2=INT20+IATOM
           AH00=AHGTF(IADR00+IADD)
           AH0T=AHGTF(IADR0T+IADD)
           AH0U=AHGTF(IADR0U+IADD)
           AH0V=AHGTF(IADR0V+IADD)
           AHTT=AHGTF(IADRTT+IADD)
           AHTU=AHGTF(IADRTU+IADD)
           AHTV=AHGTF(IADRTV+IADD)
           AHUU=AHGTF(IADRUU+IADD)
           AHUV=AHGTF(IADRUV+IADD)
           AHVV=AHGTF(IADRVV+IADD)
C
           WORK1(JSCR5+INT2)=WORK1(JSCR5+INT2)
     &                      +FACT*EEE*(AHUU+AHVV-TWO*AHTT)
           WORK1(JSCR6+INT2)=WORK1(JSCR6+INT2)
     &                      -FACT*EEE*AHTU*THREE
           WORK1(JSCR7+INT2)=WORK1(JSCR7+INT2)
     &                      -FACT*EEE*AHTV*THREE
           WORK1(JSCR8+INT2)=WORK1(JSCR8+INT2)
     &                      +FACT*EEE*(AHTT+AHVV-TWO*AHUU)
           WORK1(JSCR9+INT2)=WORK1(JSCR9+INT2)
     &                      -FACT*EEE*AHUV*THREE
           WORK1(JSCR10+INT2)=WORK1(JSCR10+INT2)
     &                       +FACT*EEE*(AHTT+AHUU-TWO*AHVV)
500       CONTINUE
400      CONTINUE
         IADRAU = IADRAU + ISTEPU
300     CONTINUE
        IADRAV = IADRAV + ISTEPV
200    CONTINUE
100   CONTINUE
C
C  NOW FORM THE CONTRACTIONS
C
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
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR5+IATOM),NATOMC,
     &             WORK1(IS00XX+IOFF+IATOM),ISKIP)
c      write(*,*) '1',(work1(is00xx+ioff+iatom-1+j),
c     &         j=1,khkta*khktb,iskip)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR6+IATOM),NATOMC,
     &             WORK1(IS00XY+IOFF+IATOM),ISKIP)
c      write(*,*) '2',(work1(is00xy+ioff+iatom-1+j),
c     &         j=1,khkta*khktb,iskip)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR7+IATOM),NATOMC,
     &             WORK1(IS00XZ+IOFF+IATOM),ISKIP)
c      write(*,*) '3',(work1(is00xz+ioff+iatom-1+j),
c     &         j=1,khkta*khktb,iskip)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR8+IATOM),NATOMC,
     &             WORK1(IS00YY+IOFF+IATOM),ISKIP)
c      write(*,*) '4',(work1(is00yy+ioff+iatom-1+j),
c     &         j=1,khkta*khktb,iskip)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR9+IATOM),NATOMC,
     &             WORK1(IS00YZ+IOFF+IATOM),ISKIP)
c      write(*,*) '5',(work1(is00yz+ioff+iatom-1+j),
c     &         j=1,khkta*khktb,iskip)
         CALL SAXPY(KHKTA*KHKTB,CONT,WORK1(JSCR10+IATOM),NATOMC,
     &             WORK1(IS00ZZ+IOFF+IATOM),ISKIP)
c      write(*,*) '6',(work1(is00zz+ioff+iatom-1+j),
c     &         j=1,khkta*khktb,iskip)
710     CONTINUE
        ENDIF
       IOFF=IOFF+NATOMC
700   CONTINUE
      RETURN
      END
