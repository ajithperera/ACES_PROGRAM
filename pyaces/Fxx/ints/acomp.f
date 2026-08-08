










      SUBROUTINE ACOMP(SAAB,NUCA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

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

      parameter (nh4=4*nht-3, nh2=nht+nht+1)
      parameter (mxp2=maxprim*maxprim)
      COMMON/NULL_COM/  FACT(nh4),RFACT(nh4),FACTM(nh4),RFACTM(nh4)
     & ,MAA,IFD1,IFD2,KCD,KBCD,NHCD,NHBCD,NNC,IFPL(3),
     & NN1,NN2,NHKTA,NHKTB,NHKTC,NHKTD,KHKTA,KHKTB,KHKTC,KHKTD,NNB
      DIMENSION SAAB(NUCA,NUCA),AINT(NUC1,NUC1,NUC34),
     $ BINT(NUC34,NUC1,NUC1),AAINT(NUC1,NUC1,NUC34,2),
     $ BBINT(NUC34,NUC1,NUC1,2)
C...  TRIANGULARIZE THE RADIAL OVERLAP VECTOR IF IT IS SYMMETRIC !
C
      DO 20 I=2,NUCA
      DO 10 J=1,I-1
      SAAB(J,I) = SAAB(J,I)*2.
      SAAB(I,J) = 0.
   10 CONTINUE
   20 CONTINUE
      RETURN
      ENTRY BCOMP(SAAB,NUCA)
      DO 21 I=2,NUCA
      DO 11 J=1,I-1
      SAAB(J,I) = SAAB(J,I)*0.5
      SAAB(I,J) = SAAB(J,I)
   11 CONTINUE
   21 CONTINUE
      RETURN
      ENTRY AEXP(AINT,NUC1,NUC34)
C...  EXPAND THE PRIMITIVE INTEGRAL SET IF THE RADIAL OVERLAP
C     VECTOR HAD PREVIOUSLY BEEN TRIANGULARIZED.
C     THIS IS NECESSARY ONLY IN CASE OF GENERAL CONTRACTION !
C
C...  TRIANGULARITY IN P AND Q (INNER PAIR-LOOP)
      DO 50 IP=2,NUC1
      DO 40 IQ=1,IP-1
      DO 29 IRS=1,NUC34
   29 AINT(IQ,IP,IRS) = 0.5*(AINT(IQ,IP,IRS)+AINT(IP,IQ,IRS))
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 30 IRS=1,NUC34
   30 AINT(IP,IQ,IRS) = AINT(IQ,IP,IRS)
   40 CONTINUE
   50 CONTINUE
      RETURN
      ENTRY BEXP (BINT,NUC1,NUC34)
C     TRIANGULAR IN IR AND IS (OUTER PAIR-LOOP)
      DO 90 IR=2,NUC1
      DO 80 IS=1,IR-1
      DO 69 IPQ=1,NUC34
   69 BINT(IPQ,IR,IS) = 0.5*(BINT(IPQ,IS,IR)+BINT(IPQ,IR,IS))
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 70 IPQ=1,NUC34
   70 BINT(IPQ,IS,IR) = BINT(IPQ,IR,IS)
   80 CONTINUE
   90 CONTINUE
      RETURN
      ENTRY AAEXP (AAINT,NUC1,NUC34)
C...  TRIANGULARITY IN P AND Q (INNER PAIR-LOOP)
      DO 101 LOP=1,MAA
      DO 51 IP=2,NUC1
      DO 41 IQ=1,IP-1
      DO 26 IRS=1,NUC34
   26 AAINT(IQ,IP,IRS,LOP) = 0.5*(AAINT(IQ,IP,IRS,LOP)
     $ + AAINT(IP,IQ,IRS,LOP))
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 31 IRS=1,NUC34
   31 AAINT(IP,IQ,IRS,LOP) = AAINT(IQ,IP,IRS,LOP)
   41 CONTINUE
   51 CONTINUE
  101 CONTINUE
      RETURN
      ENTRY BBEXP (BBINT,NUC1,NUC34)
C     TRIANGULAR INIR AND IS (OUTER PAIR-LOOP)
      DO 111 LOP=1,MAA
      DO 91 IR=2,NUC1
      DO 81 IS=1,IR-1
      DO 65 IPQ=1,NUC34
   65 BBINT(IPQ,IR,IS,LOP) = 0.5*(BBINT(IPQ,IS,IR,LOP)
     $ + BBINT(IPQ,IR,IS,LOP))
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 71 IPQ=1,NUC34
   71 BBINT(IPQ,IS,IR,LOP) = BBINT(IPQ,IR,IS,LOP)
   81 CONTINUE
   91 CONTINUE
  111 CONTINUE
      RETURN
      END
