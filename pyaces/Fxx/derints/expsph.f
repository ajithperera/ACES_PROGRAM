










      SUBROUTINE EXPSPH(CMO,SCR)
C
C THIS ROUTINE EXPANDS A GIVEN SET OF MO COEFFICIENTS
C FROM SPHERICAL TO CARTESIAN GAUSSIAN BASIS FUNCTIONS.
C
CEND
C
C CODED  NOV/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL SHARE,DOPERT,DCORD,DCORGD,NOORBT
      INTEGER AND,OR,XOR
C
      DIMENSION CMO(1),SCR(1)
C

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
C
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
      COMMON/SYMMET/FMULT(0:7),PT(0:7),
     &              MAXLOP,MAXLOT,MULT(0:7),ISYTYP(3),
     &              ITYPE(8,36),NPARSU(8),NPAR(8),NAOS(8),
     &              NPARNU(8,8),IPTSYM(MXCORB,0:7),
     &              IPTCNT(3*MXCENT,0:7),NCRREP(0:7),
     &              IPTCOR(MXCENT*3),NAXREP(0:7),IPTAX(3),
     &              IPTXYZ(3,0:7)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/BASSPH/NBAS5(8),NBASIS5,NBASSQ5,NBASTT5
C
      DATA TWO/2.D0/,THREE/3.D0/,FOUR/4.D0/
C
      IBTAND(I,J)=AND(I,J) 
      IBTOR(I,J)=OR(I,J)
      IBTXOR(I,J)=XOR(I,J)
C
      IND10=0
      IND20=0
C
C   LOOP OVER ALL IRREPS
C
      DO 100 IRREP=0,MAXLOP
C
       IND1=0
       IND2=0
C
       NMO=NBAS5(IRREP+1)
       NAO=NBAS(IRREP+1) 
C
C   LOOP OVER ALL SHELLS
C
       DO 200 ISHELL=1,KMAX
C
        KHKTA=KHKT(ISHELL)
        NHKTA=NHKT(ISHELL)
        NRCA=NRCO(ISHELL)
        MULA=ISTBAO(ISHELL)
C
        IF(NHKTA.LE.2) THEN
C 
C S- AND P-TYPE
C
         DO 300 NA=1,KHKTA
          IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,NA))).EQ.0) THEN
           DO 400 IRCA=1,NRCA
C
           IND1=IND1+1
           IND2=IND2+1     
C
           DO 500 IMO=1,NMO 
C
            SCR(IND10+IND1+(IMO-1)*NAO)=
     &      CMO(IND20+IND2+(IMO-1)*NMO)
C
500        CONTINUE
C
400       CONTINUE
C
         ENDIF
300     CONTINUE
       ELSE IF(NHKTA.EQ.3) THEN
C
C D-TYPE
C     
        J=IND1+1-NRCA
        K=IND2+1-NRCA
        IOFFX21=0
        IOFFY21=0
        IOFFZ21=0
        IOFFXY1=0
        IOFFXZ1=0
        IOFFYZ1=0
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,1))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFX21=J
         IOFFX22=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,2))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFXY1=J
         IOFFXY2=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,3))).EQ.0) THEN 
         J=J+NRCA
         K=K+NRCA
         IOFFXZ1=J
         IOFFXZ2=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,4))).EQ.0) THEN 
         J=J+NRCA
         K=K+NRCA
         IOFFY21=J
         IOFFY22=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,5))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFYZ1=J
         IOFFYZ2=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,6))).EQ.0) THEN
         J=J+NRCA
         IOFFZ21=J
        ENDIF
        IND1=J+NRCA-1
        IND2=K+NRCA-1
        IF(IOFFX21.NE.0) THEN
         DO 401 IRCA=1,NRCA
          DO 501 IMO=1,NMO
           SCR(IND10+IOFFX21+IRCA-1+(IMO-1)*NAO)=
     &           -CMO(IND20+IOFFX22+IRCA-1+(IMO-1)*NMO)+
     &            CMO(IND20+IOFFY22+IRCA-1+(IMO-1)*NMO)
           SCR(IND10+IOFFY21+IRCA-1+(IMO-1)*NAO)=
     &          -CMO(IND20+IOFFX22+IRCA-1+(IMO-1)*NMO)-
     &           CMO(IND20+IOFFY22+IRCA-1+(IMO-1)*NMO)
            SCR(IND10+IOFFZ21+IRCA-1+(IMO-1)*NAO)=
     &      CMO(IND20+IOFFX22+IRCA-1+(IMO-1)*NMO)*TWO
501       CONTINUE
401      CONTINUE
        ENDIF
        IF(IOFFXY1.NE.0) THEN
         DO 402 IRCA=1,NRCA
          DO 502 IMO=1,NMO
           SCR(IND10+IOFFXY1+IRCA-1+(IMO-1)*NAO)=
     &      CMO(IND20+IOFFXY2+IRCA-1+(IMO-1)*NMO)
502       CONTINUE
402      CONTINUE
        ENDIF
        IF(IOFFXZ1.NE.0) THEN
         DO 403 IRCA=1,NRCA
          DO 503 IMO=1,NMO
           SCR(IND10+IOFFXZ1+IRCA-1+(IMO-1)*NAO)=
     &      CMO(IND20+IOFFXZ2+IRCA-1+(IMO-1)*NMO)
503       CONTINUE
403      CONTINUE
        ENDIF
        IF(IOFFYZ1.NE.0) THEN
         DO 404 IRCA=1,NRCA
          DO 504 IMO=1,NMO
           SCR(IND10+IOFFYZ1+IRCA-1+(IMO-1)*NAO)=
     &      CMO(IND20+IOFFYZ2+IRCA-1+(IMO-1)*NMO)
504       CONTINUE
404      CONTINUE
        ENDIF 
C
       ELSE IF(NHKTA.EQ.4) THEN
C
C F-TYPE
C
        J=IND1+1-NRCA
        K=IND2+1-NRCA
        IOFFX31=0
        IOFFY31=0
        IOFFZ31=0
        IOFFX2Y1=0
        IOFFX2Z1=0
        IOFFXY21=0
        IOFFXZ21=0
        IOFFYZ21=0
        IOFFY2Z1=0
        IOFFXYZ1=0
C
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,1))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFX31=J
         IOFFF1=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,2))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFX2Y1=J
         IOFFF2=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,3))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFX2Z1=J
         IOFFF3=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,4))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFXY21=J
         IOFFF4=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,5))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFXYZ1=J
         IOFFF5=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,6))).EQ.0) THEN
         J=J+NRCA
         IOFFXZ21=J
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,7))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFY31=J
         IOFFF7=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,8))).EQ.0) THEN
         J=J+NRCA
         K=K+NRCA
         IOFFY2Z1=J
         IOFFF8=K
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,9))).EQ.0) THEN
         J=J+NRCA
         IOFFYZ21=J
        ENDIF
        IF(IBTAND(MULA,IBTXOR(IRREP,ITYPE(NHKTA,10))).EQ.0) THEN
         J=J+NRCA
         IOFFZ31=J
        ENDIF
        IND1=J+NRCA-1
        IND2=K+NRCA-1
        IF(IOFFX31.NE.0) THEN
         DO 602 IRCA=1,NRCA
          DO 702 IMO=1,NMO
           SCR(IND10+IOFFX31+IRCA-1+(IMO-1)*NAO)=
     &      -CMO(IND20+IOFFF1+IRCA-1+(IMO-1)*NMO)
     &      +CMO(IND20+IOFFF4+IRCA-1+(IMO-1)*NMO)
           SCR(IND10+IOFFXY21+IRCA-1+(IMO-1)*NAO)=
     &      -CMO(IND20+IOFFF1+IRCA-1+(IMO-1)*NMO)
     &      -THREE*CMO(IND20+IOFFF4+IRCA-1+(IMO-1)*NMO)
           SCR(IND10+IOFFXZ21+IRCA-1+(IMO-1)*NAO)=
     &      FOUR*CMO(IND20+IOFFF1+IRCA-1+(IMO-1)*NMO)
702       CONTINUE
602      CONTINUE
        ENDIF
        IF(IOFFX2Y1.NE.0) THEN
         DO 603 IRCA=1,NRCA
          DO 703 IMO=1,NMO
           SCR(IND10+IOFFX2Y1+IRCA-1+(IMO-1)*NAO)=
     &      -CMO(IND20+IOFFF2+IRCA-1+(IMO-1)*NMO)
     &      +THREE*CMO(IND20+IOFFF7+IRCA-1+(IMO-1)*NMO)
           SCR(IND10+IOFFY31+IRCA-1+(IMO-1)*NAO)=
     &      -CMO(IND20+IOFFF2+IRCA-1+(IMO-1)*NMO)
     &      -CMO(IND20+IOFFF7+IRCA-1+(IMO-1)*NMO)
           SCR(IND10+IOFFYZ21+IRCA-1+(IMO-1)*NAO)=
     &      FOUR*CMO(IND20+IOFFF2+IRCA-1+(IMO-1)*NMO)
703       CONTINUE
603      CONTINUE
        ENDIF
        IF(IOFFX2Z1.NE.0) THEN
         DO 604 IRCA=1,NRCA
          DO 704 IMO=1,NMO
           SCR(IND10+IOFFX2Z1+IRCA-1+(IMO-1)*NAO)=
     &      -THREE*CMO(IND20+IOFFF3+IRCA-1+(IMO-1)*NMO)
     &      +CMO(IND20+IOFFF8+IRCA-1+(IMO-1)*NMO)
           SCR(IND10+IOFFY2Z1+IRCA-1+(IMO-1)*NAO)=
     &      -THREE*CMO(IND20+IOFFF3+IRCA-1+(IMO-1)*NMO)
     &      -CMO(IND20+IOFFF8+IRCA-1+(IMO-1)*NMO)
           SCR(IND10+IOFFZ31+IRCA-1+(IMO-1)*NAO)=
     &      TWO*CMO(IND20+IOFFF3+IRCA-1+(IMO-1)*NMO)
704       CONTINUE
604      CONTINUE
        ENDIF
        IF(IOFFXYZ1.NE.0) THEN
         DO 601 IRCA=1,NRCA
          DO 701 IMO=1,NMO
           SCR(IND10+IOFFXYZ1+IRCA-1+(IMO-1)*NAO)=
     &      CMO(IND20+IOFFF5+IRCA-1+(IMO-1)*NMO)
701       CONTINUE
601      CONTINUE
        ENDIF 
C
       ENDIF
200   CONTINUE
C
      IND10=IND10+NAO*NMO
      IND20=IND20+NMO*NMO
100   CONTINUE
C
c YAU : old
c     CALL ICOPY(IINTFP*IND10,SCR,1,CMO,1)
c YAU : new
      CALL DCOPY(IND10,SCR,1,CMO,1)
c YAU : end
C
      RETURN
      END
