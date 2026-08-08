










      SUBROUTINE DMGAVR(KHKTA,KHKTB,IDENA,IDENB,
     &                  ISYMOP,NATOMC,
     &                  LDIAG,WORK1,MSZ,CSH,CSHSCF,
     &                  NCOORD)
C
C     CALCULATES THE EXPECTATION VALUE OF SECOND DERIVATIVE
C     ONE INETGRALS WITH RESPECT TO MAGNETIC PERTURBATION,
C     THAT MEANS THE CORRESPONDING CONTRIBUTIONS TO THE
C     DIAMAGNETIC SUSZEPTIBILITY AND CHEMICAL SHIELDING
C     TENSOR
C
CEND
C
C  UNIVERSITY OF KARLSRUHE,    AUGUST/91 JG
C
C  SYMMETRY INCLUDED FEBRUARY/92 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL ONECEN,LDIAG,SCF,NONHF
      DOUBLE PRECISION MSZ
C

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

      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN*
     &           MXCONT*MXCONT)
C
      DIMENSION WORK1(1),MSZ(3,3),CSH(3,NCOORD),
     &          CSHSCF(3,NCOORD)
C
      COMMON/DF/DSHELL(MXAQNS),DSHELL2(MXAQNS)
      COMMON/SYMMET/FMULT(0:7), PT(0:7),
     &              MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &              ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &              NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &              IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &              IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &              IPTXYZ(3,0:7)
      COMMON/CENTC/SIGNCX(MXCENT),SIGNCY(MXCENT),SIGNCZ(MXCENT),
     &             NCENTC(MXCENT),JSYMC(MXCENT),JCENTC(MXCENT),
     &             ICXVEC(MXCENT),ICYVEC(MXCENT),ICZVEC(MXCENT)
      COMMON/ADER/ADER0 (MXAQNS)
      COMMON/POINTER/IS0000, IS000X, IS000Y, IS000Z,
     &               IS00XX, IS00XY, IS00XZ, IS00YY,
     &               IS00YZ, IS00ZZ, IT0000, IT000X,
     &               IT000Y, IT000Z, IT00XX, IT00XY,
     &               IT00XZ, IT00YY, IT00YZ, IT00ZZ,
     &               ID0000, ID000X, ID000Y, ID000Z,
     &               ID00XX, ID00XY, ID00XZ, ID00YX,
     &               ID00YY, ID00YZ, ID00ZX, ID00ZY,
     &               ID00ZZ,
     &               IA0000, IA0X00, IA0Y00, IA0Z00,
     &               IAXX00, IAXY00, IAXZ00, IAYY00,
     &               IAYZ00, IAZZ00, IA000X, IA000Y,
     &               IA000Z, IA00XX, IA00XY, IA00XZ,
     &               IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     &               IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     &               IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON/DENPOIN/KDEN,KFOC
      COMMON/DSCFPOI/KDENSCF
      COMMON/GENCON/NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
      COMMON/METHOD/IUHF,SCF,NONHF
C
      DATA AZERO,HALF,TWO /0.0D0,0.5D0,2.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTXOR(I,J) = XOR(I,J)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
      IX=IPTAX(1)
      IY=IPTAX(2)
      IZ=IPTAX(3)
C
C GET DENSITY MATRIX ELEMENTS
C
      IOFF=0
      ISKIP=NRCA*NRCB
      IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
      DO 200 IRCA=1,NRCA
       MAXB=NRCB
       IF(LDIAG) MAXB=IRCA
       DO 200 IRCB=1,MAXB
        IOFF=IOFF+1
        ICOMP=IOFF
        DO 150 IORBA=IDENA+IRCA,IDENA+NRCA*KHKTA,NRCA
        DO 150 IORBB=IDENB+IRCB,IDENB+NRCB*KHKTB,NRCB
         DSHELL(ICOMP)=WORK1(KDEN-1+ITRI(IORBA,IORBB))
        ICOMP=ICOMP+ISKIP
150    CONTINUE
       IF(LDIAG.AND.(IRCA.EQ.IRCB)) THEN
        CALL SSCAL(KHKTA*KHKTB,HALF,DSHELL(IOFF),ISKIP)
        CALL SSCAL(KHKTA,TWO,DSHELL(IOFF),ISKIP*(KHKTA+1))
       ENDIF
200   CONTINUE
C
      MAXCMP = ISKIP*KHKTA*KHKTB
C
      DO 250 IATOM=1,NATOMC
C
       XAVRX=SDOT(MAXCMP,DSHELL,1,WORK1(IA0X0X+IATOM),NATOMC)
       XAVRY=SDOT(MAXCMP,DSHELL,1,WORK1(IA0X0Y+IATOM),NATOMC)
       XAVRZ=SDOT(MAXCMP,DSHELL,1,WORK1(IA0X0Z+IATOM),NATOMC)
       YAVRX=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Y0X+IATOM),NATOMC)
       YAVRY=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Y0Y+IATOM),NATOMC)
       YAVRZ=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Y0Z+IATOM),NATOMC)
       ZAVRX=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Z0X+IATOM),NATOMC)
       ZAVRY=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Z0Y+IATOM),NATOMC)
       ZAVRZ=SDOT(MAXCMP,DSHELL,1,WORK1(IA0Z0Z+IATOM),NATOMC)
C 
       ICENTC=NCENTC(IATOM)
       KCENTC=JCENTC(IATOM)
       NCX=3*KCENTC-2
       NCY=3*KCENTC-1
       NCZ=3*KCENTC
       SCX=SIGNCX(IATOM)
       SCY=SIGNCY(IATOM)
       SCZ=SIGNCZ(IATOM)
C
       ISYMPC=JSYMC(IATOM)
C
       DO 500 IREP = 0, MAXLOP
C
        CHI=PT(IBTAND(ISYMPC,IREP)) 
C
        CSCX = CHI*SCZ*SCY
        CSCY = CHI*SCX*SCZ
        CSCZ = CHI*SCX*SCY
        ICX = IPTCNT(NCX,IREP)
        ICY = IPTCNT(NCY,IREP)
        ICZ = IPTCNT(NCZ,IREP)
        IF (IBTXOR(ISYTYP(3),ISYTYP(2)) .EQ. IREP) THEN
         IF (ICX.GT.0)
     &   CSH(IX,ICX) = CSH(IX,ICX) + XAVRX
         IF (ICY.GT.0)
     &   CSH(IX,ICY) = CSH(IX,ICY) + XAVRY*CSCY
         IF (ICZ.GT.0)
     &   CSH(IX,ICZ) = CSH(IX,ICZ) + XAVRZ*CSCZ
       END IF
       IF (IBTXOR(ISYTYP(1),ISYTYP(3)) .EQ. IREP) THEN
        IF (ICX.GT.0)
     &  CSH(IY,ICX) = CSH(IY,ICX) + YAVRX*CSCX
        IF (ICY.GT.0)
     &  CSH(IY,ICY) = CSH(IY,ICY) + YAVRY
        IF (ICZ.GT.0)
     &  CSH(IY,ICZ) = CSH(IY,ICZ) + YAVRZ*CSCZ
       END IF
       IF (IBTXOR(ISYTYP(1),ISYTYP(2)) .EQ. IREP) THEN
        IF (ICX.GT.0)
     &  CSH(IZ,ICX) = CSH(IZ,ICX) + ZAVRX*CSCX
        IF (ICY.GT.0)
     &  CSH(IZ,ICY) = CSH(IZ,ICY) + ZAVRY*CSCY
        IF (ICZ.GT.0)
     &  CSH(IZ,ICZ) = CSH(IZ,ICZ) + ZAVRZ
       END IF
C
  500 CONTINUE
C 
  250 CONTINUE
C
      IF(.NOT.SCF) THEN
       IOFF=0
       ISKIP=NRCA*NRCB
       IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
       DO 210 IRCA=1,NRCA
        MAXB=NRCB
        IF(LDIAG) MAXB=IRCA
        DO 210 IRCB=1,MAXB
         IOFF=IOFF+1
         ICOMP=IOFF
         DO 160 IORBA=IDENA+IRCA,IDENA+NRCA*KHKTA,NRCA
         DO 160 IORBB=IDENB+IRCB,IDENB+NRCB*KHKTB,NRCB
          DSHELL2(ICOMP)=WORK1(KDENSCF-1+ITRI(IORBA,IORBB))
         ICOMP=ICOMP+ISKIP
160     CONTINUE
        IF(LDIAG.AND.(IRCA.EQ.IRCB)) THEN
         CALL SSCAL(KHKTA*KHKTB,HALF,DSHELL2(IOFF),ISKIP)
         CALL SSCAL(KHKTA,TWO,DSHELL2(IOFF),ISKIP*(KHKTA+1))
        ENDIF
210    CONTINUE
C
      MAXCMP = ISKIP*KHKTA*KHKTB
C
      DO 1250 IATOM=1,NATOMC
C
       XAVRX=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0X0X+IATOM),NATOMC)
       XAVRY=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0X0Y+IATOM),NATOMC)
       XAVRZ=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0X0Z+IATOM),NATOMC)
       YAVRX=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0Y0X+IATOM),NATOMC)
       YAVRY=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0Y0Y+IATOM),NATOMC)
       YAVRZ=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0Y0Z+IATOM),NATOMC)
       ZAVRX=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0Z0X+IATOM),NATOMC)
       ZAVRY=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0Z0Y+IATOM),NATOMC)
       ZAVRZ=SDOT(MAXCMP,DSHELL2,1,WORK1(IA0Z0Z+IATOM),NATOMC)
C 
       ICENTC=NCENTC(IATOM)
       KCENTC=JCENTC(IATOM)
       NCX=3*KCENTC-2
       NCY=3*KCENTC-1
       NCZ=3*KCENTC
       SCX=SIGNCX(IATOM)
       SCY=SIGNCY(IATOM)
       SCZ=SIGNCZ(IATOM)
C
       ISYMPC=JSYMC(IATOM)
C
       DO 1500 IREP = 0, MAXLOP
C
        CHI=PT(IBTAND(ISYMPC,IREP)) 
C
        CSCX = CHI*SCY*SCZ
        CSCY = CHI*SCX*SCZ
        CSCZ = CHI*SCX*SCY
        ICX = IPTCNT(NCX,IREP)
        ICY = IPTCNT(NCY,IREP)
        ICZ = IPTCNT(NCZ,IREP)
        IF (IBTXOR(ISYTYP(3),ISYTYP(2)) .EQ. IREP) THEN
         IF (ICX.GT.0)
     &   CSHSCF(IX,ICX) = CSHSCF(IX,ICX) + XAVRX
         IF (ICY.GT.0)
     &   CSHSCF(IX,ICY) = CSHSCF(IX,ICY) + XAVRY*CSCY
         IF (ICZ.GT.0)
     &   CSHSCF(IX,ICZ) = CSHSCF(IX,ICZ) + XAVRZ*CSCZ
       END IF
       IF (IBTXOR(ISYTYP(1),ISYTYP(3)) .EQ. IREP) THEN
        IF (ICX.GT.0)
     &  CSHSCF(IY,ICX) = CSHSCF(IY,ICX) + YAVRX*CSCX
        IF (ICY.GT.0)
     &  CSHSCF(IY,ICY) = CSHSCF(IY,ICY) + YAVRY
        IF (ICZ.GT.0)
     &  CSHSCF(IY,ICZ) = CSHSCF(IY,ICZ) + YAVRZ*CSCZ
       END IF
       IF (IBTXOR(ISYTYP(1),ISYTYP(2)) .EQ. IREP) THEN
        IF (ICX.GT.0)
     &  CSHSCF(IZ,ICX) = CSHSCF(IZ,ICX) + ZAVRX*CSCX
        IF (ICY.GT.0)
     &  CSHSCF(IZ,ICY) = CSHSCF(IZ,ICY) + ZAVRY*CSCY
        IF (ICZ.GT.0)
     &  CSHSCF(IZ,ICZ) = CSHSCF(IZ,ICZ) + ZAVRZ
       END IF
C
 1500 CONTINUE
C 
 1250 CONTINUE
C
      ENDIF
      RETURN
      END
