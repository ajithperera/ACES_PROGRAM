










      
      SUBROUTINE PDCC_FORMH4(ICORE,MAXCOR,IUHF)
C
C   THIS ROUTINE CALCULATE L(MN,EF) T(NJ,FB) => (MB,EJ) for
C   Distinguish cluster approximation. 
C
C  CCD
C   
C   H(ME,JB) = - SUM N,F  L(MN,EF) T(NJ,FB)
C
C  QCISD
C
C   H(ME,JB) = - SUM N,F L(MN,EF) T(NJ,FB)
C
C  CCSD
C
C   G(IA,JB) = - SUM N,F L(MN,EF) (T(NJ,FB) + ...)
C
C  THIS TERM IS VERY SIMILAR TO THE T1(IJ,AB) CONTRIBUTION TO
C  THE W-RING INTERMEDIATE
C
C  THE SPIN CASES ARE
C
C    AAAA : =  - SUM M,E T1(IM,BE) T1(MJ,EA) - SUM m,e T1(Im,Be) T1(Jm,Ae)
C
C    ABAB : =  - SUM m,E T1(Im,Eb) T1(Jm,Ea)
C
C    ABBA : =  - SUM M,E T1(IM,BE) T1(Mj,Ea) - SUM m,e T1(Im,Be) T1(mj,eb)
C
C FOR UHF IN ADDITION THE BBBB BABA AND BABA SPIN CASES HAS TO CALCULATES
C
C
C  AAAA : LIST 254  (CALCULATED FROM ABAB AND ABBA IN QUIKAA)
C  BBBB : LIST 255  (UHF ONLY)
C  ABBA : LIST 258
C  BAAB : LIST 259  (UHF ONLY)
C  ABAB : LIST 256
C  BABA : LIST 257  (UHF ONLY)
C
CEND
C 
C  CODED AUGUST/90  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DEBUG
      DIMENSION ICORE(MAXCOR)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

C
      write(6,"(a)") " ---Entered PDCC_FORMH4---"

      Imode  = 0
      Irrepx = 1

      Call Inipck(Irrepx,14,14,253,Imode,0,1)
      Call Inipck(Irrepx,9,9,254,Imode,0,1)
      Call Inipck(Irrepx,9,10,256,Imode,0,1)
      Call Inipck(Irrepx,11,11,258,Imode,0,1)

      If (Iuhf .ne. 0) Then
         Call Inipck(Irrepx,10,10,255,Imode,0,1)
         Call Inipck(Irrepx,10,9,257,Imode,0,1)
         Call Inipck(Irrepx,12,12,259,Imode,0,1)
      Endif

      Call Inipck(Irrepx,13,14,216,Imode,0,1)
      Call Inipck(Irrepx,13,14,217,Imode,0,1)
      If (Iuhf .ne. 0) Then
         Call Inipck(Irrepx,1,3,214,Imode,0,1)
         Call Inipck(Irrepx,2,4,215,Imode,0,1)
      Endif 
 

      CALL PDCC_H4ALL(ICORE,MAXCOR,'ABBA',IUHF)
      CALL PDCC_H4ALL(ICORE,MAXCOR,'ABAB',IUHF)
      CALL PDCC_H4ALL(ICORE,MAXCOR,'AAAA',IUHF)

CSSS      CALL PDCC_QUIKAA2(ICORE,MAXCOR,IUHF)

      IF (IUHF .NE. 0) THEN
         CALL PDCC_H4ALL(ICORE,MAXCOR,'BBBB',IUHF)
         CALL PDCC_H4ALL(ICORE,MAXCOR,'BABA',IUHF) 
      ENDIF 

      CALL POST_MODH4(ICORE,MAXCOR,IUHF)

      RETURN
      END
