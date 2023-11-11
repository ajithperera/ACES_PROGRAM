      SUBROUTINE PCCD_QUIKAA2(ICORE,MAXCOR)
C
C THIS SUBROUTINE FORMS THE ALL ALPHA H(MEJB) INTERMEDIATE
C  FROM THE ABAB AND ABBA PIECES FOR RHF REFERENCE FUNCTIONS.
C
CEND
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER TOTSIZ
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

      DATA ONEM /-1.0/
      IOFFLISTH=0
      IF(IFLAGS(3).EQ.2) IOFFLISTH=200
      TOTSIZ=ISYMSZ(ISYTYP(1,23),ISYTYP(1,23))
      I000=1
      I010=I000+TOTSIZ*IINTFP
      I020=I010+TOTSIZ*IINTFP 
      IF(I020.GT.MAXCOR)CALL INSMEM('QUIKAA2',I020,MAXCOR)
      CALL GETALL(ICORE(I000),TOTSIZ,1,158+IOFFLISTH)

      Write(6,"(a)") "@-QUIKAA2:"
      call checksum("158    :", ICORE(I000), TOTSIZ)
      CALL GETALL(ICORE(I010),TOTSIZ,1,156+IOFFLISTH)
      call checksum("156    :", ICORE(I010), TOTSIZ)
      CALL SAXPY(TOTSIZ,ONEM,ICORE(I010),1,ICORE(I000),1)
      CALL PUTALL(ICORE,TOTSIZ,1,154+IOFFLISTH)

      call checksum("154    :", ICORE, TOTSIZ)
      RETURN
      END
