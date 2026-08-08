











      SUBROUTINE POST_MODH4(ICORE,MAXCOR)
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


      DATA ONE /1.0/

      TOTSIZ=ISYMSZ(ISYTYP(1,23),ISYTYP(1,23))
      I000=1
      I010=I000+TOTSIZ*IINTFP
      I020=I010+TOTSIZ*IINTFP 
      I030=I020+TOTSIZ*IINTFP

      IF(I020.GT.MAXCOR)CALL INSMEM('QUIKAA2',I020,MAXCOR)

      CALL GETALL(ICORE(I000),TOTSIZ,1,258)
      CALL GETALL(ICORE(I010),TOTSIZ,1,256)
      CALL GETALL(ICORE(I020),TOTSIZ,1,254)

C      CALL SAXPY(TOTSIZ,ONE,ICORE(I000),1,ICORE(I010),1)
C      CALL PUTALL(ICORE(I010),TOTSIZ,1,256)
C      CALL SAXPY(TOTSIZ,ONE,ICORE(I000),1,ICORE(I020),1)
C      CALL PUTALL(ICORE(I020),TOTSIZ,1,254)
C      CALL DZERO(ICORE(I000),TOTSIZ)
C      CALL PUTALL(ICORE(I000),TOTSIZ,1,258)

      call checksum("List-256:",ICORE(I010),TOTSIZ)
      call checksum("List-258:",ICORE(I000),TOTSIZ)
      call checksum("List-254:",ICORE(I020),TOTSIZ)
      RETURN
      END
