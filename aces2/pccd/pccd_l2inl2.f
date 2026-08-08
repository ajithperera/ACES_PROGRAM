










      SUBROUTINE PCCD_L2INL2(ICORE,MAXCOR,IUHF,IDOPPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION ICORE(MAXCOR)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ iFlags2(500)
C

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
      IBOT=3
      ITOP=6
      IF(IDOPPL.NE.0) ITOP=1

      DO 10 ITYPE=1,ITOP,5
         IF(ITYPE.EQ.6.AND.IFLAGS(93).EQ.2)THEN
            CALL DRAOLAD(ICORE,MAXCOR,IUHF,.TRUE.,1,0,143,60,
     &                   243,260)
         ELSE
            CALL PCCD_L2LAD(ICORE,MAXCOR,IUHF,ITYPE)
        ENDIF

10    CONTINUE

      Write(6,*) "The lambda residuals after l2lads"
      call check_leom(Icore,Maxcor,Iuhf)
      CALL PCCD_L2RNG(ICORE,MAXCOR,3,ITYPE,IUHF)

      Write(6,*) "The lambda residuals after l2rng"
      call check_leom(Icore,Maxcor,Iuhf)

      RETURN
      END
