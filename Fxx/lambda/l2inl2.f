










      SUBROUTINE L2INL2(ICORE,MAXCOR,IUHF,IDOPPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
C
C THIS ROUTINE CALCULATES THE L2 CONTRIBUTION TO THE L2 INCREMENT
C
      DIMENSION ICORE(MAXCOR)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ iFlags2(500)
      COMMON /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &             CC2
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
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
      IF(IUHF.EQ.1)IBOT=1
      IF(IUHF.EQ.0)IBOT=3

      IF (.NOT. CC2) THEN

      ITOP=6
      IF(IDOPPL.NE.0)ITOP=1
       DO 10 ITYPE=1,ITOP,5
        IF(ITYPE.EQ.6.AND.IFLAGS(93).EQ.2)THEN
         CALL DRAOLAD(ICORE,MAXCOR,IUHF,.TRUE.,1,0,143,60,
     &                243,260)
        ELSE

CSSS#ifdef _DCC_FLAG
        If (Ispar) Then
           IF (ITYPE .EQ. 1) then
               LISTIN = 250
               CALL PDCC_L2LAD(ICORE,MAXCOR,IUHF,ITYPE,LISTIN)
           Else
               CALL L2LAD(ICORE,MAXCOR,IUHF,ITYPE)
           Endif 
        Else
CSSS#else
          CALL L2LAD(ICORE,MAXCOR,IUHF,ITYPE)
        Endif 
CSSS#endif 
        ENDIF
10     CONTINUE

      ENDIF

      Write(6,*) "The lambda residuals after l2lads"
      call check_leom(Icore,Maxcor,Iuhf)
     
         DO 15 ISPIN=IBOT,3
          CALL L2RNG(ICORE,MAXCOR,ISPIN,ITYPE,IUHF)
15       CONTINUE

      RETURN
      END
