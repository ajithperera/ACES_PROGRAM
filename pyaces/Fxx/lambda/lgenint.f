












      SUBROUTINE LGENINT(ICORE,MAXCOR,IUHF)
C
C DRIVER FOR THE FORMATION OF THE V AND G INTERMEDIATES IN SOLVING THE
C LAMBDA-CC EQUATIONS
C
CEND
C
C  CODED AUGUST/90  JG
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TCPU,TSYS,ONE
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      DOUBLE PRECISION  TIMEIN, TIMENOW, TIMETOT, TIMENEW
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

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
      DATA ONE /1.0D0/
C
      IF(IFLAGS(1).GE.10)THEN
       WRITE(LUOUT,100)
100    FORMAT(T3,'@LGENINT-I, Calculating intermediates.')
       CALL TIMER(1)
      ENDIF
C
C HOLE-HOLE LADDER INTERMEDIATE : V(MN,IJ)
C
CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
         CALL PDCC_FORMV1(ICORE,MAXCOR,IUHF)
CSSS#else
      ELse
         CALL FORMV1(ICORE,MAXCOR,IUHF) 
      Endif 
CSSS#endif 

      IF(IFLAGS(1).GE.10)CALL TIMER(1)
C
C G-INTERMEDIATES
C
CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
         IRREPX      = 1
         IHHA_LENGTH = IRPDPD(IRREPX,21)
         IHHB_LENGTH = IRPDPD(IRREPX,22)
         IPPA_LENGTH = IRPDPD(IRREPX,19)
         IPPB_LENGTH = IRPDPD(IRREPX,20)

         Call zerlst(ICORE,IHHA_LENGTH,1,1,1,191)
         Call zerlst(ICORE,IPPA_LENGTH,1,1,1,192)
         If (iuhf .ne. 0) then
            call zerlst(ICORE,IHHB_LENGTH,1,1,2,191)
            call zerlst(ICORE,IPPB_LENGTH,1,1,2,192)
         Endif

         CALL FORMG1(ICORE,MAXCOR,IUHF,Gae_scale)
         CALL FORMG2(ICORE,MAXCOR,IUHF,Gmi_scale)
       Else
CSSS#else
         CALL FORMG1(ICORE,MAXCOR,IUHF,ONE)
         CALL FORMG2(ICORE,MAXCOR,IUHF,ONE)
       Endif 
CSSS#endif 
      IF(IFLAGS(1).GE.10)THEN
       CALL TIMER(1)
       WRITE(LUOUT,102)TIMENEW
102    FORMAT(T3,'@GENINT-I, G intermediates required ',F9.3,
     &           ' seconds.')
      ENDIF
      RETURN
      END
