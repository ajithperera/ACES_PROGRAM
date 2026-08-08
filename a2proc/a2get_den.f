
















      SUBROUTINE A2GET_DEN(WORK, MAXCOR, IREORDER, DENSITY_TYPE, 
     &                     SPIN_D, NMBR_OF_PERTS, NBFNS, NAOBFNS, 
     &                     ISCF_TDEN, ICOR_TDEN, ISCF_DDEN, 
     &                     ICOR_DDEN, IBEGIN_P_DENS, IUHF) 
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPHERICAL, SPIN_D
C
c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      CHARACTER*7 DENSITY_TYPE
      CHARACTER*8 LABELPSCF, LABELDSCF, LABELDENS, LABELSDEN,
     &            RECNAMEA, RECNAMEB
      DIMENSION WORK(MAXCOR), IREORDER(NAOBFNS), ISCR1(MAXBASFN)
C
      DATA IONE /1/
C


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end








































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C
      SPHERICAL = (IFLAGS(62) .EQ. 1)
      LABELPSCF ='HFDENSTY'
      LABELDSCF ='HDDENSTY'
      LABELDENS ='TDENSITY'
      LABELSDEN ='DDENSITY'
      NAOBFNS2  = NAOBFNS*NAOBFNS
      NBFNS2    = NBFNS*NBFNS
C
      IBEGIN_P_DENS = IONE
      ISCF_TDEN     = IONE
      ISCF_DDEN     = IONE
      ICOR_TDEN     = IONE
      ICOR_DDEN     = IONE
C
      IF (DENSITY_TYPE .EQ. "0-ORDER") THEN
      
C
C Get the SCF and correlated (relaxed) density matirices in AO basis. For 
C open shell calculations, we also read the density diferences. Note that
C the "relaxed" density read here only contain the correlated contributions.
C
      ISCF_TDEN = IONE
      ICOR_TDEN = ISCF_TDEN + NAOBFNS2
      ISP2CART  = ICOR_TDEN + NAOBFNS2
      ITMP1     = ISP2CART  + NAOBFNS*NBFNS
      ISCR2     = ITMP1     + NAOBFNS*NAOBFNS
      INEXT     = ISCR2     + NAOBFNS
C
      IF (INEXT .GT. MAXCOR) CALL INSMEM("A2GET_DEN", INEXT, MAXCOR)

C The SCF density here does not have the trace property right 
C (multiplication from the S is needed), but this is minor issue 
C in the current context (the correlated density has the correct trace).

      CALL Getrec(-20, 'JOBARC', LABELPSCF, NBFNS2*IINTFP,
     &             WORK(ISCF_TDEN))
      CALL Getrec(-20, 'JOBARC', LABELDENS, NBFNS2*IINTFP, 
     &             WORK(ICOR_TDEN))
C  
      IF (SPHERICAL) THEN
         CALL Getrec(20, 'JOBARC','CMP2CART', NAOBFNS*NBFNS*IINTFP,
     &               WORK(ISP2CART))
      Write(6,*) "The spherical to Cartesian transformation matrix"
      CALL OUTPUT(WORK(ISP2CART), 1, NAOBFNS, 1, NBFNS, NAOBFNS, 
     &            NBFNS, 1)
      Write(6,*) "The SCF/Cor. tot. density matrices in SP AO basis"
      CALL OUTPUT(WORK(ISCF_TDEN), 1, NBFNS, 1, NBFNS, NBFNS,
     &            NBFNS, 1)
      CALL OUTPUT(WORK(ICOR_TDEN), 1, NBFNS, 1, NBFNS, NBFNS,
     &            NBFNS, 1)

C
        CALL TRANS_SPH2CART(WORK(ISCF_TDEN), WORK(ISP2CART),  
     &                      WORK(ITMP1), NAOBFNS, NBFNS)
C
        CALL TRANS_SPH2CART(WORK(ICOR_TDEN), WORK(ISP2CART),  
     &                      WORK(ITMP1), NAOBFNS, NBFNS)
      ENDIF
C
      Write(6,*) "The SCF/Cor. tot. density matrices in Car AO basis"
      CALL OUTPUT(WORK(ISCF_TDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      CALL OUTPUT(WORK(ICOR_TDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      
      IF (SPIN_D .AND. (IUHF .GT. 0)) THEN
         
          ISCF_DDEN = INEXT
          ICOR_DDEN = ISCF_DDEN + NBFNS2
          INEXT     = ICOR_DDEN + NBFNS2
          IF (INEXT .GT. MAXCOR) CALL INSMEM("A2GET_DEN", INEXT, 
     &                                        MAXCOR)
 
          CALL Getrec(-20, 'JOBARC', LABELDSCF, NBFNS2*IINTFP, 
     &                 WORK(ISCF_DDEN))
          CALL Getrec(-20, 'JOBARC', LABELSDEN, NBFNS2*IINTFP, 
     &                 WORK(ICOR_DDEN))
      ENDIF
C      
      Write(6,*) "The SCF and Cor. tot. density matrices in AO basis"
      CALL OUTPUT(WORK(ISCF_TDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
      CALL OUTPUT(WORK(ICOR_TDEN), 1, NAOBFNS, 1, NNAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
C
      IF (SPIN_D .AND. IUHF .GT. 0) THEN
      Write(6,*) "The SCF and Cor. spin. density matrices in AO basis"
      CALL OUTPUT(WORK(ISCF_DDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
      CALL OUTPUT(WORK(ICOR_DDEN), 1, NAOBFNS, 1, NAOBFNS, NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
      ENDIF
C
      ELSE IF (DENSITY_TYPE .EQ. "1-ORDER") THEN
C
           MAX_MEM_NEED = 2*NBFNS2*NMBR_OF_PERTS + 2*NBFNS2
           IF (MAXCOR .LE. MAX_MEM_NEED) CALL INSMEM("A2GET_DEN", 
     &                                                MAX_MEM_NEED, 
     &                                                MAXCOR)
C
           IKEEP_ONE_DENA = IONE 
           IKEEP_ONE_DENB = IKEEP_ONE_DENA + NBFNS2
           IKEEP_P_TDEN_N = IKEEP_ONE_DENB + NBFNS2
           IBEGIN_P_DENS  = IKEEP_P_TDEN_N  

           DO IPERT = 1, NMBR_OF_PERTS
C
              IKEEP_P_DDEN_N =  IKEEP_P_TDEN_N + NBFNS2
              INEXT          =  IKEEP_P_DDEN_N + NBFNS2
C
              WRITE(RECNAMEA, "(A,I2)") "PTDENA", IPERT
              CALL Getrec(-20, "JOBARC", RECNAMEA, NBFNS2*IINTFP, 
     &                    WORK(IKEEP_ONE_DENA))
C
              WRITE(RECNAMEB, "(A,I2)") "PTDENB", IPERT
              CALL Getrec(-20, "JOBARC", RECNAMEB, NBFNS2*IINTFP,
     &                     WORK(IKEEP_ONE_DENB))
C
              CALL ZERO(WORK(IKEEP_P_TDEN_N), NBFNS2)
              CALL ZERO(WORK(IKEEP_P_DDEN_N), NBFNS2)
C
C Build perturbed total density for n perturbation 
C
              CALL SAXPY(NAOBFNS2, 1.0D0, WORK(IKEEP_ONE_DENA), 1, 
     &                   WORK(IKEEP_P_TDEN_N), 1) 
              CALL SAXPY(NAOBFNS2, 1.0D0, WORK(IKEEP_ONE_DENB), 1, 
     &                   WORK(IKEEP_P_TDEN_N), 1) 
C
C Build perturbed spin density for n perturbation 
C
              CALL SAXPY(NAOBFNS2, 1.0D0, WORK(IKEEP_ONE_DENA), 1, 
     &                   WORK(IKEEP_P_DDEN_N), 1) 
              CALL SAXPY(NAOBFNS2, -1.0D0, WORK(IKEEP_ONE_DENB), 1, 
     &                   WORK(IKEEP_P_DDEN_N), 1) 
C
              IKEEP_P_TDEN_N =  INEXT
C
           ENDDO
C
C endif for density type
C
      ENDIF
C
      RETURN
      END


