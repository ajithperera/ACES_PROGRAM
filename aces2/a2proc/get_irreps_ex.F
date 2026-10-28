C
C DRIVER FOR DETERMINATION OF ORBITAL SYMMETRIES. From ACES II. 
C 
      SUBROUTINE GET_IRREPS_EX(SCFVEC,SCFEVAL,ICORE,MAXCOR,NBAS,NBASX,
     &                         ISPIN, NOCC, IUHF, IFLAG, O_MOS, 
     &                         T_MOS, N_o_mos, N_t_mos, STATE)
C
      IMPLICIT INTEGER (A-Z)

      DIMENSION ICORE(MAXCOR), IREPS(9), NBFIRR(8), NOCC(8,2)
      INTEGER O_MOS(NBASX), T_MOS(NBASX)
      DOuble Precision Scfvec(Nbasx*Nbas), Scfeval(Nbas)
      CHARACTER*1 SP(2)
      CHARACTER*7 STATE
C
#include "machsp.com"
C
      DATA IONE /1/
      DATA SP /'A','B'/

C Write a Table of Eigenvalues and their symmetry. 

      IREPS(1) = 1
      DO IRREP = 1, NIRREP
         IREPS(IRREP+1) = IREPS(IRREP) + NBFIRR(IRREP)
      ENDDO 

      I020 = 1
      I030 = I020 + NBASX
      I040 = I030 + NBASX

      CALL GETREC(-1,'JOBARC','ANMOMBF0',NBASX,ICORE(I030))

      CALL EVLOUT_EX(SCFEVAL, ICORE(I020), ICORE(I030), IREPS, 
     &               NBAS, NBASX, NIRREP, NOCC, ISPIN, IFLAG, 
     &               O_MOS, T_MOS, N_o_mos, N_t_mos, STATE)

      RETURN
      END

