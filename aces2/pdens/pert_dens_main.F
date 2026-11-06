      SUBROUTINE PERT_DENS_MAIN(ICORE, MAXCOR, IUHF)
C
#include "flags.h"
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      INTEGER POP, VRT, DIRPRD
      LOGICAL SPN_SPN, POLAR, JFC, JPSO, JSD, LINTDMP
C
#include "machsp.com"
#include "mxatms.par"
#include "jodaflags.com"
C
      DIMENSION ICORE(MAXCOR), LENGTH(8), LENGTH_SD(8, 2, 3), 
     &          POL_TENSOR(3, 3), SPNC_TENSOR(6*MXATMS*6*MXATMS),
     &          NPERT(8), NPRTRB(4)
C
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SPN_SPN_TYPE/ JFC, JPSO, JSD
      COMMON /PDDMP/ IUDMP, LINTDMP
C
      DATA IONE /1/      
C
#ifdef _MKLST
      CALL MAKE_LISTS4_PERTRB_DENS(ICORE, MAXCOR, IUHF)
#endif 
      DO IRREP = 1, NIRREP    
C
        LENGTH(IRREP) = 0
C
        DO ISPIN = 1, 1 + IUHF
                        LENGTH(IRREP) = LENGTH(IRREP) + 
     &                                  IRPDPD(IRREP, 8 + ISPIN)
           LENGTH_SD(IRREP, 1, ISPIN) = IRPDPD(IRREP, 8 + ISPIN)
        END DO
C
        DO ISPIN = 3, 3-2*IUHF, -1
                               LSTTYP = 43 + ISPIN
                        LENGTH(IRREP) = LENGTH(IRREP) + 
     &                                  IDSYMSZ(IRREP,
     &                                  ISYTYP(1, LSTTYP), 
     &                                  ISYTYP(2, LSTTYP))
           LENGTH_SD(IRREP, 2, ISPIN) = IDSYMSZ(IRREP,
     &                                  ISYTYP(1, LSTTYP), 
     &                                  ISYTYP(2, LSTTYP))
        END DO    
      END DO
C 
#ifdef _DEBUG_LVLM2
      print*, "Length of the perturbed T2 vector"  
      print*, (length(i), i=1, nirrep)
      print*, "Length of the perturbed singles"
      print*, ((length_sd(irrep, 1, i), irrep =1, nirrep), 
     &          i=1, iuhf+1)
      print*, "Length of perturrbed doubles"
      print*, ((length_sd(irrep, 2, i), irrep =1, nirrep), 
     &          i=3, 3-2*IUHF, -1 )
#endif
#ifdef _TDENS
      CALL GETREC(0, 'JOBARC', 'NTOTPERT', IRECLEN, NTPERT)  
      IF (IRECLEN .GT. 0) THEN
          CALL GETREC(20, 'JOBARC', 'NTOTPERT', 1, NTPERT)
          LISTPT2 = 373 
          CALL ZERO(SPNC_TENSOR, 9*MXATMS*MXATMS)
      ELSE
          NTPERT  = 3
          LISTPT2 = 373
          CALL ZERO(POL_TENSOR, 9)
      ENDIF 
      CALL GETREC(-1, 'JOBARC', 'NPERTIRR', 8, NPERT)
          LISTR1 = 490
       LISTR2_AA = 444 
       LISTR2_BB = 445
       LISTR2_AB = 446
#endif
#ifdef _DEBUG_LVLM
       print*, "The total number of perturbations and symmetries", 
     &          NTPERT
       Print*,  (NPERT(I), I=1, NIRREP)
       Print*, "The perturbed amplitudes are read from list:", 
     &          LISTPT2
#endif
#ifdef _TEST1
         LISTPT2 = 374
          LISTR1 = 90      
       LISTR2_AA = 44
       LISTR2_BB = 45
       LISTR2_AB = 46
#endif
C
C
C JPERT RUNS OVER ALL PERTURBATIONS, IRREP BY IRREP (VCCEH ORDER);
C IPERT IS ONLY THE INDEX WITHIN IRREPX.
C
      POLAR   = .FALSE.
      SPN_SPN = .FALSE.
      JPERT   = 0
C
C ANALYSIS FILE (SEE PDDUMP): RESPONSE DENSITIES, PROPERTY INTEGRALS
C AND BASIS INFORMATION.
C
      IUDMP   = 77
      LINTDMP = .TRUE.
      OPEN(UNIT=IUDMP, FILE='pdens_dump.bin', FORM='UNFORMATTED',
     &     ACCESS='STREAM', STATUS='REPLACE')
      NPRTRB(1) = IFLAGS(18)
      NPRTRB(2) = NTPERT
      NPRTRB(3) = NIRREP
      NPRTRB(4) = IUHF
      CALL PDDUMPI('HEADER  ', 0, 4, NPRTRB)
      DO IRREPX = 1, NIRREP
C
         IF (NPERT(IRREPX) .GT. 0) CALL NEWLST(IRREPX, ICORE, MAXCOR,
     &                                         IUHF)
C
#ifdef _DEBUG_LVL0
         IF (NPERT(IRREPX) .GT. 0 .AND. IUHF .NE. 0)
     &      CALL VCCEH_XPROP(ICORE, MAXCOR/IINTFP, LENGTH(IRREPX),
     &                       NPERT(IRREPX), IRREPX)
#endif
         DO IPERT = 1,  NPERT(IRREPX)
            JPERT = JPERT + 1

            IOFF_T1_A  = IONE
            IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX, 1, 1))*IINTFP
            IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX, 1, 1) +
     &                           LENGTH_SD(IRREPX, 1, 2))*IINTFP

C
C THE VCCEH VECTOR ON LIST 373 IS T1(A), T1(B), T2(AB), T2(BB),
C T2(AA) (SEE LOADVEC1/LOADVEC2).
C
            IF (IUHF .NE. 0) THEN
              IOFF_T1_B  = IONE + (LENGTH_SD(IRREPX,1,1))*IINTFP
              IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2))*IINTFP
              IOFF_T2_BB = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2) +
     &                             LENGTH_SD(IRREPX,2,3))*IINTFP
              IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2) +
     &                             LENGTH_SD(IRREPX,2,3) +
     &                             LENGTH_SD(IRREPX,2,2))*IINTFP
            ENDIF
C          
            IOFFSET = IPERT
C
            ILOC_PERT = IONE
            IEND_PERT = ILOC_PERT + LENGTH(IRREPX)*IINTFP
              MEMLEFT = (MAXCOR - IEND_PERT)/IINTFP 
C
            IF (IEND_PERT .GE. MAXCOR) CALL INSMEM("@-WAVEFN_ANALYSE",
     &                                             ILOC_PERT, MAXCOR)
#ifdef _GETPTS
#ifdef _DEBUG_LVLM0
        print*, "Reading the perturbed amplitudes"
#endif
            cALL GETLST(ICORE(ILOC_PERT), IOFFSET, 1, 1, IRREPX,
     &                  LISTPT2)

            CALL GET_RNORM(IRREPX, LENGTH(IRREPX), ICORE(ILOC_PERT),
     &                     ICORE(IEND_PERT), MEMLEFT, IUHF, RNORM)
#endif
C
#ifdef _DEBUG_LVLM2
            nsizec = length(irrepx)
            print*, "The length of the vector: ", nsizec, ILOC_PERT,
     &      ioffset
            call checkSum("Comp_pert_dens", icore(iloc_pert), nsizec)
            print*, "The pert. T amps"
            call output(icore(iloc_pert),1,nsizec, 1, 1,nsizec,1,1)
#endif
C 
            ILEN_T1_A = LENGTH_SD(IRREPX, 1, 1)
            ILEN_T2_AA = LENGTH_SD(IRREPX, 2, 1)
            ILEN_T2_AB = LENGTH_SD(IRREPX, 2, 3)

            IF (IUHF .NE. 0) THEN
               ILEN_T1_B = LENGTH_SD(IRREPX, 1, 2)
              ILEN_T2_BB = LENGTH_SD(IRREPX, 2, 2)
             ENDIF
C
#ifdef _DEBUG_LVL0
            print*, "The singles length AI and ai: ",
     &      ilen_T1_A, ilen_T1_B
            write(6,*)
            print*, "The offsets for T 1s: ", Ioff_t1_a,
     &      ioff_t1_b
            print*, "The perturbed T1s AI"
            call output(icore(ioff_t1_a), 1,ilen_t1_a,1,1,ilen_t1_a,1,1)
            call checksumP("T1AA", icore(ioff_t1_a),ilen_t1_a)
            print*, "The perturbed T1s ai"
            call output(icore(ioff_t1_b), 1,ilen_t1_b,1,1,ilen_t1_b,1,1)
            call checksumP("T1BB", icore(ioff_t1_b),ilen_t1_b)
#endif
C
C PUT T(1) WHERE PTDENS READS THE RIGHT VECTOR: R1 ON 490, R2 ON
C 444-446 (RETYPED FOR IRREPX BY NEWLST ABOVE).
C
            CALL PUTLST(ICORE(IOFF_T1_A), 1, 1, 1, 1, LISTR1)
            IF (IUHF .NE. 0) CALL PUTLST(ICORE(IOFF_T1_B), 1, 1, 1,
     &                                   2, LISTR1) 
            CALL PUTALL(ICORE(IOFF_T2_AB), LENGTH_SD(IRREPX,2,3),
     &                  IRREPX, LISTR2_AB)
            IF (IUHF .NE. 0) THEN
               CALL PUTALL(ICORE(IOFF_T2_AA), LENGTH_SD(IRREPX,2,1),
     &                     IRREPX, LISTR2_AA)
               CALL PUTALL(ICORE(IOFF_T2_BB), LENGTH_SD(IRREPX,2,2),
     &                     IRREPX, LISTR2_BB)
            ENDIF
C
#ifdef _TEST
            CALL GENRSPNSDEN(ICORE, MAXCOR, IUHF) 
#endif
C
#ifdef _TDENS
#ifdef _DEBUG_LVLM2
            print*, "entering  GEN_PERTRB_RSPNS_DENS"
#endif
            CALL GEN_PERTRB_RSPNS_DENS(ICORE, MAXCOR, IRREPX, IUHF)
#endif
            CALL PRP_DEN4PROPS(ICORE, MAXCOR, NAO, NATOMS, NCENTER, 
     &                         JPERT, IRREPX, IUHF)
C
            IF (JPERT .EQ. 1) THEN
               CALL GETREC(20, 'JOBARC', 'AOOVRLAP', NAO*NAO*IINTFP,
     &                     ICORE(1+NAO*NAO*IINTFP))
               CALL PDDUMP('OVERLAP ', 0, NAO*NAO,
     &                     ICORE(1+NAO*NAO*IINTFP))
               CALL GETREC(20, 'JOBARC', 'CENTERBF', NAO,
     &                     ICORE(1+NAO*NAO*IINTFP))
               CALL PDDUMPI('CENTERBF', 0, NAO,
     &                      ICORE(1+NAO*NAO*IINTFP))
               CALL GETREC(20, 'JOBARC', 'ANGMOMBF', NAO,
     &                     ICORE(1+NAO*NAO*IINTFP))
               CALL PDDUMPI('ANGMOMBF', 0, NAO,
     &                      ICORE(1+NAO*NAO*IINTFP))
C SCF MO COEFFICIENTS AND ENERGIES (ALPHA) FOR MO-PAIR ANALYSES; THE
C AOS ARE THOSE OF THE DENSITIES, AS IN MO2AO3.
               CALL GETREC(20, 'JOBARC', 'SCFEVECA', NAO*NAO*IINTFP,
     &                     ICORE(1+NAO*NAO*IINTFP))
               CALL PDDUMP('SCFEVECA', 0, NAO*NAO,
     &                     ICORE(1+NAO*NAO*IINTFP))
               CALL GETREC(20, 'JOBARC', 'SCFEVALA', NAO*IINTFP,
     &                     ICORE(1+NAO*NAO*IINTFP))
               CALL PDDUMP('SCFEVALA', 0, NAO,
     &                     ICORE(1+NAO*NAO*IINTFP))
            ENDIF
            CALL PDDUMP('DENSITY ', JPERT, NAO*NAO, ICORE(1))
C
               I000 = IONE
              ILEFT = I000 + NAO*NAO*IINTFP
            MEMLEFT = MAXCOR - ILEFT
                                     
            CALL COMP_2NDORD_PROPS(ICORE(I000), ICORE(ILEFT), NTPERT,
     &                            MEMLEFT, NAO, NATOMS, JPERT, IUHF,
     &                            SPN_SPN, POLAR, SPNC_TENSOR, 
     &                            POL_TENSOR)
            LINTDMP = .FALSE.
C
C RESPONSE NATURAL ORBITAL DECOMPOSITION OF THE FC RESPONSE (NO
C SYMMETRY: THE PERTURBATIONS ARE THEN THE INDIVIDUAL NUCLEI).
C
            IF (NIRREP .EQ. 1 .AND. JFC)
     &         CALL RSPNS_NO(ICORE(I000), ICORE(ILEFT), MEMLEFT/IINTFP,
     &                       NAO, NATOMS, JPERT, IUHF)
#ifdef _BLKCHK
C BLOCK-BY-BLOCK CHECK AGAINST VCCEH (REBUILDS THE DENSITY SIX TIMES
C PER PERTURBATION; RUN WITH EOMPROP=LINEAR SO LIST 374 IS M(L)).
            IF (NIRREP .EQ. 1 .AND. (SPN_SPN .OR. POLAR))
     &         CALL BLKCHK(ICORE, MAXCOR, IUHF, IRREPX, IPERT,
     &                     NPERT(IRREPX), NTPERT)
#endif
         ENDDO
      ENDDO
C
      CLOSE(IUDMP)
      IUDMP = 0
C
      IF (SPN_SPN) THEN
         IF (NIRREP .NE. 1) THEN
            I000 = IONE
            I010 = I000 + NTPERT*NTPERT*IINTFP
            I020 = I010 + NTPERT*NTPERT*IINTFP
            IF (I020 .GT. MAXCOR) CALL INSMEM("@-PERT_DENS_MAIN",
     &                                        I020, MAXCOR)
#ifdef _DEBUG_LVL0
            Write(6,"(a)") " Raw atom x SALC response matrix"
            CALL OUTPUT(SPNC_TENSOR, 1, NTPERT, 1, NTPERT,
     &                  NTPERT, NTPERT, 1)
#endif
            CALL SALC2ATM(SPNC_TENSOR, ICORE(I000), ICORE(I010),
     &                    NTPERT)
         ENDIF
#ifdef _DEBUG_LVL0
         Write(6,"(a)") " Atom x atom response matrix before symmet."
         CALL OUTPUT(SPNC_TENSOR, 1, NTPERT, 1, NTPERT,
     &               NTPERT, NTPERT, 1)
#endif
#ifdef _DEBUG_LVLM2
         CALL OUTPUT(SPNC_TENSOR, 1, NTPERT, 1, NTPERT, 
     &               NTPERT, NTPERT, 1)
#endif
         CALL SYMMET2(SPNC_TENSOR, NTPERT)
         CALL DSCAL(NTPERT*NTPERT, 2.0D0, SPNC_TENSOR, 1)
#ifdef _DEBUG_LVL0
         CALL OUTPUT(SPNC_TENSOR, 1, NTPERT, 1, NTPERT, 
     &               NTPERT, NTPERT, 1)
#endif
         CALL PRNT_SPNSPN_TNSR(SPNC_TENSOR, NTPERT, JFC, JSD, JPSO,
     &                         .FALSE., 1)
      ELSE IF (POLAR) THEN
CSSS         CALL SYMMET2(POL_TENSOR, 3) 
         CALL DSCAL(9, 2.0D0, POL_TENSOR, 1)
         CALL NEATPRINT(6, POL_TENSOR, 3, 3, 3, 3)
      ENDIF

      RETURN 
      END
















