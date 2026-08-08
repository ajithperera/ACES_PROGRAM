










      SUBROUTINE PERT_DENS_MAIN(ICORE, MAXCOR, IUHF)
C



































































































































































































C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      INTEGER POP, VRT, DIRPRD
      LOGICAL SPN_SPN, POLAR, JFC, JPSO, JSD
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



C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)







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
      DIMENSION ICORE(MAXCOR), LENGTH(8), LENGTH_SD(8, 2, 3), 
     &          POL_TENSOR(3, 3), SPNC_TENSOR(6*MXATMS*6*MXATMS),
     &          NPERT(8)
C
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SPN_SPN_TYPE/ JFC, JPSO, JSD
C
      DATA IONE /1/      
C
      CALL MAKE_LISTS4_PERTRB_DENS(ICORE, MAXCOR, IUHF)
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
C
      DO IRREPX = 1, NIRREP
C
         DO IPERT = 1,  NPERT(IRREPX)

            IOFF_T1_A  = IONE
            IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX, 1, 1))*IINTFP
            IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX, 1, 1) +
     &                           LENGTH_SD(IRREPX, 1, 2))*IINTFP

            IF (IUHF .NE. 0) THEN
              IOFF_T1_B  = IONE + (LENGTH_SD(IRREPX,1,1))*IINTFP
              IOFF_T2_AB = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2))*IINTFP
              IOFF_T2_AA = IONE + (LENGTH_SD(IRREPX,1,1) +
     &                             LENGTH_SD(IRREPX,1,2) +
     &                             LENGTH_SD(IRREPX,2,3))*IINTFP
              IOFF_T2_BB = IONE + (LENGTH_SD(IRREPX,1,1) +
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
            cALL GETLST(ICORE(ILOC_PERT), IOFFSET, 1, 1, IRREPX,
     &                  LISTPT2)

            CALL GET_RNORM(IRREPX, LENGTH(IRREPX), ICORE(ILOC_PERT),
     &                     ICORE(IEND_PERT), MEMLEFT, IUHF, RNORM)
C
C 
            ILEN_T1_A = LENGTH_SD(IRREPX, 1, 1)
            ILEN_T2_AA = LENGTH_SD(IRREPX, 2, 1)
            ILEN_T2_AB = LENGTH_SD(IRREPX, 2, 3)

            IF (IUHF .NE. 0) THEN
               ILEN_T1_B = LENGTH_SD(IRREPX, 1, 2)
              ILEN_T2_BB = LENGTH_SD(IRREPX, 2, 2)
             ENDIF
C
            print*, "The singles length AI and ai: ",
     &      ilen_T1_A, ilen_T1_B
            write(6,*)
            print*, "The offsets for T 1s: ", Ioff_t1_a,
     &      ioff_t1_b
            print*, "The perturbed T1s AI"
            call output(icore(ioff_t1_a), 1,nt(1), 1,1,nt(1),1,1)
            call checksumP("T1AA", icore(ioff_t1_a),nt(1))
            print*, "The perturbed T1s ai"
            call output(icore(ioff_t1_b), 1,nt(2), 1,1,nt(2),1,1)
            call checksumP("T1BB", icore(ioff_t1_b),nt(2))
     
            IOFF_T1_A = IOFF_T1_A + ILEN_T1_A*IINTFP    
            IF (IUHF .NE . 0) IOFF_T1_B = IOFF_T1_B + ILEN_T1_B*
     &                                    IINTFP
            DO IRREP = 1, NIRREP
               IRREPR = IRREP
               IRREPL = DIRPRD(IRREPR, IRREPX) 

               LENGTH_OOAA = IRPDPD(IRREPR, 3)
               LENGTH_VVAA = IRPDPD(IRREPL, 1)
               LENGTH_AAAA = LENGTH_OOAA*LENGTH_VVAA
               IF (IUHF .NE . 0) THEN
                   LENGTH_OOBB = IRPDPD(IRREPR, 4)
                   LENGTH_VVBB = IRPDPD(IRREPL, 2)
                   LENGTH_BBBB = LENGTH_OOBB*LENGTH_VVBB
               ENDIF
               LENGTH_OOAB = IRPDPD(IRREPR, 14)
               LENGTH_VVAB = IRPDPD(IRREPL, 13)
               LENGTH_ABAB = LENGTH_OOAB*LENGTH_VVAB
C
C
              IOFF_T2_AB = IOFF_T2_AB +  LENGTH_ABAB*IINTFP
              IOFF_T2_AA = IOFF_T2_AA +  LENGTH_AAAA*IINTFP
              IF (IUHF .NE. 0) IOFF_T2_BB = IOFF_T2_BB + 
     &                                      LENGTH_BBBB*IINTFP

            END DO
C
C
            CALL GEN_PERTRB_RSPNS_DENS(ICORE, MAXCOR, IRREPX, IUHF)
            CALL PRP_DEN4PROPS(ICORE, MAXCOR, NAO, NATOMS, NCENTER, 
     &                         IPERT, IRREPX, IUHF)
C
               I000 = IONE
              ILEFT = I000 + NAO*NAO*IINTFP
            MEMLEFT = MAXCOR - ILEFT
                                     
            CALL COMP_2NDORD_PROPS(ICORE(I000), ICORE(ILEFT), NTPERT,
     &                            MEMLEFT, NAO, NATOMS, IPERT, IUHF,
     &                            SPN_SPN, POLAR, SPNC_TENSOR, 
     &                            POL_TENSOR)
         ENDDO
      ENDDO
C
      IF (SPN_SPN) THEN
         CALL SYMMET2(SPNC_TENSOR, NTPERT)
         CALL DSCAL(NTPERT*NTPERT, 2.0D0, SPNC_TENSOR, 1)
         CALL OUTPUT(SPNC_TENSOR, 1, NTPERT, 1, NTPERT, 
     &               NTPERT, NTPERT, 1)
         CALL PRNT_SPNSPN_TNSR(SPNC_TENSOR, NTPERT, JFC, JSD, JPSO,
     &                         .FALSE., 1)
      ELSE IF (POLAR) THEN
CSSS         CALL SYMMET2(POL_TENSOR, 3) 
         CALL DSCAL(9, 2.0D0, POL_TENSOR, 1)
         CALL NEATPRINT(6, POL_TENSOR, 3, 3, 3, 3)
      ENDIF

      RETURN 
      END
















