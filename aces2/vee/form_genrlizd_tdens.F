      Subroutine Form_genrlizd_tdens(Scr,Maxcor,Iuhf)
C
C Generalized (state-to-state) EOM-CC transition densities
C
C     D(pq; k,l) = <L_k| {p+ q} |R_l>
C
C for every pair of converged EOM-EE states (R state l, L state k),
C including pairs from different irreps. The R and L vectors are
C saved by nextdav/edens_4gtdens on list 472 (slot 3): R of state n
C on record 17+2n, L on record 18+2n; RNORM(2n-1) is r0 of state n
C and LNORM(2n-1) the factor that biorthonormalizes L_n to R_n.
C
C Output: one JOBARC record per R state l, "GTDENSll", holding the
C spin-summed AO-basis transition densities for all L states k
C (NTOTAL blocks of NAO*NAO, in the same state order as the loops
C below: irreps, then roots within an irrep). Spin summed means
C 2*D(alpha) for RHF and D(alpha)+D(beta) for UHF, so that
C <L_k|O|R_l> = Sum_pq O(pq) D(pq; k,l) for any one-electron
C spin-free operator O in the AO basis.
C
C 11/2026 rewrite (Ajith Perera): per-pair list shapes (the resorted
C L/R lists and the G intermediates were left in the shape of the
C last Davidson irrep), densities sized and expanded with the pair
C symmetry IRREPY x IRREPX, overlap term ZNORM = delta(k,l), MO -> AO
C transformation and spin sum added, and one record per R state
C (the old one-record-per-pair-and-spin layout overran the 32
C labels once there were more than 5 RHF / 4 UHF states).
C
      Implicit Double Precision(A-H,O-Z)
C
      DIMENSION SCR(MAXCOR)
      DOUBLE PRECISION LNORM(200),RNORM(200)
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      CHARACTER*2 SUFFIX(32)
      CHARACTER*8 LABEL
C
#include "flags.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"
#include "machsp.com"
C
      COMMON /CALCINFO/ NROOT(8)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
C
      DATA SUFFIX /"01","02","03","04","05","06","07","08","09","10",
     +             "11","12","13","14","15","16","17","18","19","20",
     +             "21","22","23","24","25","26","27","28","29","30",
     +             "31","32"/
      DATA ONE,ZILCH /1.0D0,0.0D0/
C
      I00=100
      CALL GETREC(20,"JOBARC","RNORM   ",I00*IINTFP,RNORM)
      CALL GETREC(20,"JOBARC","LNORM   ",I00*IINTFP,LNORM)
C
      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      NMO =NOCCO(1)+NVRTO(1)
      NAO2=NAO*NAO
C
      NTOTAL=0
      DO IRREP=1,NIRREP
         NTOTAL=NTOTAL+NROOT(IRREP)
      ENDDO
      IF (NTOTAL .EQ. 0) RETURN
C
C Largest packed density blocks and EOM vector over all irreps
C
      MXOO =0
      MXVV =0
      MXVO =0
      MAXLEN=0
      DO IRREP=1,NIRREP
         MXOO=MAX(MXOO,IRPDPD(IRREP,21)+IUHF*IRPDPD(IRREP,22))
         MXVV=MAX(MXVV,IRPDPD(IRREP,19)+IUHF*IRPDPD(IRREP,20))
         MXVO=MAX(MXVO,IRPDPD(IRREP, 9)+IUHF*IRPDPD(IRREP,10))
         MAXLEN=MAX(MAXLEN,LENVEC(IRREP,IUHF))
      ENDDO
C
      IDOO =1
      IDVV =IDOO +MXOO
      IDVO =IDVV +MXVV
      IDOV =IDVO +MXVO
      IDMO =IDOV +MXVO
      IGT  =IDMO +NMO*NMO
      ITMP1=IGT  +NTOTAL*NAO2
      ITMP2=ITMP1+NAO2
      ITMP3=ITMP2+NAO*MAX(NAO,NMO)
      ITOP =ITMP3+NAO2
      IF (ITOP+MAXLEN .GT. MAXCOR) CALL
     +   INSMEM("@-form_genrlizd_tdens",ITOP+MAXLEN,MAXCOR)
      MXCOR=MAXCOR-ITOP+1
C
C RHF: the spin-adapted density is the alpha block; total = 2*alpha
C
      FACTSP=DFLOAT(2-IUHF)
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The R and L norms "
      Write(6,"(5(1x,F15.10))") (RNORM(I),I=1,NTOTAL*2)
      Write(6,*)
      Write(6,"(5(1x,F15.10))") (LNORM(I),I=1,NTOTAL*2)
#endif
C
      IDONE=0
      DO IRREPX=1,NIRREP
      DO IROOT =1,NROOT(IRREPX)
         IDONE=IDONE+1
C
C Right vector of state IDONE: singles on 3,490 - 4,490, doubles on
C 461-463, resorted doubles on 454-459, all with symmetry IRREPX.
C
         IGET=17+2*IDONE
         R0  =RNORM(2*IDONE-1)
         IF (IRREPX .NE. 1) R0=ZILCH
         IF (MBPT2 .OR. LCCD .OR. LCCSD .OR. RCCD .OR. DRCCD
     &             .OR. CC2) R0=ZILCH
         CALL GETLST(SCR(ITOP),IGET,1,1,3,472)
#ifdef _DEBUG_LVL0
         Write(6,*)
         write(6,"(2a,(1x,i2))") " Reading the right vector",
     +                           " of state: ",iget
         Call checksum("list-472      :",scr(itop),
     +                 LENVEC(IRREPX,IUHF),s)
#endif
         CALL SETVEC(IRREPX,2,461,IUHF)
         CALL UPDATES(IRREPX,SCR(ITOP),461,2,490,IUHF)
         CALL SETRES(IRREPX,454,IUHF)
         CALL RESORT(SCR(ITOP),IINTFP*MXCOR,IUHF,IRREPX,461,454)
C
         CALL ZERO(SCR(IGT),NTOTAL*NAO2)
C
         JDONE=0
         DO IRREPY=1,NIRREP
         DO JROOT =1,NROOT(IRREPY)
            JDONE=JDONE+1
C
C Left vector of state JDONE, biorthonormalized to its own R: singles
C on 1,490 - 2,490, doubles on 444-446, resorted on 434-439 (IRREPY).
C
            JGET=18+2*JDONE
            Z   =LNORM(2*JDONE-1)
            CALL GETLST(SCR(ITOP),JGET,1,1,3,472)
            CALL DSCAL(LENVEC(IRREPY,IUHF),Z,SCR(ITOP),1)
#ifdef _DEBUG_LVL0
            write(6,"(2a,(1x,i2))") " Reading the left vector",
     +                              " of state : ", Jget
            Call checksum("list-472      :",scr(itop),
     +                    LENVEC(IRREPY,IUHF),s)
#endif
            CALL SETVEC(IRREPY,0,444,IUHF)
            CALL UPDATES(IRREPY,SCR(ITOP),444,0,490,IUHF)
            CALL SETRES(IRREPY,434,IUHF)
            CALL RESORT(SCR(ITOP),IINTFP*MXCOR,IUHF,IRREPY,444,434)
C
C G intermediates: L2*T2 (symmetry IRREPY) on 1,491/492 - 2,491/492,
C L2*R2 (symmetry IRREPXY) on 3,491/492 - 4,491/492, and the
C DOV - R0*L1 scratch (IRREPXY) on 1,493 - 2,493. All of these are
C allocated in maklst with the largest size over irreps.
C
            IRREPXY=DIRPRD(IRREPY,IRREPX)
            DO ISPIN=1,1+IUHF
               CALL ACES_LIST_RESIZE(ISPIN,  491,
     &                               IRPDPD(IRREPY, 20+ISPIN))
               CALL ACES_LIST_RESIZE(ISPIN,  492,
     &                               IRPDPD(IRREPY, 18+ISPIN))
               CALL ACES_LIST_RESIZE(ISPIN+2,491,
     &                               IRPDPD(IRREPXY,20+ISPIN))
               CALL ACES_LIST_RESIZE(ISPIN+2,492,
     &                               IRPDPD(IRREPXY,18+ISPIN))
               CALL ACES_LIST_RESIZE(ISPIN,  493,
     &                               IRPDPD(IRREPXY, 8+ISPIN))
            ENDDO
C
C The disconnected <L_k|R_l> * T1 term is delta(k,l).
C
            IF (JDONE .EQ. IDONE) THEN
               ZNORM=ONE
            ELSE
               ZNORM=ZILCH
            ENDIF
#ifdef _DEBUG_LVL0
            Write(6,"(a,2(1x,i2),3F12.8)") "k,l,Z,R0,ZNORM: ",
     +            JDONE,IDONE,Z,R0,ZNORM
#endif
            CALL GDENS(IRREPY,IRREPX,IRREPXY,SCR(IDOO),SCR(IDVV),
     &                 SCR(IDVO),SCR(IDOV),SCR(ITOP),MXCOR,IUHF,
     &                 ONE,R0,ZILCH,ZNORM,
     &                 400,400,2,0,493,0,
     &                 490,490,2,0,461,444,454,434,90,0)
C
            LOO=IRPDPD(IRREPXY,21)
            LVV=IRPDPD(IRREPXY,19)
            LVO=IRPDPD(IRREPXY, 9)
            IOFFGT=IGT+(JDONE-1)*NAO2
            DO ISPIN=1,1+IUHF
               IOFFOO=IDOO+(ISPIN-1)*LOO
               IOFFVV=IDVV+(ISPIN-1)*LVV
               IOFFVO=IDVO+(ISPIN-1)*LVO
               IOFFOV=IDOV+(ISPIN-1)*LVO
               CALL EXPDEN(SCR(IOFFOO),SCR(IOFFVV),SCR(IOFFVO),
     &                     SCR(IOFFOV),SCR(IDMO),NMO,IRREPXY,ISPIN,
     &                     .FALSE.)
               CALL MO2AO3(SCR(IDMO),SCR(ITMP1),SCR(ITMP2),
     &                     SCR(ITMP3),NAO,NMO,ISPIN)
               CALL SAXPY(NAO2,FACTSP,SCR(ITMP1),1,SCR(IOFFGT),1)
            ENDDO
#ifdef _DEBUG_LVL0
            Write(6,"(a,4(1x,i2))") "Irrepx,Iroot,Irrepy,Jroot: ",
     +                              Irrepx,Iroot,Irrepy,Jroot
            Call Checksum("GTDENS(AO)",SCR(IOFFGT),NAO2,s)
#endif
         ENDDO
         ENDDO
C
         LABEL="GTDENS"//SUFFIX(IDONE)
         CALL PUTREC(20,"JOBARC",LABEL,NTOTAL*NAO2*IINTFP,SCR(IGT))
      ENDDO
      ENDDO
C
      RETURN
      END
