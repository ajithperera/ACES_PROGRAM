











C     This program is supposed to compute fifth/sixth-order noniterative
C     triple excitation contributions beyond CCSD(T). It assumes a file
C     of second-order T3 amplitudes is available, as well, of course, as
C     T1 and T2 amplitudes.

      PROGRAM NONIT3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP1,POP2,VRT1,VRT2,DIRPRD,DISTSZ
      LOGICAL PRESNT
      LOGICAL CCSDT4,CCSDT
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
c trps.par : begin

c NUMT3 := the maximum number of triples amplitudes allowed for
c          incore processing
      INTEGER NUMT3
      PARAMETER (NUMT3=5000000)

c trps.par : end
C
      DIMENSION E45T(20)
C
      COMMON / / ICORE(1)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /CONTROL/ IPRNT,IXXX,IXXX2
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /CORENG/ ELAST
      EQUIVALENCE (IFLAGS( 2),METHOD)
      EQUIVALENCE (IFLAGS(38),INONHF)
      EQUIVALENCE (IFLAGS(39),IORBTP)
C      
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
 
c YAU - not used in this file
c#include "icdacc.com"

      COMMON/SYM/POP1(8),POP2(8),VRT1(8),VRT2(8),
     &           NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3ENRG/ E4TAAA,E4TBBB,E4TAAB,E4TBBA,
     1                E5STA, E5STB,
     1                E4STA ,E4STB,
     1                E5TD1 ,E5TD4 ,E5TD2 ,E5TD3,
     1                E5TT1 ,E5TT4 ,E5TT2 ,E5TT3,
     1                E5QT1 ,E5QT4 ,E5QT2 ,E5QT3,
     1                E6TT1 ,E6TT4 ,E6TT2 ,E6TT3
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /T3CORE/ INCORE
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /T3METH/ CCSDT4,CCSDT
      COMMON /T3FILE/ ALLT3(NUMT3)
      COMMON /AUXIO/ DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1               NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
C
      CALL CRAPSI(ICORE,IUHF,0)
      MAXCOR=ICRSIZ
C
      WRITE(LUOUT,1000)
 1000 FORMAT(/,72('-'),/,T15,
     1       ' Noniterative 5th Order Triples Calculation ',/,
     1         72('-'),/)
c1000 FORMAT(' @NONIT3-I, Welcome to NONIT3. ')
C
C     Determine what kind of calculation we are supposed to be doing :
C
C     1. Some kind of noniterative triple calculation based on CCSD,
C        QCISD, or CCSDT T1 and
C        T2 amplitudes with HF or semicanonical orbitals. This will
C        also involve quadruples, which are evaluated with another
C        module (xvcc5q).
C
C     2. An "iterative-noniterative triple" calculation, for example,
C        ROHF-CCSD(T) with standard ROHF orbitals (not currently implem-
C        ented).
C
      IF(INONHF.GT.0.AND.IORBTP.EQ.0)THEN
      WRITE(LUOUT,1010)
 1010 FORMAT(' @NONIT-I, Performing an iterative-noniterative T3',
     1       ' calculation. ')
      ENDIF
C
Calso allow for Brueckner soon !
      IF(INONHF.EQ.0)THEN
      WRITE(LUOUT,1020)
 1020 FORMAT(' @NONIT-I, Computing noniterative T3 contributions. ')
      ENDIF
C
      IF(METHOD.EQ.22) WRITE(LUOUT,1030)
      IF(METHOD.EQ.12) WRITE(LUOUT,1040)
      IF(METHOD.EQ.26) WRITE(LUOUT,1050)
      IF(METHOD.EQ.27) WRITE(LUOUT,1060)
      IF(METHOD.EQ.28) WRITE(LUOUT,1070)
      IF(METHOD.EQ.29) WRITE(LUOUT,1080)
      IF(METHOD.EQ.30) WRITE(LUOUT,1090)
      IF(METHOD.EQ.31) WRITE(LUOUT,1095)
 1030 FORMAT(' @NONIT-I, Method is CCSD(T). ')
 1040 FORMAT(' @NONIT-I, Method is CCSD+TQ*(CCSD). ')
 1050 FORMAT(' @NONIT-I, Method is QCISD(TQ). ')
 1060 FORMAT(' @NONIT-I, Method is CCSD(TQ). ')
 1070 FORMAT(' @NONIT-I, Method is CCSD+TQ(CCSD). ')
 1080 FORMAT(' @NONIT-I, Method is CCSDT+Q*(CCSDT). ')
 1090 FORMAT(' @NONIT-I, Method is CCSDT+Q(CCSDT). ')
 1095 FORMAT(' @NONIT-I, Method is CC5SD[T].')
C
      IF(METHOD.NE.22.AND.METHOD.NE.12.AND.METHOD.LT.26.AND.
     1                                     METHOD.GT.31)THEN
C
C     At the moment we handle either CCSD(T) iterative-noniterative
C     calculations or various fifth-order extensions of CCSD, QCISD,
C     and CCSDT.
C
C     Something has gone terribly wrong !
C
      WRITE(LUOUT,1200) METHOD
 1200 FORMAT(' @NONIT-F, This module does not yet support ICLLVL ',I4)
      CALL ERREX
      ENDIF
C
C     Open T3 file (but do not initialize it) and initialize
C     the various i/o arrays. Later we need logic in AUXIOI and other
C     places for handling only one set of lists on disk in simple
C     noniterative calculations (me thinks we have already done this).
C
      PRESNT = .TRUE.
      CALL AUXIOI(PRESNT,IUHF)
C
C     Create ABCI and IJKA intermediates on 107-110, 127-130.
C     For the purpose of TRPINT, we pretend that this is a CCSDT-2
C     calculation.
C
      METACT = METHOD
      METHOD = 15
      CALL TRPINT_VCC5T(ICORE(I0),MAXCOR,IUHF,METACT)
C
      METHOD = METACT
C
      CCSDT   = .FALSE.
      CCSDT4  = .FALSE.
      TRIPNI  = .FALSE.
      TRIPNI1 = .FALSE.
      TRIPIT  = .FALSE.
      T3STOR  = .TRUE.
      IF(METHOD.EQ.22)THEN
       TRIPIT  = .TRUE.
      ENDIF
      IF(METHOD.EQ.12.OR.METHOD.EQ.26.OR.METHOD.EQ.27.OR.METHOD.EQ.28
     &               .OR.METHOD.EQ.31)THEN
       TRIPNI = .TRUE.
       CCSDT4 = .TRUE.
      ENDIF
      IF(METHOD.EQ.29.OR.METHOD.EQ.30)THEN
       TRIPNI = .TRUE.
      ENDIF
      CALL SETLWI(METHOD)
C
C     Compute S1 quantities.
C
      I000 = I0
      I010 = I000 + NTAA * IINTFP
      I020 = I010 + NTBB * IINTFP
      MAXCOR= MAXCOR - IINTFP * (NTAA + NTBB)
      CALL S1WT2(ICORE(I000),ICORE(I010),ICORE(I020),MAXCOR,IUHF)
C
      CALL TRPS_VCC5T(ICORE(I020),MAXCOR,IUHF,E4T,ICORE(I000),
     &                ICORE(I010))
C
      IF(METHOD.EQ.12.OR.METHOD.EQ.26.OR.METHOD.EQ.27.OR.METHOD.EQ.28
     &               .OR.METHOD.EQ.31)THEN
C
C     The RHF branch below folds in the all-alpha (AAA) spin case as
C     E5xx1 -- but TRPS14, the block that computes those, is skipped for
C     RHF at trps.F "IF(IUHF.EQ.0.OR.ISKP14) GOTO 1000", so E5TD1/E5TT1/
C     E5QT1/E6TT1 arrive here as exact zeros and the 2*(X1+X2) below
C     silently drops half of each fifth-order term. It is not a small
C     effect: for the H2O/DZP CCSD(TQ) reference case (tests/fast 020)
C     the missing 2*E5TD1 etc. add up to 2.24e-04 hartree, and the run
C     reports -76.2737454652 where the archive (and a UHF run of the same
C     input, and a CFOUR CCSDTQ/CCSDT(Q) ladder in the same basis) all say
C     -76.2739691209. Refuse rather than write that number to JOBARC.
C
C     What a real fix needs: TRPS14 reads the AAAA OOOV intermediate
C     (list 107 = LWIC11) and the AAAA T2 amplitudes (list 44 = 43+ISPIN),
C     and an RHF run has neither -- T2-AAAA is allocated but left at zero
C     (see the @CHECKSUM T2-AAAA line in any RHF output), and several AAAA
C     lists are not even dimensioned, since trpint.F only INIPCKs them
C     under IF(IUHF.NE.0). So RHF support means initialising those lists
C     and filling them by antisymmetrising their ABAB counterparts
C     (libr/assym.f is the existing helper) before letting the block run;
C     simply deleting the guard in trps.F fails immediately with
C     "@GETTRN: Assertion failed, MOIO rows = 0". Verification targets are
C     known exactly, from a UHF run of case 020: E5TD1 = -0.000057010064,
C     E5TT1 = 0.000000815536, E5QT1 = 0.000002037487,
C     E6TT1 = -0.000000661806, total -76.273969124744.
C
C     Until then the workaround is REF=UHF, which is correct for a closed
C     shell too and reproduces the reference to 3.8e-09.
C
       IF(IUHF.EQ.0)THEN
       WRITE(6,9110)
 9110  FORMAT(/,' @NONIT3-F, Fifth-order triples/quadruples',
     &        ' corrections are not available',/,
     &        '            with an RHF reference: the all-alpha spin',
     &        ' case is never computed',/,
     &        '            (TRPS14 is skipped for IUHF=0), which would',
     &        ' silently drop about',/,
     &        '            half of each fifth-order term. Rerun with',
     &        ' REF=UHF.',/)
       CALL ERREX
       ENDIF
       IF(IUHF.EQ.0)THEN
       E5TD = 2.0D+00 * (E5TD1 + E5TD2)
       E5TT = 2.0D+00 * (E5TT1 + E5TT2)
       E5QTD= 2.0D+00 * (E5QT1 + E5QT2)
       E6TT = 2.0D+00 * (E6TT1 + E6TT2)
       ELSE
       E5TD =            E5TD1 + E5TD4 + E5TD2 + E5TD3
       E5TT =            E5TT1 + E5TT4 + E5TT2 + E5TT3
       E5QTD=            E5QT1 + E5QT4 + E5QT2 + E5QT3
       E6TT =            E6TT1 + E6TT4 + E6TT2 + E6TT3
       ENDIF
c      write(6,*) ' nonit-i, e5td, e5tt ',E5TD,E5TT
c      write(6,*) e5td1,e5tt1
c      write(6,*) e5td4,e5tt4
c      write(6,*) e5td2,e5tt2
c      write(6,*) e5td3,e5tt3
C
      CALL GETREC(20,'JOBARC','NONITT3',20*IINTFP,E45T)
      E4T     = E45T(1)
      E5ST    = E45T(2)
      E4DT    = E45T(3)
      E5TD    = E5TD - E4T
ctest E5TT    = E5TT - E5TD - E4T
      E45T(4) = E5TD
      E45T(5) = E5TT
      E5QTC   = E5TD
      E5QT    = E5QTC + E5QTD
      E45T(6) = E5QT
      E45T(7) = E6TT
      CALL PUTREC(20,'JOBARC','NONITT3',20*IINTFP,E45T)
      WRITE(LUOUT,2000) E4T,E5ST,E4DT,E5TD,E5TT,E5QTC,E5QTD,E5QT,E6TT
 2000 FORMAT(' @NONIT-I, Noniterative triple excitation energies :',/,
     1       '           E4T                         ',F20.12,/,
     1       '           E5ST                        ',F20.12,/,
     1       '           E4DT                        ',F20.12,/,
     1       '           E5TD                        ',F20.12,/,
     1       '           E5TT                        ',F20.12,/,
     1       '           E5QT(c)                     ',F20.12,/,
     1       '           E5QT(d)                     ',F20.12,/,
     1       '           E5QT                        ',F20.12,/,
     1       '           E6TT                        ',F20.12)
C
      IF(METHOD.EQ.12.OR.METHOD.EQ.27.OR.METHOD.EQ.28.OR.
     &   METHOD.EQ.31)THEN
      CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
      WRITE(LUOUT,2010) ETOT
 2010 FORMAT(' @NONIT-I, CCSD(T)  energy ',F20.12,' a.u. ')
      ETOT = ETOT + E5TD + E5TT
      WRITE(LUOUT,2020) ETOT
 2020 FORMAT(' @NONIT-I, CC5SD(T) energy ',F20.12,' a.u. ')
      CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
      ENDIF
C      
      ENDIF
C
      IF(METHOD.EQ.29.OR.METHOD.EQ.30)THEN

       IF(IUHF.EQ.0)THEN
       E5TD = 2.0D+00 * (E5TD1 + E5TD2)
c       E5TT = 2.0D+00 * (E5TT1 + E5TT2)
       E5QTD= 2.0D+00 * (E5QT1 + E5QT2)
c       E6TT = 2.0D+00 * (E6TT1 + E6TT2)
       ELSE
       E5TD =            E5TD1 + E5TD4 + E5TD2 + E5TD3
c       E5TT =            E5TT1 + E5TT4 + E5TT2 + E5TT3
       E5QTD=            E5QT1 + E5QT4 + E5QT2 + E5QT3
c       E6TT =            E6TT1 + E6TT4 + E6TT2 + E6TT3
       ENDIF
C
      CALL GETREC(20,'JOBARC','NONITT3',20*IINTFP,E45T)
c      E4T     = E45T(1)
c      E5ST    = E45T(2)
c      E4DT    = E45T(3)
c      E5TD    = E5TD - E4T
ctest E5TT    = E5TT - E5TD - E4T
      E45T(4) = E5TD
c      E45T(5) = E5TT
      E5QTC   = E5TD
      E5QT    = E5QTC + E5QTD
      E45T(6) = E5QT
c      E45T(7) = E6TT
      CALL PUTREC(20,'JOBARC','NONITT3',20*IINTFP,E45T)
      WRITE(LUOUT,2000) E4T,E5ST,E4DT,E5TD,E5TT,E5QTC,E5QTD,E5QT,E6TT
C
      CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
      WRITE(LUOUT,2030) ETOT
 2030 FORMAT(' @NONIT-I, CCSDT energy ',F20.12,' a.u. ')
      ETOT = ETOT + E5QT
      WRITE(LUOUT,2040) ETOT
 2040 FORMAT(' @NONIT-I, CCSDT + E5QT energy ',F20.12,' a.u. ')
      CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
      ENDIF
C      
      CALL AUXIOO(.TRUE.)
      call aces_fin
      STOP
      END
