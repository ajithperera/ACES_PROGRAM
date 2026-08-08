










C
C *******************************************************
C THE FOLLOWING PROCEDURES FORM THE CORE OF THE IP-EOM PROGRAM. 
C THEY EVALUATE THE PRODUCT OF THE IP-EOM MATRIX TIMES AN ARBITRARY 
C VECTOR RESIDING ON LS1IN AND LS2IN.
C *******************************************************
C
      SUBROUTINE HC_MULT(ICORE,MAXCOR,IUHF, ISIDE, KSPIN)
C
C THIS SUBROUTINE CALCULATES THE MATRIX VECTOR PRODUCT
C OF THE IP-EOM MATRIX AND A VECTOR
C
C INPUT : COMMON SLISTS
C LS1IN, 
C LS2IN : THESE LISTS CONTAIN
C         THE H AND 2HP COMPONENTS OF THE INPUT VECTOR
C         LS2IN IS LABELED BY (ISPIN,IMIXSPIN) -> S(JI,AP)
C         J AND A HAVE IMIXSPIN, I HAS ISPIN
C         ISIDE:   = 1:  IPMAT * C
C                  = 2:  C * IPMAT
C
C THE COMMONBLOCK SINFO CONTAINS FURTHER INPUT INFORMATION
C
C NS  : THE NUMBER OF INPUT VECTORS OF REPRESENTATION SIRREP
C LENS: THE LENGTH OF THE 2HP PART OF THE S-VECTOR
C
C OUTPUT: LISTS1EX, LISTS2EX: LISTS THAT CONTAIN THE RESPECTIVE 
C         COMPONENTS OF THE EA-EOM MATRIX TIMES THE INPUT VECTOR
CEND
      IMPLICIT INTEGER (A-Z)

      LOGICAL LEFTHAND, SINGONLY,DROPCORE

      DOUBLE PRECISION ONE, ONEM
      DOUBLE PRECISION TIN, TOUT,TIMDUM
      DOUBLE PRECISION TSSTGEN, TEADAVID, TEADIR
C
      DIMENSION ICORE(MAXCOR)
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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
      COMMON/IPCALC/LEFTHAND,SINGONLY, DROPCORE
      COMMON /TIMSUB/ TSSTGEN, TEADAVID, TEADIR
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW      
C
      DATA ONE, ONEM /1.0D0, -1.0D0/
C
      CALL TIMER(1)

      ISTART = KSPIN
      IEND   = KSPIN

      Write(6,"(a)") "     ----------------HC-mult----------------"
C
      DO 10 ISPIN = 1,1+IUHF
         CALL SETZERO_IP(IUHF,ISPIN,ICORE,MAXCOR,LS1OUT,
     $                   LS2OUT)
 10   CONTINUE

CSSS         Call Check_ipvecs(Iuhf,Kspin,Icore,Maxcor,Ls1in,Ls2in)

      IF ((STMODE .NE. 0) .AND. (ISIDE .EQ. 2)) THEN
         CALL SETZERO_IP(IUHF,3-KSPIN,ICORE,MAXCOR,LS1IN,
     $                   LS2IN)
      ENDIF

      DO ISPIN = ISTART, IEND 

         CALL S1HHS1(ICORE, MAXCOR, ISIDE, ISPIN)

         CALL S2HPHHS1(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

         CALL S2HPS1(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

         CALL S1PHHHS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

C Only for the left hand side (S1HPS2) 

         CALL S1HPS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

         IF (ISIDE .EQ. 1) THEN 
            CALL S2HHS2_1R(ICORE, MAXCOR, IUHF, ISPIN)
            CALL S2HHS2_2R(ICORE, MAXCOR, IUHF, ISPIN)
         ELSE
            CALL S2HHS2_1L(ICORE, MAXCOR, IUHF, ISPIN)
            CALL S2HHS2_2L(ICORE, MAXCOR, IUHF, ISPIN)
         ENDIF 

         CALL S2PPS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

         CALL S2HHHHS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

         CALL S2PHPHS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

         CALL S2HPPHS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)

         CALL S2HHPPT2S2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      ENDDO

      CALL TIMER(1)      
      TEADIR = TEADIR + TIMENEW
C
      RETURN
      END
