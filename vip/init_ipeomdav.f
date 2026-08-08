










      SUBROUTINE INIT_IPEOMDAV(NIP, SIRREP, ISPIN, IUHF, SCR, MAXCOR)
C
C  PARAMETERS FOR DAVIDSON PROCEDURE ARE INITIALIZED. ULTIMATELY THIS
C  MUST BE REWRITTEN AND SPECIFIED BY INPUT
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR, EIGVAL, THRESH, R, P, ONE, ROOT
      LOGICAL EMINFOL,EVECFOL, LEFTHAND, EXCICORE
C
      PARAMETER (MAXORD=100)
      PARAMETER (MAXROOT=100)
C
      DIMENSION SCR(MAXCOR)
C
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTINF2/ROOT
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
C
      COMMON/RMAT/ R(10000), P(10000)
      COMMON/ROOTS/EIGVAL(100,8,3), OSCSTR(100,8,3)
      COMMON/CNVRGE/EMINFOL,EVECFOL
      COMMON/IPCALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
C
      DATA ONE /1.0D0/
C
      THRESH = -100.0
      STCALC = 2
      NSIZEC = NIP
      NDIMR = 1
      IOLDEST = 1
C
      IROOT = 0
C
      CALL ZERO(R, 10000)
      CALL ZERO(P, 10000)
C
C CALCULATE THE NEEDED EXCITATION PATTERNS
C
C COLUMN 5 OF H0  CONTAINS THE INTERESTING COMPONENTS
C COLUMN 6 OF H0  CONTAINS THE ACTUAL INCLUDED COMPONENTS
C
      CALL CALCEXCP(ISPIN, IUHF, SCR, MAXCOR, 5, LISTH0)
      CALL CALCEXCP(ISPIN, IUHF, SCR, MAXCOR, 6, LISTH0)
C
      CALL CALCH0(ISPIN, IUHF, SCR, MAXCOR)
C
C LOAD GUESS FOR FIRST VECTOR
C
C  THE INITIAL GUESS IS BASED ON THE DIAGONAL MATRIX H0
C
        I000 = 1
        I010 = I000 + NSIZEC
        I020 = I010 + NSIZEC
        CALL GETLST(SCR(I000), 4,1,1,1, LISTH0)

        Write(6,"(a)") " The digonals of Hbar_0 (eigenvalues)"
        Write(6,"(6(1x,F15.9))"),(scr(i000+j-1),j=1,nsizec)
        CALL FNDMINE(NSIZEC, SCR(I000), EIGVAL(1,SIRREP,ISPIN),
     $               SCR(I010),0, ROOT, ILOC, 1.D-4, THRESH)
        SCR(ILOC) = 1.0D30
        CALL PUTLST(SCR(I000), 4, 1,1,1, LISTH0)
C
        CALL ZERO(SCR,NSIZEC)
        SCR(ILOC) = ONE
        CALL NORMVEC(SCR(I000), NSIZEC, SCR(I010), MAXCOR-I010+1,
     $     IUHF, 3, SIRREP)
        WRITE(6,'(1X,A,E12.5,I6)') 'GUESS FOR FIRST EIGENVALUE: ',
     $     ROOT, ILOC

C
        CALL PUTS_HHP(SCR(I000),NSIZEC,ISPIN, IUHF, LS1IN, LS2IN)
C
C  THIS GUESS MAY ALSO BE USED FOR LEFTHAND, AND IS PUT ON 
C  COLUMN 3 OF LISTH0
C
        CALL PUTLST(SCR(I000),3,1,1,1,LISTH0)
C
        RETURN
        END
