










      Subroutine ZeroHbar(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     


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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      UHF    = .FALSE.
      UHF    = (IUhf .EQ. 1)
      IRREPX = 1

C Hbar(AI,BC)
      IF (UHF) Then
         AAAA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
         BBBB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
         ABAB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
         LENGTH = MAX(AAAA_LENGTH_AIBC,BBBB_LENGTH_AIBC,
     +                ABAB_LENGTH_AIBC,ABBA_LENGTH_AIBC)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_AIBC, IRREPX, 27)
         Call Putall(Work, BBBB_LENGTH_AIBC, IRREPX, 28)
         Call Putall(Work, ABAB_LENGTH_AIBC, IRREPX, 29)
         Call Putall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
      Else
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
         LENGTH = ABBA_LENGTH_AIBC
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
      Endif 

C Hbar(IJ,KA)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
         LENGTH = MAX(AAAA_LENGTH_IJKA,BBBB_LENGTH_IJKA,ABAB_LENGTH_IJKA,
     +                ABBA_LENGTH_IJKA)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 7)
         Call Putall(Work, BBBB_LENGTH_IJKA, IRREPX, 8)
         Call Putall(Work, ABAB_LENGTH_AIBC, IRREPX, 9)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
      Else
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
         LENGTH = MAX(AAAA_LENGTH_IJKA,ABBA_LENGTH_IJKA)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 7)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
      Endif

C Hbar(MB,EJ)
      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
	 ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         LENGTH = MAX(AAAA_LENGTH_MBEJ,BBBB_LENGTH_MBEJ,
     +                AABB_LENGTH_MBEJ,BBAA_LENGTH_MBEJ,
     +                ABAB_LENGTH_MBEJ,BABA_LENGTH_MBEJ)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         LENGTH = MAX(AAAA_LENGTH_MBEJ,AABB_LENGTH_MBEJ,
     +                ABAB_LENGTH_MBEJ)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
      Endif

C Hbar(AB,CI)
      IF (UHF) Then
         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
         BBBB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,128),ISYTYP(2,128))
         ABAB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,129),ISYTYP(2,129))
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))
         LENGTH = MAX(AAAA_LENGTH_ABCI,BBBB_LENGTH_ABCI,
     +                ABAB_LENGTH_ABCI,ABBA_LENGTH_ABCI)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_ABCI, IRREPX, 127)
         Call Putall(Work, BBBB_LENGTH_ABCI, IRREPX, 128)
         Call Putall(Work, ABAB_LENGTH_ABCI, IRREPX, 129)
         Call Putall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)
      Else
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))
         LENGTH = ABBA_LENGTH_ABCI
         CALL DZERO(WORK,LENGTH)
         Call Getall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)
      Endif

C Hbar(IA,JK)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,108),ISYTYP(2,108))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,109),ISYTYP(2,109))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))
         LENGTH = MAX(AAAA_LENGTH_IJKA,BBBB_LENGTH_IJKA,
     +                ABAB_LENGTH_IJKA,ABBA_LENGTH_IJKA)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
         Call Putall(Work, BBBB_LENGTH_IJKA, IRREPX, 108)
         Call Putall(Work, ABAB_LENGTH_IJKA, IRREPX, 109)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
      Else
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))
         LENGTH = ABBA_LENGTH_IJKA
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
      Endif

C Hbar(AB,CD)
      IF (IFLAGS(93) .EQ. 0) THEN
      IF (UHF) Then
         AAAA_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,231),ISYTYP(2,231))
         BBBB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,232),ISYTYP(2,232))
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
         LENGTH = MAX(AAAA_LENGTH_ABCD,BBBB_LENGTH_ABCD,
     +                ABAB_LENGTH_ABCD)
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, AAAA_LENGTH_ABCD, IRREPX, 231)
         Call Putall(Work, BBBB_LENGTH_ABCD, IRREPX, 232)
         Call Putall(Work, ABAB_LENGTH_ABCD, IRREPX, 233)
      Else
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
         LENGTH = ABAB_LENGTH_ABCD
         CALL DZERO(WORK,LENGTH)
         Call Putall(Work, ABAB_LENGTH_ABCD, IRREPX, 233)
      Endif
      Endif


      Return
      End

