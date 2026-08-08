










      Subroutine Check_t(Work,Maxcor,Iuhf)

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

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

C First T1(A,I)

      IRREPX =1 
      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)
      CALL Getlst(Work,1,1,1,1,90)
      Call Checksum("T1AI", Work, HPA_LENGTH,S)
      IF (UHF) Call Getlst(Work,1,1,1,2,90)
      IF (UHF) Call Checksum("T1ai", Work, HPB_LENGTH,S)

C T2(IJ,AB)
  
      Write(*,*) 
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 44)
         Call Checksum("T2AAAA_ABIJ", Work, AAAA_LENGTH_IJAB,S)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 45)
         Call Checksum("T2BBBB_abij", Work, BBBB_LENGTH_IJAB,S)
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 46)
         Call Checksum("T2ABAB_AbIj", Work, ABAB_LENGTH_IJAB,S)
      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 46)
         Call Checksum("T2ABAB_AbIj", Work, ABAB_LENGTH_IJAB,S)
      Endif

      Write(6,*)
      Return
      End


