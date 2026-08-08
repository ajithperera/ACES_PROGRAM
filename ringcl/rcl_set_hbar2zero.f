










      Subroutine Rcl_set_hbar2zero(Work,Maxcor,Iuhf)

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

C This routine selectively set Hbar elements to zero.

C Hbar(MN,IJ)
      IF (UHF) Then
         Call Zerolist(Work,Maxcor,51)
         Call Zerolist(Work,Maxcor,52)
         Call Zerolist(Work,Maxcor,53)
      Else
         Call Zerolist(Work,Maxcor,53)
      Endif

C These are needed from doubles correction.
C Hbar(AI,BC)
CSSS      IF (UHF) Then
CSSS         Call Zerolist(Work,Maxcor,27)
CSSS         Call Zerolist(Work,Maxcor,28)
CSSS         Call Zerolist(Work,Maxcor,29)
CSSS         Call Zerolist(Work,Maxcor,30)
CSSS      Else
CSSS         Call Zerolist(Work,Maxcor,27)
CSSS         Call Zerolist(Work,Maxcor,30)
CSSS      Endif 

C Hbar(IJ,KA)
CSSS      IF (UHF) Then
CSSS         Call Zerolist(Work,Maxcor,7)
CSSS         Call Zerolist(Work,Maxcor,8)
CSSS         Call Zerolist(Work,Maxcor,9)
CSSS         Call Zerolist(Work,Maxcor,10)
CSSS      Else
CSSS         Call Zerolist(Work,Maxcor,7)
CSSS         Call Zerolist(Work,Maxcor,10)
CSSS      Endif

C Hbar(AB,CI)
      IF (UHF) Then
         Call Zerolist(Work,Maxcor,127)
         Call Zerolist(Work,Maxcor,128)
         Call Zerolist(Work,Maxcor,129)
         Call Zerolist(Work,Maxcor,130)
      Else
         Call Zerolist(Work,Maxcor,130)
      Endif

C Hbar(IA,JK)
      IF (UHF) Then
         Call Zerolist(Work,Maxcor,107)
         Call Zerolist(Work,Maxcor,108)
         Call Zerolist(Work,Maxcor,109)
         Call Zerolist(Work,Maxcor,110)
      Else
CSSS         Call Zerolist(Work,Maxcor,107)
         Call Zerolist(Work,Maxcor,110)
      Endif

      Return
      End

