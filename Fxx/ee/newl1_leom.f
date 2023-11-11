






































































































































































































      SUBROUTINE NEWL1_LEOM(Work,Length,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Length)
      Logical Sing
      


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
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

       Sing   = .False.
       Sing   = (Iflags(2) .gt. 9)

       I000 = 1
       I010 = I000 + NT(1)
       I020 = I010 + NT(1)
       I030 = I020
       I040 = I020
       If (Iuhf .Gt. 0) Then
           I030 = I020 + NT(2)
           I040 = I030 + NT(2)
       Endif 

       If (I040 .ge. Length) Call Insmem("@-NEWL1_LEOM",I040,Length)

       If (Sing) Then
      
          Call Getlst(Work(I000), 1, 1, 1, 1, 93)
          Call Getlst(Work(I010), 1, 1, 1, 3, 90)
          Call Daxpy(NT(1), 1.0D0, Work(I000), 1, Work(I010), 1)
          Call Newl1(Work(I010), Work, Length, 1)
          Call Putlst(Work(I010), 1, 1, 1, 3, 90)

          If (Iuhf .ne. 0) Then

             Call Getlst(Work(I020), 1, 1, 1, 2, 93)
             Call Getlst(Work(I030), 1, 1, 1, 4, 90)
             Call Daxpy(NT(1), 1.0D0, Work(I020), 1, Work(I030), 1)
             Call Newl1(Work(I030), Work, Length, 1)
             Call Putlst(Work(I030), 1, 1, 1, 4, 90)

          Endif 
      Endif 
C
      Return
      End
