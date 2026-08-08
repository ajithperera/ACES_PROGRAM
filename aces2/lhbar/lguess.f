













































































































































































































      Subroutine lguess(Work, Length, Iuhf, Sing)

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
    
      Irrepx = 1

      If (Iuhf .Ne. 0) Then
          NSIZE_AA = IDSYMSZ(IRREPX,ISYTYP(1,61),ISYTYP(2,61))
          NSIZE_BB = IDSYMSZ(IRREPX,ISYTYP(1,62),ISYTYP(2,62))
          NSIZE_AB = IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))
      Else
          NSIZE_AB = IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))
      Endif 

      I0SA=1
      IF(IUHF.NE.0)THEN
       I0SB =I0SA  + IRPDPD(IRREPX,9)
       I0SAA=I0SB  + IRPDPD(IRREPX,10)
       I0SBB=I0SAA + NSIZE_AA 
       I0SAB=I0SBB + NSIZE_BB 
       Iend =I0SAB + NSIZE_AB
      ELSE
       I0SAB=I0SA  + IRPDPD(IRREPX,10)
       Iend =I0SAB + NSIZE_AB
      ENDIF

      Call Zero(Work(I0SA),IRPDPD(IRREPX,9))
      If (Iuhf .ne.0) Call Zero(Work(I0SB),IRPDPD(IRREPX,9))

      If (Iuhf .ne. 0) Then
         If (SING) then
            Call Putlst(Work(I0SA),1,1,1,3,90)
            Call Putlst(Work(I0SB),1,1,1,4,90)
         Endif 
         Call Getall(Work(I0SAA),NSIZE_AA,1,14)
         Call Putall(Work(I0SAA),1,Irrepx,61)
         Call Getall(Work(I0SBB),NSIZE_AA,1,15)
         Call Putall(Work(I0SBB),1,Irrepx,62)
         Call Getall(Work(I0SAB),NSIZE_AA,1,16)
         Call Putall(Work(I0SAB),1,Irrepx,63)
      Else
         If (sing) Call Putlst(Work(I0SA),1,1,1,3,90)
         Call Getall(Work(I0SAB),NSIZE_AB,1,16)
         Call Putall(Work(I0SAB),1,Irrepx,63)
      Endif 

      Return
      End

      
