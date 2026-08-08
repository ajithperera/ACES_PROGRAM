










      Subroutine restore_pdcc_wmbej(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Integer AAAA_LENGTH_MBEJ,BBBB_LENGTH_MBEJ,AABB_LENGTH_MBEJ
      Integer BBAA_LENGTH_MBEJ,ABAB_LENGTH_MBEJ,BABA_LENGTH_MBEJ
      Double Precision Mone, One

      Dimension Work(Length)



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
C             -
C The P(D) CC W(mb,ej) is needed for the L2INL2 ring contributions.
C These are computed in post_cc_mods are stored in 165-169 list.
C Restore them on 54-59 lists.
C
      IRREPX = 1

      AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
      AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
      ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

      If (Iuhf .NE. 0) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
      End if



      IF (IUhf .GT. 0) Then

         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 165)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
         Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 167)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 169)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)

      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
         Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)

      Endif

      Return
      End
 
