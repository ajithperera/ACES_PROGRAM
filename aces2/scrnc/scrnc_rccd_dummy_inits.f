










      Subroutine Scrnc_rccd_dummy_inits(Work,Maxmem,Iuhf)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxmem)



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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

C Rccd and Drccd has no singles. Since all the routines that have
C been developed earier assumed CCSD (hence singles). Lets create
C singles list and set them to zero, so all of them still work
C correctly (cost consideration are minuscule). 

      Do Ispin = 1, Iuhf+1
         Call Updmoi(1,Nt(Ispin),Ispin,90,0,0)
         Call Updmoi(1,Nt(Ispin),Ispin,190,0,0)
         Call Updmoi(1,Nfmi(Ispin),Ispin,191,0,0)
         Call Updmoi(1,Nfea(Ispin),Ispin,192,0,0)

         Call Aces_list_memset(Ispin,90,0)
         Call Aces_list_memset(Ispin,190,0)
         Call Aces_list_memset(Ispin,191,0)
         Call Aces_list_memset(Ispin,192,0)
      Enddo 

      Return
      End 
