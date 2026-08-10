














































































































































































































      SUBROUTINE PROPS_(IUHF)
      Implicit none
      Integer iuhf

c COMMON BLOCKS


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end


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



       
c ----------------------------------------------------------------------
      Call Aces_init(icore,i0,icrsiz,iuhf,.true.)

      Call Props(icore(i0),Icrsiz,iuhf)

c ----------------------------------------------------------------------

      RETURN
      End


C ----------------------------------------------------------------------
C Thin PROGRAM wrapper restoring the standalone xprops entry point
C (MAIN__). PROPS_ is self-contained (calls Aces_init itself), so no
C aces_init call is needed here, only aces_fin to close out cleanly.
C ----------------------------------------------------------------------
      PROGRAM XPROPS
      IMPLICIT NONE
      INTEGER IUHF
      CALL PROPS_(IUHF)
      CALL ACES_FIN
      STOP
      END
