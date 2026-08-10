














































































































































































































      SUBROUTINE HBAR_(ICORE,ICRSIZ,IUHF)
      implicit none

      integer icore(icrsiz),icrsiz,iuhf,i0

c COMMON BLOCKS
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




c ----------------------------------------------------------------------
      I0 = 1
C
      Call Drive_Hbar(Icore(i0), icrsiz/iintfp, Iuhf)

C
c ----------------------------------------------------------------------
      RETURN
      end


C ----------------------------------------------------------------------
C Thin PROGRAM wrapper restoring the standalone xHBAR entry point
C (MAIN__), lost when this module was converted to a SUBROUTINE for
C pyaces (see project memory: ACES_PROGRAM unification, 2026-08). This
C calls the shared aces_init/aces_fin pair so the standalone binary
C manages its own memory exactly as the classic driver's own separate
C process used to, before the pyaces conversion removed that from
C inside the subroutine itself.
C ----------------------------------------------------------------------
      PROGRAM XHBAR
      IMPLICIT NONE


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
      INTEGER IUHF
      CALL ACES_INIT(ICORE, I0, ICRSIZ, IUHF, .TRUE.)
      CALL HBAR_(ICORE(I0), ICRSIZ, IUHF)
      CALL ACES_FIN
      STOP
      END
