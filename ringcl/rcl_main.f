














































































































































































































      program ringcl
      implicit none

      integer iuhf
      logical solve_4lambda,Density 

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
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      Solve_4lambda = .True.
      Density = (IFLAGS(19) .EQ. 1) 
      If (Density) Solve_4lambda = .True. 

C Respnse density always need left hand wavefunction. If only 
C excited states are requested then turn off lambda. 

      If((IFLAGS(87) .EQ.  1 .OR. 
     +   IFLAGS(87)  .EQ.  2 .OR. 
     +   IFLAGS(87)  .EQ.  4 .OR.
     +   IFLAGS(87)  .EQ. 13 .OR.
     +   IFLAGS(87)  .EQ. 14 .OR.
     +   IFLAGS(87)  .EQ. 15 .OR.
     +   IFLAGS(87)  .EQ. 16 .OR.
     +   IFLAGS(87)  .EQ. 17 .OR.
     +   IFLAGS(87)  .EQ. 18 .OR.
     +   IFLAGS(87)  .EQ. 19 .OR.
     +   IFLAGS(87)  .EQ. 20).AND. 
     +   (.NOT. Density)) Solve_4lambda = .False.

      Call Rcl_hbar(Icore(i0),icrsiz/iintfp,Iuhf,Solve_4lambda) 

      Call rcl_driver(Icore(i0),icrsiz/iintfp,Iuhf,Solve_4lambda)

CSSS        Call rcl_check_rcchbar(Icore(i0),icrsiz/iintfp,Iuhf) 

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

