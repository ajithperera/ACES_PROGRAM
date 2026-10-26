














































































































































































































C
C THIS PROGRAM DETERMINES SYMMETRY ADAPTED COORDINATES WHICH
C ARE USED FOR NUMERICAL DERIVATIVE CALCULATIONS.  IT
C CAN EASILY BE ADAPTED TO CONSTRUCT OTHER SYMMETRY ADAPTED
C QUANTITIES, SUCH AS BASIS FUNCTIONS.  IT CURRENTLY WORKS FOR
C ALL POINT GROUPS OTHER THAN THOSE WHICH CONTAIN COMPLEX
C REPRESENTATIONS (CN, SN, CNH (N>2) AND T AND TH).  FOR THESE
C "DANGEROUS GROUPS", THE ABELIAN SUBGROUP IS USED.
C WRITTEN BY J.F. STANTON, 1991-1992.

c The entire operating procedure of symcor is explained in README.

      SUBROUTINE SYMCOR(ICORE,ICRSIZ)
c
      implicit none
      INTEGER ICRSIZ
      INTEGER ICORE(ICRSIZ)
C
      LOGICAL          ENERONLY,GRADONLY,ROTPROJ,RAMAN,GMTRYOPT,
     &                 SNPTGRAD
      COMMON /CONTROL/ ENERONLY,GRADONLY,ROTPROJ,RAMAN,GMTRYOPT,
     &                 SNPTGRAD 
c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

      external aces_bd_parallel_aces




      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end







c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/





      CALL SYMCOR_ONEDIM(ICORE, ICRSIZ)

      RETURN
      END

