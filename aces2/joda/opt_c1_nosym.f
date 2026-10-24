














































































































































































































      SUBROUTINE OPT_C1_NOSYM(SCRATCH, QTMP, NEWQ, NOSILENT)
C
C Cartesian and RIC optimizations run the electronic structure in the
C point group of the starting structure. When that is C1 there is
C nothing to gain, and a structure that becomes symmetric during the
C optimization would stop the run (SYMMETRY: ascent in symmetry). So 
C for a C1 start switch to SYMMETRY=NONE (and NOREORI=ON unless it was
C set), as these optimizations always did before 10/2026. The JOBARC 
C record OPTSYMOF tells GTFLGS to do the same in the later joda runs.
C Called right after the first SYMMETRY call in FETCHZ.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL NOSILENT
      CHARACTER*4 FPGRP, BPGRP, PGRP
      COMMON /PtGp_com/ FPGrp, BPGrp, PGrp
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      DIMENSION SCRATCH(*), QTMP(*), NEWQ(*)
      DOUBLE PRECISION NEWQ
C
      IF (IFLAGS(60) .EQ. 0) RETURN
      IF (FPGRP(1:2) .NE. 'C1') RETURN
C
      WRITE(6,"(a,a)") "  The starting structure has no symmetry (C1);",
     &                 " the optimization runs with SYMMETRY=NONE."
      CALL PUTREC(20,'JOBARC','OPTSYMOF',1,1)
      IFLAGS(60) = 0
      IF (IFLAGS2(125) .EQ. 2) 
     &    IFLAGS2(125) = 1
      CALL PUTREC(20,'JOBARC','IFLAGS  ',100,IFLAGS)
      CALL PUTREC(20,'JOBARC','IFLAGS2 ',500,IFLAGS2)
C
      CALL SYMMETRY(SCRATCH, QTMP, NEWQ, NOSILENT)
C
      RETURN
      END
