










      subroutine print_internal(Rcurvy, nnx, nnx6m, nnatoms, iunit)
c
      implicit none
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end


C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end


      integer i, nnx6m, nnx, nnatoms, iunit, j, linblnk
      double precision Rcurvy(nnx), fact, bohr
c
      if (.false.) then
         write(6,*) ' Print_internal: Rcurvy on input '
         call output(Rcurvy, 1, 1, 1, nnx6m, 1, nnx6m, 1)
      endif
c
      call usqush(Rcurvy, nnx6m)
      fact = 180.0d0 / acos(-1.0d0)
      bohr = 1.88972594929722
c     
c     convert angles in radians to angles in degrees
c
      do i = 8, nnatoms*3-1, 3
         Rcurvy(i) = Rcurvy(i) * fact
      enddo
      do i = 12, nnatoms*3, 3
         Rcurvy(i) = Rcurvy(i) * fact
      enddo
c     
c      call sqush(Rcurvy, nnx)
c     
      write(iunit,880)
 880  format(' Internal ',4x, ' R (A, degrees)     ')
      do i = 1, nnx6m
         j = isquash(i)
         write(iunit,750) varnam(j)(1:linblnk(varnam(j))),
     $        Rcurvy(j)
      enddo
      write(iunit,*)
 750  format(5x, A10, '=', 2F16.8)
c     
c     convert R back to radians
c     
c      call usqush(Rcurvy, nnx6m)
c     
c     convert angles in radians to angles in degrees
c     
      do i = 8, nnatoms*3-1, 3
         Rcurvy(i) = Rcurvy(i) / fact
      enddo
      do i = 12, nnatoms*3, 3
         Rcurvy(i) = Rcurvy(i) / fact
      enddo
c     
      call sqush(Rcurvy, nnx)
c     
      return
      end
