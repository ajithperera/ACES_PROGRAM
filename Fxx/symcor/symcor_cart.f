










      subroutine symcor_cart(refq, Rint, scr,nsize,nmode)
c
c Generate cartesian coordinates from a set of internals
c   generated in setpts4_pes_normal in symcor.
c This subroutine should be in symcor, but we use common blocks 
c that are used in joda exclusively. 
c The internals in Rint are overwrtten by cartesians
c The cartesians in Rint are aligned in maximum coincidence with refq
c
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
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


C
      double precision Refq(nsize), Rint(nsize), scr(*)
      logical print
c
      nrx = nsize
      nx6m = nmode
      print = .false.
c
      if (print) then
      write(6,*) '@symcor_cart: Rint on input, natoms', natoms
      call output(Rint, 1, 1, 1, nsize, 1, nsize, 1)
      endif
c
      call usqush(Rint,nx6m)
      call DCOPY(nrx, Rint, 1, R, 1)
      call gen_cart_coord_b(scr, print)
      call mn_align_geom(natoms, atmass, q, refq, scr, mxcor)
      call dcopy(nrx, q, 1, Rint, 1)
c
      if (print) then
      write(6,*) '@symcor_cart: Rint on output'
      call output(Rint, 1, 1, 1, nsize, 1, nsize, 1)
      endif

c
      return
      end
      
