










      subroutine initgrid(angfct)

c This initializes some constants used in the grids.

      implicit none







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






c ***NOTE*** This is a genuine (though not serious) limit on what Aces3 can do.
c     12 => s,p,d,f,g,h,i,j,k,l,m,n
      integer maxangshell
      parameter (maxangshell=12)




c This contains information about each of the possible grids for
c performing the numerical integration.

c###########################################################################
c MISC
c###########################################################################
c maxgrdatm : The largest atomic number for which the Slater and Bragg-Slater
c             atomic size is known.
c atmrad    : Atomic size using Slater's' rules for the radial integration
c xbsl      : The Bragg-Slater radii (one for each atom)
c             ***NOTE*** This is a genuine constraint.  Only atoms smaller
c             then this (currently 86) may be calculated.
c numangfct : number of angular momentum functions
c minpt     : The number of points used in the integration over
c             interatomic paths to find the minimum density point
c             between two atoms
c pangfct   : a pointer to the array of x,y,z angular momentum for each
c             angular momentum function

      integer maxgrdatm
      parameter (maxgrdatm=86)
      integer minpt
      parameter (minpt=100)

      integer numangfct

      double precision
     &    atmrad(maxgrdatm),xbsl(maxgrdatm),
     &     TA(maxgrdatm),multiEX(maxgrdatm)
      common /grid/  numangfct
      save /grid/

      common /gridd/ atmrad,xbsl,TA,multiEX
      save /gridd/

      integer pangfct
      common /gridp/ pangfct
      save /gridp/

c###########################################################################
c RADGRD file
c###########################################################################
c maxanggrid: The maximum number of different angular grids which can be
c             used in any given calculation.
c             ***NOTE*** This is a genuine constraint, but it must be used
c             since we must be able to keep a record of which angular grids
c             are used (before we have any allocated memory) since we have
c             to know how many grids are used in order to determine how much
c             memory to allocate.  This is set high enough it should never
c             be a problem.
c gridlist  : A list of all grids used (see the comment on maxanggrid).  It
c             is of dimension (maxanggrid,3) to keep track of the type and
c             subtype, and the number of times each grid is used.  The type
c             refers to how the grid is arrived at.
c                1 : Lebedev
c                2 : ?
c             The subtype refers to the degree of the grid of this type.
c numgrid   : The number of different angular grids used in the calculation.
c maxangpts : The maximum number of points in any of the angular grids used.
c maxanggrd : The grid with the maximum number of angular points.
c numradpts : The number of different radial points.
c ntotrad   : The total number of points in all angular grids at all radial
c             points (i.e. the entire integration grid)
c
c iradint   : determines if the Handy method (1) or Gauss-Legendre (2)
c             radial integration is used
c autosiz   : A flag which sets whether the polyhedra are (1) equally
c             sized, (2) sized according to Bragg-Slater radii or
c             (3) automatically sized according to the minimums in
c             density.
c slater    : A flag which determines whether (0) Slater's' rules are
c             used to determine the atomic size and scale the radial
c             integration or (1) no scaling is used.
c rigid     : A flag which determines whether rigid (0) or fuzzy
c             partitioning is used.
c nitr      : The number of iterations of the equations which create the
c             'fuzzy' boundary.

      integer maxanggrid
      parameter (maxanggrid=1000)

      integer gridlist(maxanggrid,3),numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad

      common /radgrd/  gridlist,numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad
      save /radgrd/

c Memory pointers
c
c pradgrid(numradpts) : The angular grid to use at each radial point.
c pgrdangpts(numgrid) : The number of angular points in each grid.
c zgridxyz(3,maxangpts,numgrid)
c                     : The x,y,z coordinate of each angular point in each grid
c zgridwt(maxangpts,numgrid)
c                     : The weight at each point.
c pintegaxis(natoms,3): Contains information about how much of each axis to
c                       integrate over.  If integaxis(iatom,i) is set to i,
c                       integrate only over the positive half of the i^th
c                       axis.  Otherwise, integrate over the entire axis.

      integer pgrdangpts,pradgrid,zgridxyz,zgridwt,pintegaxis

      common /radgrdp/ pgrdangpts,pradgrid,zgridxyz,zgridwt,
     &    pintegaxis
      save /radgrdp/

c###########################################################################
c Old stuff
c###########################################################################

c polist    : Contains an ordered list of unique atoms
c zatmvc    : The x, y, and z distance between each pair of atoms.
c zrij      : The distance between each pair of atoms.
c zatmpth   : The cartesian coordinates for the path integration between
c              all atom pairs
c zptdis    : The distance from atom i to a point along the path between
c              atoms i and j
c zprsqrd   : The distance squared from each atom to a point along all the
c              paths between all the atoms
c zpthpt    : The cartesian coordinates with respect to each atom for
c              the points along all the paths between all the atoms
c zbslrd    : The Bragg-Slater radii.
c zaij      : Surface shifting parameter dependent on the distance between
c               pairs of atoms.

      integer polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,zptdis,
     &    zprsqrd,zpthpt,zbslrd,zaij
      common /gridold/ polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,
     &    zptdis,zprsqrd,zpthpt,zbslrd,zaij
      save /gridold/


      integer angfct(numangfct,3)

      integer i,j,k,x,y,z,ptr

c Atomic size using Slater's rules for the radial integration
      data atmrad /1.000,0.588,
     &    3.077,2.051,1.538,1.231,1.026,0.879,0.769,0.684,
     &    4.091,3.158,2.571,2.169,1.875,1.651,1.475,1.333,
     &    6.27,4.84,4.59,4.38,4.20,4.01,3.82,3.68,3.53,
     &    1.192,3.27,3.16,2.76,2.44,2.19,1.98,1.81,1.66,
     &    7.29,5.67,5.37,5.12,4.91,4.69,4.46,4.31,4.12,
     &    3.97,3.82,3.68,3.23,2.85,2.55,2.31,2.12,1.94,
     &    8.03,6.24,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     &    0.00,0.00,0.00,0.00,0.00,0.00,0.00,5.90,5.60,
     &    5.40,5.16,4.91,4.74,4.54,4.37,4.20,4.06,3.55,
     &    3.14,2.80,0.00,0.00,0.00/

c Bragg-Slater radii for determining the relative size of the
c polyhedra in the polyatomic integration scheme
      data xbsl /0.529177,0.50,
     &    1.45,1.05,0.85,0.70,0.65,0.60,0.50,0.65,
     &    1.80,1.50,1.25,1.10,1.00,1.00,1.00,0.95,
     &    2.20,1.80,1.60,1.40,1.35,1.40,1.40,1.40,1.35,
     &    0.60,0.60,1.35,1.30,1.25,1.15,1.15,1.15,1.10,
     &    2.35,2.00,1.80,1.55,1.45,1.45,1.35,1.30,1.35,
     &    1.40,1.60,1.55,1.55,1.45,1.45,1.40,1.40,1.30,
     &    2.60,2.15,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     &    0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.75,1.55,
     &    1.45,1.35,1.35,1.30,1.35,1.35,1.35,1.50,1.90,
     &    1.75,1.60,0.00,0.00,0.00/

       data TA /0.8,
     &      0.9,1.8,1.4,1.3,1.1,0.9,0.9,0.9,0.9,
     &      1.4,1.3,1.3,1.2,1.1,1.0,1.0,1.0,
     &      1.5,1.4,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.1,
     &      1.1,1.1,1.1,1.0,0.9,0.9,0.9,0.9,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0/

       data multiEX/1.30,
     &      0.5882,1.95,2.20,1.45,1.20,1.10,1.10,1.20,
     &      0.6838,2.30,2.20,2.10,1.30,1.30,1.10,1.45,
     &      1.3333,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0/


      call callstack_push('INITGRID')

c Length of angular momentum function vectors
      if (numangfct.eq.0) then
        do i=1,maxangshell
          numangfct=numangfct+i*(i+1)/2
        end do
        goto 999
      end if

c The number of angular momentum functions
c angfct(*,1) = 0  100  211000  3221110000  ...
c angfct(*,2) = 0  010  010210  0102103210  ...
c angfct(*,3) = 0  001  001012  0010120123  ...
      ptr=1
      do i=1,maxangshell
        x=i-1
        do j=1,i
          y=j-1
          z=0
          do k=1,j
            angfct(ptr,1)=x
            angfct(ptr,2)=y
            angfct(ptr,3)=z
            y=y-1
            z=z+1
            ptr=ptr+1
          end do
          x=x-1
        end do
      end do

  999 continue
      call callstack_pop
      return
      end
