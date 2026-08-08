










      subroutine gaussprak(flag,grid,type,grdangpts,gridxyz,gridwt,
     &         ntheta,nphi,glwght,glroot)

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



      integer grid,type,grdangpts(numgrid),flag

      double precision
     &    gridxyz(3,maxangpts,numgrid),gridwt(maxangpts,numgrid)

      integer maxlebgridP

      integer ntheta,nphi,itheta,iphi,ifact,ithea,itheb,iphia,iphib,
     &    i,nexppts(type)

      double precision
     &  pi,wthe,cost,sint,dphi,phii,wphi,cosp,sinp
 
      double precision
     &  glwght(ntheta,ntheta),glroot(ntheta,ntheta)

      double precision
     &    one,two,half
 
       parameter(one=1.0d0,two=2.0d0,half=0.5d0)

       call callstack_push('GAUSSPRAK')

      PI=ATAN(1.0D+00)*4.0D+00

      nexppts(type)=ntheta*nphi
      if(flag.eq.0) then
       ntotrad=ntotrad+nexppts(type)*gridlist(grid,3)
       if (nexppts(type).gt.maxangpts) then
          maxangpts=nexppts(type)
          maxanggrd=grid
      endif
        goto 999
      endif
      grdangpts(grid)=nexppts(type)

      do itheta=1,ntheta
      call findwt(-one,one,glroot(1,itheta),glwght(1,itheta),itheta)
      end do
  
      ntheta=max(ntheta,2)
      ifact=1
      ithea=1
      itheb=ntheta
      iphia=1
      iphib=nphi

      do itheta=ithea,itheb
        wthe=glwght(itheta,ntheta)
        cost=glroot(itheta,ntheta)
        sint=SQRT(one-cost**2)
        do iphi=iphia,iphib
           i=(itheta-ithea)*(iphib-iphia+1)+iphi
           dphi=two*pi/nphi
           phii=dphi*(iphi-half)
           wphi=dphi
           cosp=cos(phii)
           sinp=sin(phii)
           gridxyz(1,i,numgrid)=sint*cosp
           gridxyz(2,i,numgrid)=sint*sinp
           gridxyz(3,i,numgrid)=cost
           gridwt(i,numgrid)=wthe*wphi
        end do
      end do

 999  continue
      call callstack_pop
      return
      end

