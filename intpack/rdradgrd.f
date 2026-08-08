










      subroutine rdradgrd(radgrid)

c This will determine which grid to use for each radial point.  It will
c be called twice.  The first time (determined by numgrid=0 at the start
c of the routine) it reads constants and determine numgrid.  The second
c time (numgrid>0) it actually fills in the grids.

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




c TODO: needs proper gfname lookup
      integer radgio
      character*(*) radgfil
      parameter (radgio=40)
      parameter (radgfil='RADGRD')


c This contains flags that are set in the INTGRT namelist.  See the
c file initintgrt.F for a description of each of them.
c
c The following are exceptions:
c    int_ks           : .true. if we are doing Kohn-Sham
c    int_ks_finaliter : .true. if this is the final iteration of KS
c    int_ks_exch      : which potential to use to calculate exchange
c    int_ks_corr      : which potential to use to calculate correlation
c    int_kspot        : which hybrid functional to use to calculate potential
c                       (if equal to fun_special, use int_ks_exch and
c                       int_ks_corr)
c    int_dft_fun      : which functional to use with any SCF density
c                       Added and modified by Stan Ivanov
c                       (if fun_special the functional is user-defined
c                        if fun_hyb_name then use hybrid functional) 
c    int_printlev     : 0 if we are doing a dft calculation, 1 if we
c                       are doing the final iteration of a KS calculation,
c                       2 if we are doing a KS iteration.
c These are set in the calling routines, NOT in the namelist.
c                     : Additions by S. Ivanov
c     num_acc_ks      : .true.  if numerical accelerator is used for KS 
c                        Default is .true.
c     ks_exact_ex     : .true. if exact LOCAL exchange is used for KS
c                        Deafult is .false.
c     int_tdks        : .true. if time-dependent KS calculation is
c                        requested
c                        Default is .false.
c     int_ks_scf      : .true. if the actual KS SCF energy is being
c                        calculated and printed out. Default is .false.
c
      integer int_numradpts,int_radtyp,int_partpoly,int_radscal,
     &    int_parttyp,int_fuzzyiter,int_defenegrid,int_defenetype,
     &    int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp

      logical int_ks,num_acc_ks,ks_exact_ex,int_tdks,int_ks_scf,
     &        int_ks_finaliter

      double precision
     &    int_radlimit,coef_pot_nonlocal

      common /intgrtflags/  int_numradpts,int_radtyp,int_partpoly,
     &    int_radscal,int_parttyp,int_fuzzyiter,int_defenegrid,
     &    int_defenetype,int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_ks_finaliter,int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp
c

c  prakash added int_ksmem to the common block
      common /intgrtflagsd/ int_radlimit,coef_pot_nonlocal
      common /intgrtflagsl/ int_ks,num_acc_ks,ks_exact_ex,int_tdks,
     &                      int_ks_scf

      save /intgrtflags/
      save /intgrtflagsl/
      save /intgrtflagsd/

c The following are parameters used in the namelist

      integer int_prt_never,int_prt_dft,int_prt_ks,int_prt_always
      parameter (int_prt_never   =1)
      parameter (int_prt_dft     =2)
      parameter (int_prt_ks      =3)
      parameter (int_prt_always  =4)

      integer int_radtyp_handy,int_radtyp_gl
      parameter (int_radtyp_handy=1)
      parameter (int_radtyp_gl   =2)

      integer int_partpoly_equal,int_partpoly_bsrad,
     &    int_partpoly_dynamic
      parameter (int_partpoly_equal  =1)
      parameter (int_partpoly_bsrad  =2)
      parameter (int_partpoly_dynamic=3)

      integer int_radscal_none,int_radscal_slater
      parameter (int_radscal_none  =1)
      parameter (int_radscal_slater=2)

      integer int_parttyp_rigid,int_parttyp_fuzzy
      parameter (int_parttyp_rigid=1)
      parameter (int_parttyp_fuzzy=2)

      integer int_gridtype_leb
      parameter (int_gridtype_leb=1)




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


      integer radgrid(int_numradpts)

      integer i,radpt,type,subtype,err,j,flag
      logical exists
      character*80 line

      call callstack_push('RDRADGRD')

      flag=numgrid

c If a file containing specific radial grids (RADGRD) exists, use those.

      inquire(file=radgfil,exist=exists)
      if (exists) then
        open(unit=radgio,file=radgfil,form='formatted')
        rewind(radgio)
        err=0

c       Read the angular grid associated with each radial point.  The first
c       time this is called (numgrid.eq.0), calculate numgrid.

        read (radgio,'(a)') line
        do 30 i=1,int_numradpts
          read (radgio,9000) radpt,type,subtype
 9000     format(i4,i2,i2)

          if (flag.eq.0) then
            numgrid=0
            do 10 j=1,numgrid
              if (type.eq.gridlist(j,1) .and.
     &            subtype.eq.gridlist(j,2)) then
                gridlist(j,3)=gridlist(j,3)+1
                goto 30
              endif
   10       continue
            numgrid=numgrid+1
            gridlist(numgrid,1)=type
            gridlist(numgrid,2)=subtype
            gridlist(numgrid,3)=1
          else
            do 20 j=1,numgrid
              if (type.eq.gridlist(j,1) .and.
     &            subtype.eq.gridlist(j,2)) then
                radgrid(i)=j
                goto 30
              end if
   20       continue
          end if
   30   continue
        close(unit=radgio)

c If the file does not exist, use a default way of assigning grids.

      else
        if (flag.eq.0) then
          numgrid=1

          if (int_ks) then
            gridlist(1,1)=int_defpottype
            gridlist(1,2)=int_defpotgrid
          else
            gridlist(1,1)=int_defenetype
            gridlist(1,2)=int_defenegrid
         end if

         gridlist(1,3)=int_numradpts
        else
          do 40 i=1,int_numradpts
            radgrid(i)=1
   40     continue
        endif
      endif
      call callstack_pop
      return

  900 write (6,9010)
 9010 format(t3,'@RDRADGRD - error in RADGRD file')
      call errex
      return
      end
