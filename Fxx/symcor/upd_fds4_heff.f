






































































































































































































      subroutine upd_fds4_heff(NATOM, DOIT, IMORE, FD_POINTS, 
     &                         IPTTYPE, DSCR, NPOINT, NDSCR)
c     
c this routine determines the next point which must be run in a
c numerical frequency calculation, and the calculation type
c (single point energy or gradient), by inspecting information
c on the jobarc file.  
c     
      implicit double precision(a-h,o-z)
C
      integer ncases
      parameter (ncases=10)
      character*8 ndimheff(ncases)
      character*8 nameheff(ncases)
      character*8 namenvec(ncases)
      character*8 namehvec(ncases)
      character*8 namehdiab(ncases)
      character*8 nametran(ncases)
      character*8 nametmom(ncases)
      character*8 refheff(ncases)
      character*8 reftmom(ncases)
      character*8 pntheff(ncases)
      character*8 pnttmom(ncases)
c
      data ndimheff /'NEFF_IP2', 'NEFF_EA2', 'NEFF_EE1', 'NEFF_EE3',
     $                'NEFFDIP1', 'NEFFDIP3', 'NEFFDEA1', 'NEFFDEA3',
     $                'NEFFDEE1', 'NEFFDEE3' /
      data nameheff /'HEFF_IP2', 'HEFF_EA2', 'HEFF_EE1', 'HEFF_EE3',
     $                'HEFFDIP1', 'HEFFDIP3', 'HEFFDEA1', 'HEFFDEA3',
     $                'HEFFDEE1', 'HEFFDEE3' /
      data namenvec /'NVEC_IP2', 'NVEC_EA2', 'NVEC_EE1', 'NVEC_EE3',
     $     'NVECDIP1', 'NVECIP3', 'NVECDEA1', 'NVECDEA3',
     $     'NVECDEE1', 'NVECDEE3' /
      data namehvec /'HVEC_IP2', 'HVEC_EA2', 'HVEC_EE1', 'HVEC_EE3',
     $     'HVECDIP1', 'HVECDIP3', 'HVECDEA1', 'HVECDEA3',
     $     'HVECDEE1', 'HVECDEE3' /
c     
      data namehdiab /'HDIA_IP2', 'HDIA_EA2', 'HDIA_EE1', 'HDIA_EE3',
     $     'HDIADIP1', 'HDIADIP3', 'HDIADEA1', 'HDIADEA3',
     $     'HDIADEE1', 'HDIADEE3' /
      data nametran /'TRAN_IP2', 'TRAN_EA2', 'TRAN_EE1', 'TRAN_EE3',
     $     'TRANDIP1', 'TRANDIP3', 'TRANDEA1', 'TRANDEA3',
     $     'TRANDEE1', 'TRANDEE3' /
      data refheff /'HREF_IP2', 'HREF_EA2', 'HREF_EE1', 'HREF_EE3',
     $                'HREFDIP1', 'HREFDIP3', 'HREFDEA1', 'HREFDEA3',
     $                'HREFDEE1', 'HREFDEE3' /
      data pntheff /'HPNT_IP2', 'HPNT_EA2', 'HPNT_EE1', 'HPNT_EE3',
     $                'HPNTDIP1', 'HPNTDIP3', 'HPNTDEA1', 'HPNTDEA3',
     $                'HPNTDEE1', 'HPNTDEE3' /
      data nametmom /'TMOM_IP2', 'TMOM_EA2', 'TMOM_EE1', 'TMOM_EE3',
     $                'TMOMDIP1', 'TMOMDIP3', 'TMOMDEA1', 'TMOMDEA3',
     $                'TMOMDEE1', 'TMOMDEE3' /
      data reftmom /'TREF_IP2', 'TREF_EA2', 'TREF_EE1', 'TREF_EE3',
     $                'TREFDIP1', 'TREFDIP3', 'TREFDEA1', 'TREFDEA3',
     $                'TREFDEE1', 'TREFDEE3' /
      data pnttmom /'TPNT_IP2', 'TPNT_EA2', 'TPNT_EE1', 'TPNT_EE3',
     $                'TPNTDIP1', 'TPNTDIP3', 'TPNTDEA1', 'TPNTDEA3',
     $                'TPNTDEE1', 'TPNTDEE3' /
c
c Effective hamiltonians could be stored on JOBARC for a variety of methods.
c IP / EA / EE / DIP / DEA /DEE indicate the IP/EA-eomcc methods or the STEOM variants. DEE is eomee
c The last digit 1, 2 or 3 indicates the multiplicity of the manifold.
c
c The mrcc program writes the # of states included (ndimheff), and it writes the
c records 'nameheff' (ndim*ndim) dimensional and a transition moment record nametmom,
c of dimension 3*ndim. they contain the effective hamiltonian and the tmoms in 
c the diabatic basis.
c
c The rest of these records is internal to symcor. The reference values (refheff, reftmom)
c need to be
c stored and the finite difference records for all geometry points (pntheff, pnttmom)
c
c finally we will calculate heffgrd, tmomgrd, and possibly hefffcm and tmomfcm later on.
c These will be written to files HGRD_method and TGRD_method etc.
c





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




      LOGICAL          ENERONLY,GRADONLY,ROTPROJ,RAMAN,GMTRYOPT
      COMMON /CONTROL/ ENERONLY,GRADONLY,ROTPROJ,RAMAN,GMTRYOPT
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
 
      character*4 doit
C
      double precision  DSCR(NDSCR), FD_POINTS(3*NATOM*NPOINT)
      logical print
C
      integer ipttype(NPOINT)
c     
      data tol   /1.d-8/
      data one   /1.0/
      data zilch /0.0/
c     
      ione=1
      print = .true.
      nsize=3*natom
c   o read point type record and the list of displacements
      write(6,*) ' @upd_fds4_heff: npoint : ', npoint
      CALL GETREC(20,'JOBARC','FDCALCTP',NPOINT,IPTTYPE)
      if (print) then
      write(6,*) ' @upd_fds4_heff, ipttype'
      do i = 1, npoint
         write(6,*) i, ipttype(i)
      enddo
      endif

c Point types can be one of the following:
c  > 0 : this point will be done (1+iRank will be done by this process)
c  = 0 : this point is skipped
c  < 0 : this point is done (-1-iRank was done by this process)
c
c   o find first entry which must be calculated
      inext = 1
      do while ((ipttype(inext).ne.1+irank).and.(inext.le.npoint))
         inext = inext + 1
      end do

c   o find last entry which was calculated
      ilast = inext-1
      do while ((ipttype(ilast).ne.-1-irank).and.(ilast.gt.0))
         ilast = ilast - 1
      end do

      write(6,*) ' @upd_fds4_heff, inext  : ', inext
      if (print) then
      write(6,*) ' @upd_fds4_heff, ilast  : ', ilast
      write(6,*) ' @upd_fds4_heff, npoint : ', npoint
      endif
c
c   o prepare info for the next calculation
      IF (INEXT.NE.NPOINT+1) THEN
         CALL GETREC(20,'JOBARC','FDCOORDS',IINTFP*NSIZE*INEXT,
     &               FD_POINTS)
         ILOC=1+(INEXT-1)*NSIZE
         CALL PUTREC(20,'JOBARC','NEXTGEOM',IINTFP*NSIZE,
     &               FD_POINTS(ILOC))
         if (print) then
            write(6,*) ' NEXTGEOM in upd_fds4_heff'
            call output(FD_POINTS(iloc), 1, 1, 1, nsize, 1,
     $           nsize, 1)
         endif
         IPTTYPE(INEXT)=-IPTTYPE(INEXT)
         CALL PUTREC(20,'JOBARC','FDCALCTP',INEXT,IPTTYPE)
         IMORE=1
c      o tag the last displacement
         inext = inext+1
         do while ((ipttype(inext).ne.1+irank).and.(inext.le.npoint))
            inext = inext + 1
         end do
         if (inext.eq.npoint+1) then
            CALL PUTREC(20,'JOBARC','LASTGEOM',1,1)
         end if
      ELSE
         IMORE = 0 
      END IF

      if (Ilast .NE. 0) Then
c
c if ACES2 calculation:
c
         write(6,*) ' @upd_fds4_heff: Program flag ',
     $        IFLAGS2(132)
         IF (IFLAGS2(132).NE.3) THEN
c      o update energy vector
         call getrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)
         IF (IFLAGS(87).EQ.0) THEN
            CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,DSCR(iloc))
         ELSE
            CALL GETREC(20,'JOBARC','TOTENER2',IINTFP,DSCR(iloc))
         END IF
         call putrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
         write(6,*) ' @upd_fds4_heff: update ACES2 PNTENERG '
         call output(dscr, 1, 1, 1, npoint, 1, npoint, 1)
         go to 999
         endif

c    
c     process information from last calculation
c     
c process the parent energy
c
         call getrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)
         call getrec(20, 'JOBARC', 'PARENERG',
     $        iintfp, DSCR(iloc))
         call putrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
c     
c loop over a variety of effective Hamiltonians and t-moments
c     
         do icase = 1, ncases
            ndim = 0
            call getrec(20, 'JOBARC', ndimheff(icase), ione, ndim)

            if (ndim .ne. 0) then
c     
c process heff part
c
               write(6,*) ' Process icase ', icase, ndimheff(icase),
     $              ndim
               write(6,*) ' Process icase ', icase, ndimheff(icase)
               ndim2 = ndim*ndim
               call getrec(20,'JOBARC',pntheff(icase),
     $              ndim2*npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)*ndim2
               call getrec(20,'JOBARC',nameheff(icase),
     $              ndim2*iintfp, DSCR(iloc))
               write(6,880)
     $              nameheff(icase), ilast, iloc
 880           format(' Read Heff from Jobarc ',A12, 2I10)
               call output(DSCR(iloc), 1, ndim, 1, ndim, ndim,
     $              ndim, 1)
               call putrec(20,'JOBARC',pntheff(icase),
     $              ndim2*npoint*iintfp, DSCR) 
c     
c now process transition moments
c     
               call getrec(20,'JOBARC',pnttmom(icase),
     $              ndim*3*npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)*ndim*3
               call getrec(20,'JOBARC',nametmom(icase),
     $              ndim*3*iintfp, DSCR(iloc))
               write(6,881) pnttmom(icase), ilast, iloc
 881           format(' Read Tmom from Jobarc ',A12, 2I10)
               call output(DSCR(iloc), 1, 3, 1, ndim, 3,
     $              ndim, 1)
               call putrec(20,'JOBARC',pnttmom(icase),
     $              ndim*3*npoint*iintfp, DSCR)
c     
            endif
         enddo
      endif
c     
c markzero indicates the position from where Jobarc can be zeroed out
c
 999  continue
      CALL PUTREC(20,'JOBARC','MARKZERO',ione, ione)
c
      return
      end
