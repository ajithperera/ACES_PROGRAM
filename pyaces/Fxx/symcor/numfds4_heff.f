










      Subroutine numfds4_heff(icore, icrsiz)
C   
c loop over a variety of effective Hamiltonians and t-moments
c At this point we need to reserve space on jobarc for various quantities.
c Also we need to keep everything in core. Run through this twice.
C
      implicit none
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





      character*8 irrnmf(14), irrnmc(14) 
      integer icore(icrsiz)
      integer ioff, icase, nbasx, nstart, nirrep, irreps, itop, 
     &        dirprd, nrow, ndim, icrsiz, ione, ndim2, npoints
      double precision eref
C
      ioff = 1
      ione = 1
c        
      call getrec(20,'JOBARC','PARENERG',iintfp, icore(ioff))         
      ioff = ioff + iintfp
c        
c diabatic orbitals
c
      call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABOCC', 1, ndim)
      if (ndim .ne. 0) then
      call getrec(20, 'JOBARC', 'DIABCPI0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABVRT', 1, ndim)
      if (ndim .ne. 0) then
      call getrec(20, 'JOBARC', 'DIABCPA0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif
      
      do icase = 1, ncases
         ndim = 0
         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
         if (icase .eq. 1 .or. icase .eq. 2) then 
            nrow = nbasx
         else
            nrow = nbasx*nbasx
         endif
C
         if (ndim .ne. 0) then
            write(6,*) '@numfds4_heff ', icase, ndimheff(icase), ndim
c           
c process heff part
c     
            ndim2 = ndim*ndim
            call getrec(20,'JOBARC',nameheff(icase),
     $                  ndim2*iintfp, icore(ioff))
            ioff = ioff + ndim2*iintfp
c              
            call getrec(20,'JOBARC',namehdiab(icase),
     $                  nrow*ndim*iintfp, icore(ioff))
         write(6,*) ' namehdiab in heff, icase, ndim ',
     $           icase, ndim
         call output(icore(ioff), 1, ndim, 1, ndim, ndim,
     $               ndim, 1)
            ioff = ioff + nrow*ndim*iintfp
c
            call getrec(20,'JOBARC',nametran(icase),
     $                ndim2*iintfp, icore(ioff))
            ioff = ioff + ndim2*iintfp
c
c now process transition moments
c
            call getrec(20,'JOBARC',nametmom(icase),
     $                  ndim*3*iintfp, icore(ioff))
            ioff = ioff + 3*ndim*iintfp
         endif
      enddo
c
c do same for extra set of diabatization vectors
c
      do icase = 1, ncases
         ndim = 0
         call getrec(-1, 'JOBARC', namenvec(icase), ione, ndim)
         call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
         if (icase .eq. 1 .or. icase .eq. 2) then 
            nrow = nbasx
         else
            nrow = nbasx*nbasx
         endif
C
         if (ndim .ne. 0) then
            write(6,*) '@numfds4_heff ', icase, namenvec(icase), ndim
            write(6,*) ' extra vec ', namehvec(icase), ioff
            call getrec(20,'JOBARC',namehvec(icase),
     $                  nrow*ndim*iintfp, icore(ioff))
            ioff = ioff + nrow*ndim*iintfp
         endif
      enddo
c
c zero out JOBARC records: This is not necessary as far as I can
c tell. Revisit during Debuggiing. Ajith Perera
c
       itop = ioff
       call aces_ja_truncate('JOBARCST', 1)
c
c now set everything up for the calculation
c
      CALL GETREC(20,'JOBARC','NUMPOINT',IONE,npoints)
      write(6,*) ' npoints when initializing JOBARC', npoints
      itop = ioff
      ioff = 1
C
      call putrec(20,'JOBARC','REFPAR_E', iintfp, icore(ioff))
      call getrec(20,'JOBARC','REFPAR_E', iintfp, eref)
      write(6,*) ' @numfds4_heff, REFPAR_E to JOBARC ', eref
c
      ioff = ioff+iintfp
      call zero(icore(itop), npoints)
      call putrec(20,'JOBARC','PNTENERG', npoints*iintfp, icore(itop))
c
c        
c diabatic orbitals
c
      call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABOCC', 1, ndim)
      if (ndim .ne. 0) then
      call putrec(20, 'JOBARC', 'DIABCPI0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABVRT', 1, ndim)
      if (ndim .ne. 0) then
      call putrec(20, 'JOBARC', 'DIABCPA0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif

       do icase = 1, ncases
          ndim = 0
          call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
          call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
          if (icase .eq. 1 .or. icase .eq. 2) then
              nrow = nbasx
          else
              nrow = nbasx*nbasx
          endif
          if (ndim .ne. 0) then
c
c process heff part
c
             ndim2 = ndim*ndim
             call putrec(20,'JOBARC',refheff(icase),
     $                   ndim2*iintfp, icore(ioff))
             ioff = ioff + ndim2*iintfp
c
             call putrec(20,'JOBARC',namehdiab(icase),
     $                   nrow*ndim*iintfp, icore(ioff))
C
            write(6,*) ' namehdiab in heff putrec, nrow*ndim',
     $            nrow*ndim
C
             ioff = ioff + nrow*ndim*iintfp
c
             call putrec(20,'JOBARC',nametran(icase),
     $                   ndim2*iintfp, icore(ioff))
             ioff = ioff + ndim2*iintfp
c
c also create the pnt records
c
             call zero(icore(itop), ndim2*npoints)
             call putrec(20,'JOBARC',pntheff(icase),
     $                   ndim2*npoints*iintfp, icore(itop))
c
c now process transition moments
c
             call putrec(20,'JOBARC',reftmom(icase),
     $                    ndim*3*iintfp, icore(ioff))
c            call putrec(20,'JOBARC',nametmom(icase),
c    $                    ndim*3*iintfp, icore(ioff))
             ioff = ioff + ndim*3*iintfp
             call zero(icore(itop), ndim*3*npoints)
             call putrec(20,'JOBARC',pnttmom(icase),
     $                    ndim*3*npoints*iintfp, icore(itop))
c
          endif
      enddo
c
c same for extra vectors
c
      do icase = 1, ncases
         ndim = 0
         call getrec(-1, 'JOBARC', namenvec(icase), ione, ndim)
         call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
         if (icase .eq. 1 .or. icase .eq. 2) then 
            nrow = nbasx
         else
            nrow = nbasx*nbasx
         endif
C
         if (ndim .ne. 0) then
            write(6,*) ' extra vec 2 ', namehvec(icase), ioff
            call putrec(20,'JOBARC',namehvec(icase),
     $                  nrow*ndim*iintfp, icore(ioff))
            ioff = ioff + nrow*ndim*iintfp
         endif
      enddo         
c
      return
      end  

