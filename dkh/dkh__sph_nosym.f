













































































































































































































      Subroutine Dkh__sph_nosym(Work,Mxdcor,Iuhf,Spherical,Cartesian,
     +                          Symmetry,Contract)

      Implicit Double Precision (A-H, O-Z)



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
      parameter (max_centers = 300)
      parameter (max_shells  = 5000)
      parameter (max_prims   = max_shells)
      Parameter (max_cbf     = 1000)
      parameter (max_primcc  = Max_prims*max_cbf)
      parameter (Ndi4 = 550, Ndi9 = max_shells, Ndi10 = max_centers,
     &           Ndi13 = 350, Ndico = 10, ndi14 = 120, ndi27 = 400)
     &
      parameter (Maxang = 6)
      parameter (Maxproj =5)
      parameter (ndilmx = Maxang+1)
      parameter (nh4=4*(ndilmx-1)-3)
      parameter (maxjco = 10)
      Parameter(Max_ecpmem = 50000)
C
C


cvw   ===================================================================================
cvw
cvw   basic parameters:
cvw   maximum angular momentum of basis functions     (mxang)
cvw   maximum angular momentum of ECP projector + 1   (mxproj)
cvw
cvw   currently, up to i functions in basis   
cvw              up to h projectors          (note that this means 'lmax=5')
cvw
cvw   ===================================================================================
cvw      
      integer   mxang, mxproj
      parameter (mxang =6)
C if mxproj is to be changed, do remember to change the same variable in
C ecppar and ecpcab !!
      parameter (mxproj=5)
cvw      
cvw   ===================================================================================
cvw
cvw   derived parameters
cvw
cvw   lmaxcomb: maximum 'combined' angular momentum
cvw             either from basis function and ECP projector
cvw             or from the product of two basis functions.
cvw             This is the maximum angular momentum for which
cvw             angular momentum gymnastics have to be performed
      integer lmaxcomb
      parameter(lmaxcomb=mxang+max(mxang,mxproj-1))
cvw
cvw   nftmax:   maximum number of cartesian mononomials
cvw             (i. e. the dimension of the n_cart, l_cart, m_cart arrays)
cvw   nftmxg:   the same for gradient calculations
cvw
      integer nftmax,nftmxg
cvw
cvw   ndegen(l) = (l+1)(l+2)/2, number of cartesian components for ang. mom. l
cvw   nftmax = Sum{l=0 ... mxang} ndegen(l)
cvw   nftmxg = Sum{l=0 ... mxang} ndegen(l-1)+ndegen(l+1)
cvw   nftmxe = Sum{l=0 ... mxang} ndegen(l+2)
cvw
cvw   nftmxg primitives arise in geo gradient
cvw   nftmxe primitives arise in basis gradient (exponent derivatives)
cvw
      parameter (nftmax=(mxang+1)*(mxang+2)*(mxang+3)/6)
      parameter (nftmxg=(mxang+1)*(9+5*mxang+mxang*mxang)/3)
      parameter (nftmxe=(mxang+1)*(36+11*mxang+mxang*mxang)/6)
cvw 
cvw   ===================================================================================
cvw 
cvw   TABLES FOR ANGULAR MOMENTUM GYMNASTICS:
cvw
cvw   angular momentum data is indexed by l, m
cvw   l=1,2,3,.... ,lfdim
cvw   m = 1,2,3,... 2*l-1
cvw
cvw   and stored in a linear array at position lf(l)+m
cvw   linear dimension is lmfdim = lfdim**2
cvw
cvw   Real spherical harmonics Y(l,m) are expressed as a linear combination
cvw   of cartesian products ZLM(i,j,k) x**i y**j z**k
cvw
cvw   The coefficients are in zlm(istart ... iend)
cvw   with istart=lmf(lf(l)+m), iend=lml(lf(l)+,m)
cvw   and the corresponding values of i,j,k in the arrays lmx, lmy, lmz
cvw
cvw   Note that for a given l,m, all cartesian terms have the same
cvw   parity pattern of the ijk. This pattern is stored in lmm(lf(l)+m)
cvw
cvw   bit 0 of lmm is set if k is EVEN
cvw   bit 1 of lmm is set if j is EVEN
cvw   bit 2 of lmm is set if i is EVEN
cvw
cvw   The total length of the zlm, lmx, lmy, lmz arrays is lmxdim, a number
cvw   which must be calculated by a rather complex formula
cvw
cvw   There are some auxiliary tables such as double factorial (dfac),
cvw   its inverse (dfaci) and binomial coefficients (binom_coef)
cvw   which are also stored in this data structure.
cvw
      integer lfdim,lmfdim,lmxdim
      parameter(lfdim=lmaxcomb+1)
      parameter(lmfdim=lfdim**2)
      parameter(lmxdim=(
     &    lmaxcomb*(lmaxcomb+2)*(lmaxcomb+3)*(lmaxcomb+4)/3 +
     &    (lmaxcomb+2)**2 * (lmaxcomb+4)
     &    ) /16)
cvw      
cvw   ===================================================================================
cvw
cvw   cvw_ecp_r/cvw_ecp_i: readonly data, loaded by subr. tab_ecp
cvw                        note that flmtx, mc, mr, mrclo, mrchi
cvw                        are only needed for spin-orbit integrals
cvw                        (compressed matrices of Lx, Ly, Lz between
cvw                         real spherical harmonics of same l)
cvw
cvw   common cvw_ecp_i contains the integer/logical data
cvw   common cvw_ecp_r contains the floating point data
cvw
      double precision zlm(lmxdim)
      double precision flmtx(mxproj*mxproj,3)
      double precision binom_coef(lfdim*(lfdim+1)/2)
      integer          lf(lfdim),lmf(lmfdim),lml(lmfdim)
      integer          lmx(lmxdim),lmy(lmxdim),lmz(lmxdim)
      integer          lmm(lmfdim)
      integer          n_cart(nftmax),l_cart(nftmax),m_cart(nftmax)
      integer          mc(mxproj*mxproj,3),mr(mxproj*mxproj,3)
      integer          mrclo(mxproj),mrchi(mxproj)
      double precision dfac (4*mxang+2*mxproj+9)
      double precision dfaci(4*mxang+2*mxproj+9)
      double precision fprod(lfdim,lfdim)
cvw      
      integer          nprim_shell(mxang+1)
      integer          nprim_shellg(mxang+1)
      integer          nlm_grd(nftmxg,7)
      integer          nlm_min(mxang+1)
      integer          nlm_max(mxang+1)
      integer          nlm_ming(mxang+1)
      integer          nlm_maxg(mxang+1)
cvw
cvw   mrclo, mrchi, mc, mr:    only used for SO integrals
cvw   nlm_grd:                 only used in gradient code
cvw
      common /cvw_ecp_r/ zlm,binom_coef,flmtx,fprod,
     &                   dfac,dfaci
      common /cvw_ecp_i/ lf,lmf,lml,lmx,lmy,lmz,lmm,
     &                   n_cart,l_cart,m_cart,
     &                   mrclo, mrchi, mc, mr,
     &                   nlm_grd,nlm_min,nlm_max,
     &                   nlm_ming,nlm_maxg,
     &                   nprim_shell,nprim_shellg
cvw      
cvw   ===================================================================================
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




      Parameter (Maxbfns = 1000,Mxicor= 1000000)
      Character*10 INP_Fname
      Character*4 ATMNAM(Max_centers), Namat(Max_centers)
      Logical Spherical, Cartesian, Sym_adapt 
      Logical Contract,Symmetry
      Integer Nelement,Dkh_order
    
      Integer CCbeg(Maxbfns), CCend(Maxbfns), Zmax, Zhrr, 
     &        End_nfps, Val1(4), Val2(4), Angmom, Atm_4shell,
     &        Erd_index(Maxbfns)

      DIMENSION WORK(Mxdcor), Iwork(Mxicor), Erd_scale(Maxbfns)

      DIMENSION NFCT(max_centers), NAOATM(max_centers),
     &          NUFCT(max_centers), NAOUATM(max_centers),
     &          iqmstr(max_centers),jcostr(max_centers,Maxjco),
     &          nucstr(max_centers,Maxang,Maxjco),
     &          nrcstr(max_centers,Maxang,Maxjco),
     &          Jstrt(max_shells),Jrs(max_shells),
     &          NCFpS(max_shells),Angmom(max_shells),
     &          NPFpS(max_shells), Atomic_label(max_shells),
     &          Charge(max_centers),Nspc(max_centers), 
     &          NPOP(max_centers),
     &          IREORDER(max_centers), 
     &          COORD(3,max_shells),Cord(3,max_centers),
     &          Ivshloff(max_shells),Nfps(max_shells),
     &          End_nfps(max_shells),Npfs(max_shells),
     &          End_nfps_uc(max_shells),
     &          ivAngMom(max_shells),Ixalpha(max_shells),
     &          Ixpcoef(max_shells),Indx_cc(max_shells),
     &          Ixshells(max_shells),Atm_4shell(max_shells),
     &          Temp(max_shells)
C
      Iecp = 1
      Ispherical = 1
      INP_Fname = "MOL_ERDOED"
      Dkh_order =   IFLAGS2(167)

      Call SIMPLE_INSPECT_MOL(INP_Fname, max_centers, max_shells,
     &                        ncenters, nshells, nspc, cartesian, 
     &                        ITFCT, LNP1, lnpo, nfct, nufct, 
     &                        nbasis, NAOBASIS, nCFpS, 
     &                        nPFpS, NAOATM, angmom, atomic_label, 
     &                        vnn, Maxang, Maxjco, Iecp,
     &                        NUcstr, Nrcstr, Iqmstr, Jcostr, 
     &                        Jstrt, Jrs, Atmnam, Charge)



      nalpha    = 0
      npcoef    = 0
      npcoef_uc = 0 
      do i = 1, nshells
         nalpha    = nalpha + npfps(i)
         npcoef    = npcoef + npfps(i) * ncfps(i)
         npcoef_uc = npcoef_uc + npfps(i) * npfps(i) 
         Ivangmom(i) = Angmom(i)
      enddo

      IALPHA    = 1
      IPCOEF    = IALPHA + nalpha
      IPCOEF_UC = IPCOEF + npcoef
      INEXT     = IPCOEF_UC + npcoef_uc

      do i = 1, ncenters
         npop(i)     = 1
         ireorder(i) = i
         Namat(i) = ATMNAM(i)
      enddo
      Nshll = 0
      do iat=1,Ncenters
        do ilq=1,iqmstr(iat)
          do ijco=1,jcostr(iat,ilq)
             nshll=nshll+1
             Atm_4shell(nshll)=iat
C      write(*,*)'Atm(nshll)',Atm_4shell(nshll)
           Enddo
        Enddo
      Enddo
C
      Call READ_BASIS_INFO(INP_FNAME, Ncenters, Ncenters, NPOP,
     &                     IREORDER,
     &                     CARTESIAN, ITFCT, LNP1, LNPO,
     &                     NFCT, NBASIS, WORK(IALPHA), IXALPHA,
     &                     Work(IPCOEF), IXPCOEF, MAX_CENTERS,
     &                     ATMNAM, COORD, CORD, NAOATM)


      Ixpcoef(1) = 1
      do i = 1, Nshells 
          Ixpcoef(i+1) = Npfps(i) * Npfps(i) + Ixpcoef(i)
      enddo 
C
C Here the number of contracted functions are the same as primitive 
C functions. 

      Call int_gen_init(Nshells, Ivshloff, Ivangmom, Nfps, Npfs, 
     &                  Npfps,
     &                  Ixshells, Ixalpha, Ixpcoef, End_nfps, Temp, 
     &                  Ispherical)


      call dzero(work(ipcoef_uc), npcoef_uc)

      Index = 0
      Do ishell=1, Nshells 
         ishel_off = Ixpcoef(ishell) 
         iprim_off = npfps(ishell)
         index = 0
         Do i=1, iprim_off
            work(ipcoef_uc+ishel_off+index-1) = 1.0D0
            Index = Index + iprim_off + 1
         Enddo 
      Enddo 
C
C Here the number of contracted functions are the same as primitive 
C functions. 
C
      Call Gen_oed_to_vmol(Nfps, Ivangmom, Nshells, Spherical,
     &                     Erd_index,ERD_Scale)

      Call setup_ccbeg(WORK(Ialpha), IXalpha, Work(Ipcoef_uc),Ixpcoef,
     &                 Npfps, Npfps, Nshells, CCbeg, CCend, Indx_cc)
C
 
      Call Erd_scratch_mem_calc(nshells, ivangmom, npfps, npfps,
     &                          atom, Coord, Work(ialpha),
     &                          Work(ipcoef_uc), 
     &                          ixalpha, ixpcoef,
     &                          Ccbeg, Ccend, indx_cc,
     &                          spherical, Ncenters, .true., 
     &                          intmax, zmax)

      Call Erd__memory_hrr_correction(ivAngMom, nshells,
     &                                spherical, ihrr, zhrr)
     
      Imem   = iflags(36)
      Intmax = Intmax + Ihrr
      Zmax   = Zmax   + Zhrr

       
      NBFNS   = NAOBASIS
      NAOBFNS = NAOBASIS

      N_Cart_fns = Naobasis
      N_Sphe_fns = Nbasis

      Nprims = 0
      ncnfns = 0
      Do Ishell = 1, Nshells 
         If (Spherical) Then
            lvalue = Ivangmom(Ishell) 
            Nprims = Nprims + (2*lvalue+1) *
     &               npfps(ishell) 
            Ncnfns = ncnfns  + (2*lvalue+1) * ncfps(ishell)
         Else 
            lvalue = Ivangmom(Ishell) 
            Nprims = Nprims + (lvalue+1)*(lvalue+2)/2 * 
     &               npfps(ishell) 
            Ncnfns = Ncnfns + (lvalue+1)*(lvalue+2)/2 * 
     &               ncfps(ishell) 
         Endif 
     
      Enddo 
      write(6,*)
      Write(6,"(a,1x,I4,1x,I4)") " Nprims,Ncnfns: ", Nprims,Ncnfns
      NBFNS = Nprims 

      Ioed_dbuf     = INext
      Ioed_ovl      = Ioed_dbuf   + Zmax 
      Ioed_kin      = Ioed_ovl    + Nprims**2
      Ioed_nai      = Ioed_kin    + Nprims**2
      Ioed_pvp      = Ioed_nai    + Nprims**2
      Icoreham      = Ioed_pvp    + Nprims**2
      Idkh          = Icoreham    + Nprims**2 
      Ireord        = Idkh        + Nprims**2
      IEnd          = Ireord      + 40*Nprims**2 
      Ileft         = Mxdcor - Iend 

      If (Iend .Ge. Mxdcor) Call insmem("Dkh:", Iend, Mxdcor)
      If (Intmax .Ge. Mxicor) Call insmem("Dkh:",Intmax,Mxicor)
C
C  There is one segment in serial work (not necessary, but lets work
C  with that in mind).
C
      Val1(1) = 1
      Val1(2) = 1
      Val1(3) = 1
      Val1(4) = 1
      Val2(1) = Nprims
      Val2(2) = Nprims 
      Val2(3) = Nprims 
      Val2(4) = Nprims 

   
      Call compute_oed_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ioed_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ipcoef_uc),
     &                       Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, 
     &                       Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,Contract,dkh_order,
     &                       Work(Ioed_kin),
     &                       Work(Ioed_ovl),
     &                       Work(Ioed_nai),
     &                       Work(Ioed_pvp), 
     &                       Work(Icoreham),
     &                       Work(Ireord),
     &                       Erd_index,Erd_scale,
     &                       Nsend,
     &                       Nalpha,Nprims,Ncnfns,Npcoef,Npcoef_uc,
     &                       Ncenters,
     &                       charge,
     &                       max_centers,Work(Idkh),Ileft)
C
C Contracted DKH integrals are returned in Work(Ioed_pvp). First transform
C them to ACES II ordering (and normalizations). Reuse the memory locations.
C Also reset the Nbfns to the number of basis functions (it was set to
C number of primitives.

      Nbfns = N_Sphe_fns

      Write(6,"(a)") "Unscaled contracted DKH integrals"
      Call output(Work(Ioed_pvp), 1, Nbfns, 1, Nbfns, Nbfns,
     &            Nbfns, 1)
       Call rescale_2ints_op2(Work(Ioed_pvp),Work(Ioed_kin),
     +                        Nbfns,Erd_index,Erd_scale,Symmetry)

      Write(6,"(a)") "Scaled/reordered contracted DKH integrals"
      Call output(Work(Ioed_kin), 1, Nbfns, 1, Nbfns, Nbfns,
     &            Nbfns, 1)
C
       Ndim = (Nbfns * (Nbfns+1))/2
       Call Squez2(Work(Ioed_kin), Work(Icoreham), Nbfns)
       Call Putrec(20,"JOBARC","DKH_INTS",Ndim*IINTFP,Work(Icoreham))

      Return
      End
