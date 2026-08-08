













































































































































































































      Program A3ecp
      
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
C#include "ecp.par"
C#include "ecp.com"
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

      Parameter (Maxbfns = 1000, Mxdcor=80000000, Mxicor= 1000000)
      Character*4 INP_Fname
      Character*4 ATMNAM(Max_centers), Namat(Max_centers)
      Logical Spherical, Cartesian 
      Logical Pvp
    
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
     &          ivAngMom(max_shells),Ixalpha(max_shells),
     &          Ixpcoef(max_shells),Indx_cc(max_shells),
     &          Ixshells(max_shells),Atm_4shell(max_shells),
     &          Temp(max_shells),CordM(3,max_centers),
     &          CordP(3,max_centers) 

cPV 
       logical delta_int, twoCent_int,FourCent_soi, angmom_int,
     &         efield_int, Dshield_int, Nattrac, fourcent_dshield,
     &         fourcent_repl,Fourcent_spnspn_cpling,pvp_ints,
     &         secnd_der_nai,Two_cent_fld_grad,ovl_der_int
cendPV
C
      Call aces_init(icore, i0,icrsiz, iuhf, .true.) 
C
      If (IFLAGS(71) .EQ. 1) Iecp = 1
      Ispherical = 0
      Spherical  = .False.
      cartesian  = .False.
      Iecp       = 1
      Eps        = 0.000001D0
      If (IFLAGS(62) .EQ. 1) Then
          Ispherical = 1
          Spherical  = .True.
      Else
          Cartesian  = .True.
      Endif 
C
      INP_Fname = "MOL"
      Call SIMPLE_INSPECT_MOL(INP_Fname, max_centers, max_shells,
     &                        ncenters, nshells, nspc, cartesian, 
     &                        ITFCT, LNP1, lnpo, nfct, nufct, 
     &                        nbasis, NAOBASIS, nCFpS, 
     &                        nPFpS, NAOATM, angmom, atomic_label, 
     &                        vnn, Maxang, Maxjco, Iecp,
     &                        NUcstr, Nrcstr, Iqmstr, Jcostr, 
     &                        Jstrt, Jrs, Atmnam, Charge)

C---
      Write(6,"(a,a,7(1x,I4))") "natoms,nshells,ITFCT,LNP1,lnpo,",
     &                          "nbasis,NAOBASIS:",ncenters,nshells,
     &                                         ITFCT,LNP1,lnpo,
     &                                        nbasis,NAOBASIS
      Write(6,*) 
      Write(6,"(a)") "The number of con. function per shell"
      Write(6,"(6(1x,I4))") (nCFpS(i), i=1, nshells)
      Write(6,"(a)") "The number of prim. function per shell"
      Write(6,"(6(1x,I4))") (nPFpS(i), i=1, nshells)
      Write(6,"(a)") "The Iqmstr"
      Write(6,"(6(1x,I4))") (Iqmstr(i), i=1, ncenters)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
      Enddo
      Write(6,"(a)") "The Nucstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nucstr(i,j,k),j=1,4),k=1,ITFCT)
      Enddo
      Write(6,"(a)") "The Nrcstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nrcstr(i,j,k),j=1,4),k=1,ITFCT)
      Enddo
      Write(6,"(a)") "The jstrt"
      Write(6,"(6(1x,I4))") (Jstrt(i), i=1, Nshells)
      Write(6,"(a)") "The jrs"
      Write(6,"(6(1x,I4))") (Jrs(i), i=1, Nshells)

      nalpha = 0
      npcoef = 0
      do i = 1, nshells
         nalpha = nalpha + npfps(i)
         npcoef = npcoef + npfps(i) * ncfps(i)
         Ivangmom(i) = Angmom(i)
      enddo

      IALPHA = 1
      IPCOEF = IALPHA + nalpha
      INEXT  = IPCOEF + npcoef

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

      Write(6,*)
      Write(6,"(a)") "The exponents"
      Write(6,"(6(1x,F12.8))") (Work(Ialpha+i), i=0, Nalpha-1)
      Write(6,*)
      Write(6,"(a)") "The contraction coefs"
      Write(6,"(6(1x,F12.8))") (Work(Ipcoef+i), i=0, Npcoef-1)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
      Enddo

      Write(6,*)
      Write(6,*) "Atom Coordiantes of each shell"
      Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
      Enddo

      Call int_gen_init(Nshells, Ivshloff, Ivangmom, Nfps, Npfs, Ncfps,
     &                  Ixshells, Ixalpha, Ixpcoef, End_nfps, Temp, 
     &                  Ispherical)
      Write(6,*)
      Write(6,*) "The Ivangmom array" 
      Write(6,"(8(1x,I4))") (Ivangmom(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixshells array" 
      Write(6,"(8(1x,I4))") (Ixshells(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Nfps array" 
      Write(6,"(8(1x,I4))") (Nfps(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The End_Nfps array" 
      Write(6,"(8(1x,I4))") (End_Nfps(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The IVshloff array" 
      Write(6,"(8(1x,I4))") (Ivshloff(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixpcoef array" 
      Write(6,"(8(1x,I4))") (Ixpcoef(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixalpha array" 
      Write(6,"(8(1x,I4))") (Ixalpha(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The End_nfps array" 
      Write(6,"(8(1x,I4))") (End_nfps(i), i=1, Nshells)
      
      Call aces_to_erd(Nfps, Ivangmom, Nshells, ispherical, Erd_index,
     &                 ERD_Scale)

      Call setup_ccbeg(WORK(Ialpha), IXalpha, Work(Ipcoef), Ixpcoef,
     &                 Ncfps, Npfps, Nshells, CCbeg, CCend, Indx_cc)
C
      Write(6,*)
      Write(6,*)"The CCbeg and CCend arrays" 
      Write(6,"(10(1x,I4))") (CCbeg(i), i=1, Itfct)
      Write(6,"(10(1x,I4))") (CCend(i), i=1, Itfct)
      Write(6,*)"The Erd_index and Erd_scalae arrays" 
      Write(6,"(10(1x,I4))") (Erd_index(i), i=1,Nbasis)
      Write(6,"(10(1x,F5.2))") (Erd_Scale(i), i=1,Nbasis)
      Write(6,*)
 
      Call Erd_scratch_mem_calc(nshells, ivangmom, ncfps, npfps,
     &                          atom, Coord, Work(ialpha),
     &                          Work(ipcoef), ixalpha, ixpcoef,
     &                          Ccbeg, Ccend, indx_cc,
     &                          spherical, Ncenters, .true., 
     &                          intmax, zmax)

      Call Erd__memory_hrr_correction(ivAngMom, nshells,
     &                                spherical, ihrr, zhrr)
     
      Imem   = iflags(36)
      Intmax = Intmax + Ihrr
      Zmax   = Zmax   + Zhrr

      Write(6,"(a,1x,I6,1x,I10)") "ERD memmory requirments Ibuf & Dbuf", 
     &                            Intmax, Zmax
      NBFNS = NAOBASIS
      If (Spherical) NBFNS = Nbasis 

      Ierd_dbuf          = INext
      Ierd_d4cent_Jz     = Ierd_dbuf + Zmax
      Ierd_d4cent_Jx     = Ierd_d4cent_Jz + Nbfns**4
      Ierd_d4cent_Jy     = Ierd_d4cent_Jx + Nbfns**4
      Ioed_delta_int     = Ierd_d4cent_Jy + Nbfns**4
      Ioed_d2cent_intefx = Ioed_delta_int + (Nbfns**2)*Ncenters
      Ioed_d2cent_intefy = Ioed_d2cent_intefx + Nbfns**2
      Ioed_d2cent_intefz = Ioed_d2cent_intefy + Nbfns**2
      Ioed_d2cent_intlx  = Ioed_d2cent_intefz + Nbfns**2
      Ioed_d2cent_intly  = Ioed_d2cent_intlx  + Nbfns**2
      Ioed_d2cent_intlz  = Ioed_d2cent_intly  + Nbfns**2
      Iattrac            = Ioed_d2cent_intlz  + Nbfns**2
      Ioed_d2cent_intdxx = Iattrac            + Nbfns**2
      Ioed_d2cent_intdxy = Ioed_d2cent_intdxx + Nbfns**2
      Ioed_d2cent_intdxz = Ioed_d2cent_intdxy + Nbfns**2
      Ioed_d2cent_intdyx = Ioed_d2cent_intdxz + Nbfns**2
      Ioed_d2cent_intdyy = Ioed_d2cent_intdyx + Nbfns**2
      Ioed_d2cent_intdyz = Ioed_d2cent_intdyy + Nbfns**2
      Ioed_d2cent_intdzx = Ioed_d2cent_intdyz + Nbfns**2
      Ioed_d2cent_intdzy = Ioed_d2cent_intdzx + Nbfns**2
      Ioed_d2cent_intdzz = Ioed_d2cent_intdzy + Nbfns**2
      Ioed_d4cent_intdxx = Ioed_d2cent_intdzz + Nbfns**2 
      Ioed_d4cent_intdxy = Ioed_d4cent_intdxx + Nbfns**4
      Ioed_d4cent_intdxz = Ioed_d4cent_intdxy + Nbfns**4
      Ioed_d4cent_intdyx = Ioed_d4cent_intdxz + Nbfns**4
      Ioed_d4cent_intdyy = Ioed_d4cent_intdyx + Nbfns**4
      Ioed_d4cent_intdyz = Ioed_d4cent_intdyy + Nbfns**4
      Ioed_d4cent_intdzx = Ioed_d4cent_intdyz + Nbfns**4
      Ioed_d4cent_intdzy = Ioed_d4cent_intdzx + Nbfns**4
      Ioed_d4cent_intdzz = Ioed_d4cent_intdzy + Nbfns**4
      Ioed_4eri_int      = Ioed_d4cent_intdzz + Nbfns**4
      Ioed_4repls_int    = Ioed_4eri_int      + Nbfns**4
      Iaces_ord          = Ioed_4repls_int    + Nbfns**4
      Ioed_pvp_int       = Iaces_ord          + Nbfns**4
      Ioed_hess_int      = Ioed_pvp_int       + Nbfns**2
      Ioed_aces_ord      = Ioed_hess_int      + Nbfns**2 
      Ioed_nai_x1        = Ioed_aces_ord      + Nbfns**2
      Ioed_nai_y1        = Ioed_nai_x1        + Nbfns**2
      Ioed_nai_z1        = Ioed_nai_y1        + Nbfns**2
      Ioed_nai_x2        = Ioed_nai_z1        + Nbfns**2
      Ioed_nai_y2        = Ioed_nai_x2        + Nbfns**2
      Ioed_nai_z2        = Ioed_nai_y2        + Nbfns**2
      Ioed_nai           = Ioed_nai_z2        + Nbfns**2
      Ioed_ovl_x         = Ioed_nai           + Nbfns**2
      Ioed_ovl_y         = Ioed_ovl_x         + Nbfns**2
      Ioed_ovl_z         = Ioed_ovl_y         + Nbfns**2
      Iend               = Ioed_ovl_z         + Nbfns**2
      Write(6,"(a,1x,i5,1x,i5)") "Iend and Mxdcor:", Iend, Mxdcor

      If (Iend .Ge. Mxdcor) Call insmem("Main dcore:", Iend, Mxdcor)
      If (Intmax .Ge. Mxicor) Call insmem("Main icore:",Intmax,Mxicor)
C
C  There is one segment in serial work (not necessary, but lets work
C  with that in mind).
C
      Val1(1) = 1
      Val1(2) = 1
      Val1(3) = 1
      Val1(4) = 1
      Val2(1) = Nbfns
      Val2(2) = Nbfns
      Val2(3) = Nbfns
      Val2(4) = Nbfns
      Write(6,*)
      Write(6,*) "Atom Coordinates of each shell"
      Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
      Enddo

      Write(6,*)
      Write(6,*) "The atom label for each shell in main",Nshells
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)
        
      Nattrac=.False. 
      delta_int=.false. 
      twoCent_int=.false.
      FourCent_soi=.false.
      FourCent_Repl=.true.
      angmom_int=.false.
      efield_int=.false.
      dshield_int=.false.
      Fourcent_dshield=.false.
      Fourcent_spnspn_cpling=.false. 
      pvp_ints=.false.
      secnd_der_nai=.true.
      two_cent_fld_grad=.false.
      ovl_der_int=.true.

      if(ovl_der_int) then
      Call compute_1derovl_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_ovl_x),Work(Ioed_ovl_y),
     &                       Work(Ioed_ovl_z),
     &                       Nsend,
     &                       Nalpha, Npcoef,
     &                       max_centers)
      end if

      if(Nattrac) then

      Call compute_kin_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Iattrac),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)

      Call compute_nattrac_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Iattrac),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)

       end if

      if(FourCent_soi) then
      Call compute_4cent_integrals_SOI(val1(1),val2(1), val1(2),val2(2),
     &                       val1(3),val2(3),val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef, 
     &                       Ispherical,
     &                       Work(Ierd_d4cent_Jz), 
     *                       work(Ierd_d4cent_Jx),
     *                       work(Ierd_d4cent_Jy),Nsend, 
     &                       Nalpha, Npcoef,ERD_index,ERD_scale,
     &                       Work(Iaces_ord))
      end if

c      Write(*,*) 'Time to calculate 1electron SOI'
c      write(*,*) 'delta_int=',delta_int
c      write(*,*) 'twoCent_int=',twoCent_int

      if(delta_int) then
      Call compute_delta_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &   Work(Ioed_delta_int),
     *   Nsend, Nalpha, Npcoef, Ncenters,  max_centers)
      end if

       if(angmom_int) then
      Call compute_angmom_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intlx), 
     &                       Work(Ioed_d2cent_intly), 
     &                       Work(Ioed_d2cent_intlz),
     &                       Nsend, 
     &                       Nalpha, Npcoef, Ncenters,
     &                       max_centers)
       end if

      if(pvp_ints) then
      Call compute_pvp_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_pvp_int),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge, 
     &                       max_centers,ERD_index, 
     &                       ERD_scale,Work(Iaces_ord))
       end if

      if(secnd_der_nai) then

      Do Iatm = 1, 2
      Call Dcopy(3*Ncenters,Cord,1,CordM,1)
      Do Ixyz  = 3, 3
         CordM(Ixyz,Iatm) = Cord(Ixyz,Iatm) - EPs
      Do Jatm = 1, Iatm 
      Call Dcopy(3*Ncenters,Cord,1,CordP,1)
      Do Jxyz  = 3, 3
         CordP(Ixyz,Jatm) = Cord(Ixyz,Jatm) + EPs

         Write(6,"(a,4(1x,i3))") "Iatm,jatm,Ixyz,Jxyz: ",Iatm,Jatm,
     &                            Ixyz,Jxyz
      Write(6,"(a)") "-------------------------------------"        

      Call compute_1dernai_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       CordM,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_nai_x1),Work(Ioed_nai_y1),
     &                       Work(Ioed_nai_z1),Work(Ioed_nai),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge,
     &                       max_centers,ERD_index,
     &                       ERD_scale)

      Call compute_1dernai_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       CordP,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_nai_x2),Work(Ioed_nai_y2),
     &                       Work(Ioed_nai_z2),Work(Ioed_nai),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge,
     &                       max_centers,ERD_index,
     &                       ERD_scale)

      Call Perform_numder(Work(Ioed_nai_x1),Work(Ioed_nai_y1),
     &                    Work(Ioed_nai_z1),Work(Ioed_nai_x2),
     &                    Work(Ioed_nai_y2),Work(Ioed_nai_z2),
     &                    Work(Ioed_nai),Nbfns,Eps,Ixyz,Jxyz)
      Enddo
      Enddo
      Enddo
      Enddo
C#ifdef _NOSKIP
      Call compute_2dernai_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_pvp_int),Work(Ioed_hess_int),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,Charge,
     &                       max_centers,ERD_index,
     &                       ERD_scale)
C#endif 
      End if
C#ifdef _4CENTDSH
      if(efield_int) then
      Call compute_efield_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intefx), 
     &                       Work(Ioed_d2cent_intefy), 
     &                       Work(Ioed_d2cent_intefz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)

      Call compute_efield_integrals_test(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intefx),
     &                       Work(Ioed_d2cent_intefy),
     &                       Work(Ioed_d2cent_intefz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)


       end if
C#endif 

      if(dshield_int) then
      Call compute_dshield_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intdxx), 
     &                       Work(Ioed_d2cent_intdxy), 
     &                       Work(Ioed_d2cent_intdxz), 
     &                       Work(Ioed_d2cent_intdyx), 
     &                       Work(Ioed_d2cent_intdyy), 
     &                       Work(Ioed_d2cent_intdyz), 
     &                       Work(Ioed_d2cent_intdzx), 
     &                       Work(Ioed_d2cent_intdzy), 
     &                       Work(Ioed_d2cent_intdzz), 
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers)
       end if

      if(Fourcent_dshield) then
      Call compute_4cent_dshield_integrals(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       val1(3),val2(3),
     &                       val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d4cent_intdxx),
     &                       Work(Ioed_d4cent_intdxy),
     &                       Work(Ioed_d4cent_intdxz),
     &                       Work(Ioed_d4cent_intdyx),
     &                       Work(Ioed_d4cent_intdyy),
     &                       Work(Ioed_d4cent_intdyz),
     &                       Work(Ioed_d4cent_intdzx),
     &                       Work(Ioed_d4cent_intdzy),
     &                       Work(Ioed_d4cent_intdzz),
     &                       Work(Ioed_4eri_int), Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers,ERD_index,ERD_scale,
     &                       Work(Iaces_ord))
       end if

      if(Fourcent_spnspn_cpling) then
      Call compute_4cent_spnspn_cpling_integrals(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       val1(3),val2(3),
     &                       val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d4cent_intdxx),
     &                       Work(Ioed_d4cent_intdxy),
     &                       Work(Ioed_d4cent_intdxz),
     &                       Work(Ioed_d4cent_intdyx),
     &                       Work(Ioed_d4cent_intdyy),
     &                       Work(Ioed_d4cent_intdyz),
     &                       Work(Ioed_d4cent_intdzx),
     &                       Work(Ioed_d4cent_intdzy),
     &                       Work(Ioed_d4cent_intdzz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers,ERD_index,ERD_scale,
     &                       Work(Iaces_ord))
       end if

       if(fourCent_Repl) then
      Call compute_4cent_integrals(val1(1),val2(1),val1(2),val2(2),
     &                       val1(3),val2(3),
     &                       val1(4),val2(4),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &   Work(Ioed_4repls_int),
     &   Nalpha, Npcoef)
       end if

       if(twoCent_int) then
      Call compute_2cent_integrals_SOI(val1(1),val2(1),val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &   Work(Ioed_d2cent_intefx), Work(Ioed_d2cent_intefy),
     &   Work(Ioed_d2cent_intefz),Nsend,
     &   Nalpha, Npcoef, Ncenters,charge,max_centers)
       end if

      if(two_cent_fld_grad) then
      Call compute_2cent_fld_grad(val1(1),val2(1),
     &                       val1(2),val2(2),
     &                       Work(Ierd_dbuf),Zmax,
     &                       IWork(1),
     &                       Cord,Work(ipcoef),Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_d2cent_intdxx),
     &                       Work(Ioed_d2cent_intdxy),
     &                       Work(Ioed_d2cent_intdxz),
     &                       Work(Ioed_d2cent_intdyx),
     &                       Work(Ioed_d2cent_intdyy),
     &                       Work(Ioed_d2cent_intdyz),
     &                       Work(Ioed_d2cent_intdzx),
     &                       Work(Ioed_d2cent_intdzy),
     &                       Work(Ioed_d2cent_intdzz),
     &                       Nsend,
     &                       Nalpha, Npcoef, Ncenters,charge,
     &                       max_centers,ERD_index,ERD_scale,
     &                       Work(Ioed_aces_ord))
       end if
C
      Call aces_fin
      
      Stop
      End
