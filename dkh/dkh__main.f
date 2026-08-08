














































































































































































































      Program Dkh

      Implicit Double Precision (A-H, O-Z)

      Logical Spherical, Cartesian
      Logical Contract,Symmetry



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




      Call aces_init(icore, i0,icrsiz, iuhf, .true.)

      Contract  = (IFLAGS2(168) .EQ.0)
      Symmetry  = (IFLAGS(60)  .EQ. 2) 
      Spherical = (IFLAGS(62) .EQ. 1) 
      If (Spherical) Cartesian  = .False.

      If (Symmetry .and. Spherical .and. Contract) Then
         Write(6,"(20x,a)") "      -----Warning-------"
         Write(6,"(a,a,a)") " When contracted spherical harmonic",
     +                      " basis set with Abelian symmetry"
         Write(6,"(a,a)")   " is requested the DKH transformtion",
     +                      " is carried out using Cartesian basis."
         Write(6,"(a,a)")   " As a result sym=on or off with",
     +                      " the spherical harmonic basis sets",
     +                      " can give slightly different results."
      Endif

C If no-symmetry, spherical or Cartesian much easier to handle
 
      If (Spherical) Then 
          If (symmetry) then
             Call Dkh__sph_sym(Icore(I0),Icrsiz,Iuhf,Symmetry,
     +                         Contract)
          Else
             Call Dkh__sph_nosym(Icore(I0),Icrsiz,Iuhf,Spherical,
     +                           Cartesian,Symmetry,Contract)
          Endif 
      Else
          Call Dkh__cart(Icore(I0),Icrsiz,Iuhf,Spherical,Cartesian,
     +                   Symmetry,Contract)
      Endif 
    
      Call Aces_fin

      Stop
      End
       
      
       

