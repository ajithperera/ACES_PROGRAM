










      subroutine heffdiag(iuhf,icycle)
      implicit integer (i-n)
      implicit real*8  (a-h,o-z)


cjp
cjp data for multireference state specific Brillouin-Wigner CC method
cjp coded by Jiri Pittner 1998-2000
cjp
      logical isbwcc,masik,isactive,bwgossip,useeq429,scfrefread
      logical bwwarning
      character*256 bwwarntext
      real*8 ecorrbw0,ecorrbw,epsilon0,fockcontr,denomblow,fockcd
      real*8 heff,heffevalr,heffevali,heffevecl,heffevecr,hdiagcontr
      real*8 fock2elcontr,enerscf,hcore,lambdahomotop,hfakt,diishonset
      real*8 fock2elcontr0,enerscf0,cbwstate,totmaxdenom,heffevecrold
      real*8 intruder,hfaktmax
      integer maxorb,maxref,nref, iref, iocc,iocc0
      integer iphnum,invpnum,invhnum,nbwstates,ibwstate
      integer ibwconvg,internfrom,internto,internnum,internindex
      integer ibwpass
      integer nactive,numactive,ihubaccorr,ihomotop
      integer internfrom1,internto1,internindex1,internnum1
      integer ihefferank, iheffefrom, iheffeto, iheffespin, maxexcit
      integer correctiontype
      integer maxbwwarnings, nproc, myproc
c
      parameter(maxorb=512,maxref=32,maxexcit=9,maxbwwarnings=10)
c      NOTE!!! change of maxorb parameter requires format change
c              and character* change in bwread routine!!!
      parameter(denomblow=1d250)
c

cjp common has been splitted in order to avoid problems
cjp with padding on different 32 and 64 bit architectures

      common/bwccint/isbwcc, masik, nref, iref,iocc(maxorb,maxref,2),
     +     iphnum(maxorb,maxref,2),invpnum(maxorb,maxref,2),
     +     invhnum(maxorb,maxref,2),
     +     isactive(maxorb,2),nbwstates,ibwstate(maxref+1),
     +     internfrom(maxref*(maxref-1)/2,maxref,3),
     +     internto(maxref*(maxref-1)/2,maxref,3),
     +     internindex(maxref*(maxref-1)/2,maxref,3),
     +     internnum(maxref,3),
     +     internfrom1(maxref*(maxref-1)/2,maxref,2),
     +     internto1(maxref*(maxref-1)/2,maxref,2),
     +     internindex1(maxref*(maxref-1)/2,maxref,2),
     +     internnum1(maxref,2),
     +     ibwpass,ibwconvg(maxref),bwgossip,useeq429,
     +     nactive(2),numactive(maxorb,2),ihubaccorr, ihomotop,
     +     iocc0(maxorb,2),scfrefread,
     +     ihefferank(maxref,maxref),iheffefrom(maxexcit,maxref,maxref),
     +     iheffeto(maxexcit,maxref,maxref),
     +     iheffespin(maxexcit,maxref,maxref),
     +     correctiontype,bwwarning(maxbwwarnings),
     +     bwwarntext(maxbwwarnings),nproc,myproc

      common/bwccreal/ecorrbw,epsilon0,cbwstate(maxref+1),
     +     fockcontr(maxorb*(maxorb+1)/2,2),fockcd(maxorb,maxref,2),
     +     heff(maxref,maxref),heffevalr(maxref),heffevali(maxref),
     +     heffevecl(maxref,maxref),heffevecr(maxref,maxref),
     +     hdiagcontr(maxref), fock2elcontr(maxorb,2),
     +     enerscf(maxref), hcore(maxorb,2),
     +     lambdahomotop,hfakt,diishonset,enerscf0,
     +     fock2elcontr0(maxorb,2),ecorrbw0,totmaxdenom,
     +     heffevecrold(maxref,maxref),intruder,hfaktmax

c
c
cjp BRIEF DESCRIPTION OF VARIABLES INTRODUCED FOR THE MR-BWCC ROUTINES
cjp IN FACT, A LOT OF THAT COULD BE USEFUL FOR ANY HILBERT-SPACE MR-CC
c
c
c isbwcc ... flag for doing bwcc calculation
c maxbwwarnings, bwwarning, bwwarntext ... serious warnings will be
c    summarized at the end of xvcc output for the user's' convenience
c ihefferank(jref,iref) ... degree of excitation between jref and iref
c iheffefrom(maxexcit,jref,iref) , iheffeto, iheffespin ... list of
c    indices of that excitation, sorted according to spin and then the
c    indices, numbers stored are defined as effective particle-hole
c    indices of reference iref
c ihomotop ... whether to use homotopic transition to the
c    size-extensivity correction, after which iteration (if .ne.0)
c lambdahomotop scaling factor of the geometrical series of
c    lambda 1->0 transition
c hfakt ... current value of the homotopy parameter
c hfaktmax ... maximal homotopy parameter allowed to consider cc
c    equations converged
c diishonset ... at which value of hfact restart diis convergence acceleration
c masik ... prepare sorted integral file for the program by Masik and stop
c nref ... number of reference configurations
c iref ... current reference configuration and fermi vacuum
c bwgossip ... switch on debugging output
c ibwpass ... routines like newt2 have to be splitted in two passes -
c    construction of Heff and amplitude update after heff is diagonalized
c    for backw. compatibility, instead of introducing a new routine
c    the same routine does different things being called twice with
c    different ibwpass value
c ecorrbw ... correlation energy from BWCC - Heff(iref,iref) ...
c    denominator correction
c ecorrbw0 .... ecorrbw, but not scaled by the homotopic factor hfact
c denomblow ... huge number to cause division underflow - used for
c    zeroing out the internal amplitudes automatically
c nactive(spin): total count of active spinorbitals
c numactive(i=1..nactive,spin): number of i-th active spinorbital
c    in sequential numbering
c isactive(maxorb,spin): belongs given orbital to the active space?
c    for RHF, the beta ones must be initialized to be identical with alpha ones
c iocc(maxorb,1..nref,spin): defines the nref reference configurations
c    for both RHF and UHF:  iocc(i,iref,spin)=0 or 1
c iphnum(orbital no, iref, spin): gives the effective number of orbital
c    (both particle and hole ones are counted starting from 1)
c invpnum(eff.p.orb.,iref,spin): gives true orbital no. from the
c    effective particle one
c invhnum(eff.h.orb.,iref,spin): gives true orbital no. from the
c    effective hole one
c    all these three ones must be in RHF case initialised to be equal in the
c    alpha and beta parts to keep the code unique and simple
c internfrom(sequence counter n,iref,ispin) is for the ab spin case, the
c    other ones have to be iuhf-indexed
c    internfrom, internto: they are first and second index of n-th internal
c    excitation when processing reference given by second index to the array
c
c internindex(sequence counter,iref,ispin) is the position of
c    corresponding denominator in the denominator list
c internnum(iref,ispin) - number of internal excitation in that
c    category = max sequence counter here ispin=1,2,3 for AA,BB,AB
c internfrom1 etc. are analogous quantities for monoexcitations, here
c    ispin =1,2 note for later: all intern.... quantities are irrep-specific!
c fockcontr(findex(i,j),ispin) ... addition to the fock matrix of
c    reference no.1 to obtain the fock matrix of current reference
c    (fermi vacuum)
c fockcd(i,iref,ispin) ... diagonal part of that correction for ref. no. iref.
c hcore(i,ispin) ... one electron diagonal hamiltonian elements
c fock2elcontr(i,ispin) ... 2el contribution to the diagonal fock element
c    used temporarily
c hdiagcontr(iref) ... contribution of differences of HF energies of
c    different Fermi vacua to diagonal Heff elements
c enerscf(iref) ... HF energy of iref-th fermi vacuum
c ihubaccorr ... =1 ... calculate the size extenzivity correction for BWCC
c                =2,3 ... second and third pass of that calculation
c iocc0(maxorb,2) ... like iocc, but for dummy reference configuration
c    corresponding to SCF WF
c fock2elcontr0(maxorb,2) ... like fock2elcontr0 but for dummy reference
c    configuration
c enerscf0 ... like enerscf, but for dummy reference config
c scfrefread ... tells to bwprep that SCF reference has been read from input
c    and should not be generated automatically from nocc(ispin)
c nbwstates ... how many states to average
c ibwstates(1..nbwstates),cbwstates() ... their numbers and coefficients
c correctiontype ... 0=DC,L T2 term is removed/scaled, 1=DC/L term is not
c    removed/scaled
c totmaxdenom ... max 1/denom found for given reference's' fermi vacuum -
c    as indication of possible intruder problem
c intruder ... limit of 1/denom to be considered intruder and its
c    amplitude zeroed
c for parallelization
c nproc ... number of processors (counted from 1)
c myproc ... number of the processor currently executing the code
c    (counted from 1)


cjp
cjp this routine diagonalizes the effective hamiltonian from common bwcc
cjp selects the root wanted and stores the energy in common bwcc
cjp
c
      logical o
      dimension heff2(maxref,maxref)
      dimension work(8*maxref)
      dimension ip(maxref)
      dimension wr(maxref),wi(maxref),vl(maxref,maxref),
     &          vr(maxref,maxref)
cjp tolerance how large the imaginary part is neglected 
      parameter(tolimag=1e-5)

cjp parallelization
cjp broadcast columns of Heff between processors
cjp in a parallel run
       o=.true.
c
cjp continue diagonalization on all processors
        if(o)write(6,*)
        if(o)write(6,*)
cjp test if heff is a normal matrix
      tmp=0d0
      do i=1,nref
      do j=1,nref
      heff2(i,j)=0d0
      do k=1,nref
        heff2(i,j)=heff2(i,j)+heff(i,k)*heff(j,k)-heff(k,i)*heff(k,j)
      enddo
      if(dabs(heff2(i,j)).gt.tmp) tmp=dabs(heff2(i,j))
      enddo
      enddo
        if(o)write(6,*) 'departure of Heff from normal matrix',tmp
        if(o)write(6,*)
      do i=1,nref
       if(o)write(6,123)  (heff2(i,j),j=1,nref)
      enddo
        if(o)write(6,*)
        if(o)write(6,*)


cjp dgeev is destructive on input matrix
      if(nref.eq.1) then
cjp trivial
      heffevalr(1)=heff(1,1)
      heffevali(1)=0d0
      heffevecl(1,1)=1d0
      heffevecr(1,1)=1d0
      else
      do i=1,nref
      do j=1,nref
       heff2(j,i)=heff(j,i)
      enddo
      enddo
cjp
cjp NOTICE!!! dgeev normalizes eigenvectors in a strange way
cjp not to have L**H . R = 1, but euclidean norm of individual eigenvectors to be one
cjp this is to be considered if they should be used in any later computation!!!
cjp
      call xgeev('V','V', nref, heff2, maxref, wr, wi, vl, maxref, vr,
     &            maxref,work, 8*maxref, info )
      if (info.ne.0) stop 'ERROR in Heff diagonalization'
cjp permute the eigenstuff in increasing order of real part of eigenvalue
      do i=1,nref
      ip(i)=i
      heffevalr(i)=wr(i)
      enddo
      call sort2(nref,heffevalr,ip)
      

cjp store permuted results
      do i=1,nref
      ii=ip(i)
      heffevalr(i)=wr(ii)
      heffevali(i)=wi(ii)
      do j=1,nref
      heffevecl(j,i)=vl(j,ii)
      heffevecr(j,i)=vr(j,ii)
      enddo
      enddo

cjp
cjp renormalize eigenvectors to become bi-unitary
cjp this calculation is valid only if they are real-valued!
cjp do only if it turns out they will be used later
cjp

cjp end of the nontrivial case
      endif

cjp print results:
      if(o)write(6,*) 
      if(o)write(6,*) 'MR-BW-CC EFFECTIVE HAMILTONIAN'
        if(o)write(6,*)

      do i=1,nref
c      if(o)write(6,*) (heff(i,j),j=1,nref)
       if(o)write(6,123)  (heff(i,j),j=1,nref)
      enddo
        if(o)write(6,*)
      if(o)write(6,*) ' results of diagonalization (for complex 
     & eivectors for storage cf. man dgeev)'
123   format(10f16.10)
       if(o)write(6,*)
      if(o)write(6,*) 'eigenvalues (real and imaginary)'
        if(o)write(6,*)
        if(o)write(6,123) ( heffevalr(i),i=1,nref)
        if(o)write(6,123) ( heffevali(i),i=1,nref)
       if(o)write(6,*)
      if(o)write(6,*) 'right eigenvectors'
        if(o)write(6,*)
        do j=1,nref
        if(o)write(6,123) ( heffevecr(j,i),i=1,nref)
        enddo
       if(o)write(6,*)
      if(o)write(6,*) 'left eigenvectors'
        if(o)write(6,*)
        do j=1,nref
        if(o)write(6,123) ( heffevecl(j,i),i=1,nref)
        enddo
        if(o)write(6,*)
      if(iuhf.gt.0) then
        if(o)write(6,*)
      if(o)write(6,*) 'NOTE for spin coupling of open shell references'
      if(o)write(6,*) 'references like ...ab, ...ba actually corr. to'
      if(o)write(6,*) '                ...hl, ...hl'
      if(o)write(6,*)
     +    'determinants |...,h-alpha,l-beta| and |...,l-alpha,h-beta|'
      if(o)write(6,*) 'not |...,h-beta,l-alpha|.Therefore triplet state'
      if(o)write(6,*) 'corr. to coupling of ...ab, ...ba with - sign'
        if(o)write(6,*)
      endif
        if(o)write(6,*)

cjp compute some simple criteria for tworeference problems
      if(nref.eq.2) then
      if(o)write(6,*)'|H01/H00| = ',abs(heff(1,2)/heff(1,1))
      if(o)write(6,*)'2*SQRT(H01*H10)/|H00-H11| = ',
     +   2d0*sqrt(heff(1,2)*heff(2,1))/abs(heff(1,1)-heff(2,2))
      if(o)write(6,*)'E0-H00 = ',heffevalr(1)-heff(1,1)
      endif

cjp test for swapped eigenvectors (only for those we are interested in)
      if(icycle.gt.1) then
      do i=1,nbwstates
        k=ibwstate(i)
        s=0.
        do j=1,nref
        s=s+heffevecrold(j,k)*heffevecl(j,k)
        enddo
        if(abs(s).lt.0.5) then
         if(o)write(6,*)'@HEFFDIAG-W: Eivec ',k,' swapped!!! icycle= ',
     &             icycle
         bwwarning(3)=.true.
         bwwarntext(3)='@HEFFDIAG-W: Eigenvector swap occured'
        endif
      enddo
      endif
cjp store eigenvectors for this test in next iteration
      do i=1,nref
      do j=1,nref
      heffevecrold(j,i)=heffevecr(j,i)
      enddo
      enddo


cjp test for complex roots
      iflag=0
      do i=1,nref
        if(abs(heffevali(i)).ge.tolimag) iflag=iflag+1
      enddo
      if(iflag.ne.0) then
          if(o)write(6,*)'@HEFFDIAG-W: ',iflag,
     &       ' complex Heff eigenvalues detected'
          bwwarning(4)=.true.
          bwwarntext(4)='@HEFFDIAG-W: complex Heff eigenvalue detected'
      endif

cjp   write the energy for later use
      epsilon0=0.0
      do i=1,nbwstates
cjp problem if the requested eigenvalue is complex
          if(abs(heffevali(ibwstate(i))).ge.tolimag) then
           if(o)write(6,*) '@HEFFDIAG-E: MR-BW-CC energy is complex'
           stop '@HEFFDIAG-E: MR-BW-CC energy is complex'
          endif
          epsilon0=epsilon0+heffevalr(ibwstate(i))*cbwstate(i)
      enddo

      if(bwgossip.and.o)write(6,*)'CORRELATION ENERGY FROM MR-BW-CC = ',
     &     epsilon0
       if(o)write(6,*)
      return
      end
