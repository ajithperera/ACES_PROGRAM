










      function qcomp(n,l)
c     ----- computes q(n,l)                                -----
c     ----- scaled by exp(-t) to prevent overflows         -----
c     ----- arguments are alpha, beta, and xval=beta**2/(4*alpha) -----
c     ----- no restriction on the magnitude of xval           -----

      implicit double precision (a-h,o-z)

CSSS      common /pifac/ pi, sqpi2
CSSS      common /dfac/ dfac(29)
CSSS      common/qstore/dum1(81),alpha,xk,t

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

c These parameters are gathered from vmol and vdint and are used by ecp
c as well. It just so happens that the vmol parameters do not exist in
c vdint and vice versa. LET'S TRY TO KEEP IT THAT WAY!

c VMOL PARAMETERS ------------------------------------------------------

C     MAXPRIM - Maximum number of primitives for a given shell.
      INTEGER    MAXPRIM
      PARAMETER (MAXPRIM=72)

C     MAXFNC  - Maximum number of contracted functions for a given shell.
C               (vmol/readin requires this to be the same as MAXPRIM)
      INTEGER    MAXFNC
      PARAMETER (MAXFNC=MAXPRIM)

C     NHT     - Maximum angular momentum
      INTEGER    NHT
      PARAMETER (NHT=7)

C     MAXATM  - Maximum number of atoms
      INTEGER    MAXATM
      PARAMETER (MAXATM=100)

C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers.
      INTEGER    MXTNPR
      PARAMETER (MXTNPR=MAXPRIM*MAXPRIM)

C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers.
      INTEGER    MXTNCC
      PARAMETER (MXTNCC=180000)

C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers.
      INTEGER    MXTNSH
      PARAMETER (MXTNSH=200)

C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions).
c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

c VDINT PARAMETERS -----------------------------------------------------

C     MXPRIM - Maximum number of primitives for all symmetry
C              inequivalent centers.
      INTEGER    MXPRIM
      PARAMETER (MXPRIM=MXTNPR)

C     MXSHEL - Maximum number of shells for all symmetry inequivalent centers.
      INTEGER    MXSHEL
      PARAMETER (MXSHEL=MXTNSH)

C     MXCORB - Maximum number of contracted basis functions.
      INTEGER    MXCORB
      PARAMETER (MXCORB=MXCBF)

C     MXORBT - Length of the upper or lower triangle length of MXCORB.
      INTEGER    MXORBT
      PARAMETER (MXORBT=MXCORB*(MXCORB+1)/2)

C     MXAOVC - Maximum number of subshells per center.
      INTEGER    MXAOVC,    MXAOSQ
      PARAMETER (MXAOVC=32, MXAOSQ=MXAOVC*MXAOVC)

c     MXCONT - ???
      INTEGER    MXCONT
      PARAMETER (MXCONT=MXAOVC)

C
C Basic parameters: Maxang set to 7 (i functions) and Maxproj set
C 5 (up to h functions in projection space).

      Parameter(Maxang=7, Maxproj=6, Lmxecp=7, Mxecpprim=Mxprim*Mxatms)
     &        

      Parameter(Maxangpwr=(Maxang+1)**2,Lmnpwr=(((Maxang*(Maxang+2)*
     &         (Maxang+4))/3)*(Maxang+3)+(Maxang+2)**2*(Maxang+4))/16)

      Parameter(Lmnmax=(Maxang+1)*(Maxang+2)*(Maxang+3)/6,
     &          Lmnmaxg=(Maxang+1)*(9+5*Maxang+Maxang*Maxang)/3)

      Parameter(Ndico=10,Ndilmx=Maxang,
     &          Ndico2=ndico*Ndico,Maxang2=((Maxang+1)**2)*
     &          ((Maxang+2)**2)/4)
C
      Parameter(Maxints_4shell=Ndico2*Maxang2)
C
C In principle Maxmem only need to be (2*Maxang+1)**2. So, the 
C current setting is very generous. 

      Parameter(Maxmem = 50000)
   
      Parameter(Rint_cutoff = 25.32838, Eps1 = 1.0D-15, Tol=46.0561)
C46.0561)

C
C This file contain all the ECP variables that need to be known
C across multiple files.
C
C
      common /ECP_INT_VARS/Zlm(Lmnpwr), Lmnval(3,Lmnmax),
     &                     Istart(0:Maxang),Iend(0:Maxang),
     &                     Ideg(0:Maxang),Lmf(Maxangpwr),
     &                     Lml(Maxangpwr),
     &                     Lmx(Lmnpwr),Lmy(Lmnpwr),Lmz(Lmnpwr),
     &                     Pi,Fpi,Sqpi2,Sqrt_Fpi,R_intcutoff
     
      Common/ECP_INTGRD_VARS/Ideg_grd(0:Maxang), 
     &                       Istart_grd(0:Maxang),Iend_grd(0:Maxang),
     &                       Lmnval_grd(7,Lmnmaxg)

      common/ECP_POT_VARS/clp(Mxecpprim),zlp(Mxecpprim),
     &                    nlp(Mxecpprim),kfirst(Maxang,Mxatms),
     &                    klast(Maxang,Mxatms),llmax(Mxatms)

      common /pseud / nelecp(Mxatms),ipseux(Mxatms),ipseud 

      common /nshel / expnt(Mxtnpr),contr(Mxtnpr,Mxtnpr),
     &                numcon(Mxtnpr),katom(Mxtnsh),ktype(Mxtnsh),
     &                kprim(Mxtnsh),kbfn(Mxtnsh),kmini(Mxtnsh),
     &                kmaxi(Mxtnsh),nprims(Mxtnsh),ndegen(Mxtnsh),
     &                nshell,nbf

      Common /Qstore/Alpha,Beta,Xval
     
      Common /RadAng_sums/Rad_Sum(Maxang,Maxang), 
     &                    Ang_sum(Maxang,Maxang)
   
      Common /Fints/Fijk(0:4*Maxang,0:4*Maxang,0:4*Maxang)

      common /factorials/Fact(0:2*Maxang),Fac2(-1:4*Maxang),
     &                   Faco(0:2*Maxang),
     &                   Bcoefs(0:2*Maxang,0:2*Maxang),
     &                   Fprod(2*Maxang, 2*maxang)
  


      dimension tmin(9)

      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,
     1   18.0d0,15.0d0/

      data am1,   a0,   accpow, accasy, a1,   a2,   a4
     1    /-1.0d0,0.0d0,1.0d-14,1.0d-10,1.0d0,2.0d0,4.0d0/

      if(mod(n+l,2).ne.0.or.n.le.l) go to 30
c     ----- use alternating series (n+l.le.22.and.l.le.10) -----

      if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(beta/(alpha+alpha))**l
      prefac=sqpi2*xkp*fac2(n+l+1-2)/
     1   ((alpha+alpha)**((n-l)/2)*dsqrt(alpha+alpha)*fac2(l+l+3-2))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
 10   if(num.eq.0) go to 20
      fnum=num
      term=term*fnum*xval/(xden*xc)
      xc=xc+am1
      sum=sum+term
      num=num+2
      xden=xden+a2
      go to 10
   20 qcomp=prefac*sum
      return

   30 if(xval.lt.tmin(min0(n,8)+1)) go to 60
c     ----- use asymptotic series (arbitrary n,l) -----
      xkp=(beta/(alpha+alpha))**(n-2)
      prefac=xkp*sqpi2/((alpha+alpha)*dsqrt(alpha+alpha))
      sum=a1
      term=a1
      faca=(l-n+2)
      facb=(1-l-n)
      xc=a1
   40 term=term*faca*facb/(a4*xc*xval)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(dabs(term/sum).lt.accasy) go to 50
      faca=faca+a2
      facb=facb+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return

c     ----- use power series (n+l.le.22.and.l.le.10) -----
   60 if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(beta/(alpha+alpha))**l
      prefac=dexp(-xval)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) prefac=prefac*sqpi2/dsqrt(alpha+alpha)
      xnum=(l+n-1)
      xden=(l+l+1)
      term=fac2(l+n+1-2)/fac2(l+l+3-2)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*xval*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
      end

