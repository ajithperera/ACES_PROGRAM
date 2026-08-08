










      Subroutine Sps(Nlq, Lprj, La, Lb, Ltot, Lamalo, Lamahi,
     &               Lamblo, Lambhi, Lamau, Lambu, Alpha2,
     &               Beta1, Beta2, prd, Dum, Rad2)

      Implicit Double Precision (A-H, O-Z)


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
  


      Dimension Rad2(0:2*Maxang,0:2*Maxang,0:2*Maxang), 
     &          Fctr(2*maxang+3), Sum(2*Maxang+3), Term(2*Maxang+3)

CSSS      common /dfac/ dfac(29)
CSSS      common/qstore/dum1(81),alpha,xk,t
CSSS
      Data a0,accrcy,a1s4,a1s2,a1 /0.0d0,1.0d-13,0.25d0,0.5d0,1.0d0/
C
      Write(6,"(a,1x,4(1x,F20.13))") "Prd,,beta1,beta2,Alp.",Prd,beta1,
     &                                beta2, Alpha
      Write(6,*) "Dum: ", Dum
      Write(6,*)

      Xka1 = Beta1
      Xkb1 = Beta2
      Alp  = Alpha

      Write(6,"(a,7(1x,I3))") "Ltot,La,Lb,Lamalo,Lamahi,Lamblo,Lambhi:",
     &    Ltot,La,Lb,Lamalo,Lamahi,Lamblo,Lambhi
      L     = Lprj + 1
      Lit   = La + 1
      Ljt   = Lb + 1
      Ljtm1 = Ljt - 1

      Ltot1 = Ltot + 1
      Lmalo = Lamalo + 1
      Lmahi = Lamahi + 1
      Lmblo = Lamblo + 1
      Lmbhi = Lambhi + 1

      Npi   = Nlq

      if(xka1.gt.xkb1) go to 10
      xka=xka1
      xkb=xkb1
      go to 12
   10 xka=xkb1
      xkb=xka1
c     ----- set up parameters for qcomp using xkb -----
   12 alpha=a1
      sqalp=dsqrt(alp)

      xk=xkb/sqalp
      t=a1s4*xk*xk
C
C Set Qstore common bolck, so that Qcomp can work.
C
      Beta = xk
      Xval = t

      prd=prd*dexp(-(dum-t))

      tk=xka*xka/(alp+alp)

      do 90 lama=lmalo,lmahi
      ldifa1=iabs(l-lama)+1
      if(xka1.gt.xkb1) go to 14
      la1=lama-1
      go to 16
   14 lb1=lama-1

   16 do 90 lamb=lmblo,lmbhi
      ldifb=iabs(l-lamb)
      nlo=ldifa1+ldifb
      nhi=(ltot1-mod(lit-ldifa1,2))-mod((ljt-1)-ldifb,2)
      if(xka1.gt.xkb1) go to 18
      lb1=lamb-1
      go to 20
   18 la1=lamb-1

c     ----- run power series using xka, obtaining initial    -----
c     ----- q(n,l) values from qcomp, then recurring up wards -----
c     ----- j=0 term in sum -----

   20 continue
      Write(6,"(a,3(1x,F10.7))") "Qstore varialbles :", Alpha, Xk, t
      Write(6, "(a,4(1x,I2))") "Nstart,Lama,Lamb: ", nlo,npi+nlo-1+la1,
     &           la1,lb1

      qold2=qcomp(npi+nlo-1+la1,lb1)/fac2(la1+la1+3-2)
      fctr(nlo)=a1
      sum(nlo)=qold2
      if(nlo.eq.nhi.and.tk.eq.a0) go to 60

c     ----- j=1 term in sum -----
      nprime=npi+nlo+la1+1
      qold1=qcomp(nprime,lb1)/fac2(la1+la1+3-2)
      if(nlo.ne.nhi) fctr(nlo+2)=fctr(nlo)
      f1=(la1+la1+3)
      fctr(nlo)=tk/f1
      term(nlo)=fctr(nlo)*qold1
      sum(nlo)=sum(nlo)+term(nlo)
      if(nlo.ne.nhi) go to 22
      qold2=fctr(nlo)*qold2
      qold1=term(nlo)
      go to 24

   22 nlo2=nlo+2
      sum(nlo2)=qold1
      if(nlo2.eq.nhi.and.tk.eq.a0) go to 60
   24 j=1
c     ----- increment j for next term -----
   30 j=j+1

      nprime=nprime+2
      f1=(nprime+nprime-5)
      f2=((lb1-nprime+4)*(lb1+nprime-3))
      qnew=(t+a1s2*f1)*qold1+a1s4*f2*qold2

      nlojj=nlo+j+j
      if(nlo.eq.nhi) go to 40
      nhitmp=min0(nlojj,nhi)

      do 38 n=nlo2,nhitmp,2
      nrev=nhitmp+nlo2-n
      fctr(nrev)=fctr(nrev-2)
   38 continue

   40 f1=(j*(la1+la1+j+j+1))
      fctr(nlo)=tk/f1
      if(nlojj.gt.nhi) go to 44
      nhitmp=nlojj-2
      term(nlojj)=qnew
      sum(nlojj)=term(nlojj)

      do 42 n=nlo,nhitmp,2
      nrev=nhitmp+nlo-n
      term(nrev)=fctr(nrev)*term(nrev+2)
       sum(nrev)=sum(nrev)+term(nrev)
   42 continue

      if(nlojj.eq.nhi.and.tk.eq.a0) go to 60
      qold2=qold1
      qold1=qnew
      go to 30

   44 qold2=fctr(nhi)*qold1
      qold1=fctr(nhi)*qnew
      term(nhi)=qold1
      sum(nhi)=sum(nhi)+term(nhi)
      if(nlo.eq.nhi) go to 47
      nhitmp=nhi-2
      do 46 n=nlo,nhitmp,2
      nrev=nhitmp+nlo-n
      term(nrev)=fctr(nrev)*term(nrev+2)
      sum(nrev)=sum(nrev)+term(nrev)
  46  continue

   47 do 48 n=nlo,nhi,2
   48 if(term(n).gt.accrcy*sum(n)) go to 30
   60 if(la1.ne.0) go to 62
      prefac=prd/sqalp**(npi+nlo+la1)
      go to 64
   62 prefac=prd*xka**la1/sqalp**(npi+nlo+la1)
   64 do 66 n=nlo,nhi,2
      Rad2(n-1,lamb-1,lama-1)=Rad2(n-1,lamb-1,lama-1)+prefac*sum(n) 

   66 prefac=prefac/alp
   90 continue
      return
      end

