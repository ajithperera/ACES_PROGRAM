










      Subroutine Make_ztab
 
      Implicit Double Precision (a-h,o-z)


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
  

C
C form of liner cartesian products. The table is generated
C using the expression provided by Richard J. Mathar. 
C
C -Note that (PRMTESTING block is used only in the early stages of
C development and no longer supported). Ajith Perera.
C
c compute tables by recursion for real spherical harmonics.  they
c are indexed by l, m and sigma.  the sequence number of the
c harmonic with quantum numbers l, m and sigma is given by
c l**2+2*m+1-sigma lmf(index) and lml(index) hold the positions of 
c the first and last terms of the harmonic in the lmx, lmy, lmz, 
c and zlm arrays. The harmonics with angular momentum l are generated 
c from those with angular momenta l-1 and l-2.
c for m = 0,1,2,...,l-1, the recursion relation
c z*Z(l-1,m,s) = sqrt(((l-m)*(l+m))/((2*l-1)*(2*l+1)))*Z(l,m,s)+
c                sqrt(((l+m-1)*(l-m-1))/((2*l-3)*(2*l-1)))*Z(l-2,m,s)
c is used.
c for m = l, the recursion relation
c x*Z(l-1,l-1,s)+(-1)**(1-s)*y*Z(l-1,l-1,1-s) =
c                sqrt((2*l))/((2*l+1)))*Z(l,l,s)
c is used.
c  l=0
      lmf(1) = 1
      lml(1) = 1
      lmx(1) = 0
      lmy(1) = 0
      lmz(1) = 0
      zlm(1) = 1.0D0/Sqrt_fpi
c l=1
      lmf(2) = 2
      lml(2) = 2
      lmx(2) = 0
      lmy(2) = 0
      lmz(2) = 1
      zlm(2) = sqrt(3.0D0)/Sqrt_fpi
      lmf(3) = 3
      lml(3) = 3
      lmx(3) = 0
      lmy(3) = 1
      lmz(3) = 0
      zlm(3) = zlm(2)
      lmf(4) = 4
      lml(4) = 4
      lmx(4) = 1
      lmy(4) = 0
      lmz(4) = 0
      zlm(4) = zlm(2)
      nterm=4

      do 270 lang=2,Maxang
        do 240 mang=0,lang-1
          anum = ((2*lang-1)*(2*lang+1))
          aden = ((lang-mang)*(lang+mang))
          coef1 = sqrt(anum/aden)
          anum = ((lang+mang-1)*(lang-mang-1)*(2*lang+1))
          aden = (2*lang-3)*aden
          coef2 = sqrt(anum/aden)
          nsigma=min(1,mang)
          do 230 isigma=nsigma,0,-1
            indexh=lang**2+2*mang+1-isigma
            lone=lang-1
            ltwo=lang-2
            ione=lone**2+2*mang+1-isigma
            itwo=ltwo**2+2*mang+1-isigma
            lmf(indexh)=lml(indexh-1)+1
            lml(indexh)=lml(indexh-1)
            nxy=(mang-isigma+2)/2
            iu=lmf(ione)+nxy-1
            do i=lmf(ione),iu
              lml(indexh)=lml(indexh)+1
              j=lml(indexh)
              lmx(j)=lmx(i)
              lmy(j)=lmy(i)
              lmz(j)=lmz(i)+1
              zlm(j)=zlm(i)*coef1
              nterm=nterm+1
            enddo
            if(ltwo.ge.mang) then
              il=iu+1
              do i=il,lml(ione)
                lml(indexh)=lml(indexh)+1
                j=lml(indexh)
                k=lmf(itwo)+i-il
                lmx(j)=lmx(k)
                lmy(j)=lmy(k)
                lmz(j)=lmz(k)
                zlm(j)=zlm(i)*coef1-zlm(k)*coef2
                nterm=nterm+1
              enddo
              il=lml(itwo)-nxy+1
              if(mod(lang-mang,2).eq.0) then
                do i=il,lml(itwo)
                  lml(indexh)=lml(indexh)+1
                  j=lml(indexh)
                  lmx(j)=lmx(i)
                  lmy(j)=lmy(i)
                  lmz(j)=lmz(i)
                  zlm(j)=-zlm(i)*coef2
                  nterm=nterm+1
                enddo
              endif
            endif
  230     enddo
  240   enddo

        anum = (2*lang+1)
        aden = (2*lang)
        coef = sqrt(anum/aden)
        mang=lang
        isigma=1
        indexh=lang**2+2*mang+1-isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lml(indexh-1)
c isig:  index of the harmonic (l-1),(m-1),sigma
c isigm: index of the harmonic (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        k=lmf(isigm)
        do i=lmf(isig),lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)+zlm(k))*coef
          k=k+1
          nterm=nterm+1
        enddo
        if(mod(mang,2).eq.1) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=zlm(k)*coef
          nterm=nterm+1
        endif
        isigma=0
        indexh=lang**2+2*mang+1-isigma
c isig:  index of the harmonic (l-1),(m-1),sigma
c isigm: index of the harmonuc (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lmf(indexh)
        j=lml(indexh)
        i=lmf(isig)
        lmx(j)=lmx(i)+1
        lmy(j)=lmy(i)
        lmz(j)=lmz(i)
        zlm(j)=zlm(i)*coef
        nterm=nterm+1
        k=lmf(isigm)
        do  i=lmf(isig)+1,lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)-zlm(k))*coef
          k=k+1
          nterm=nterm+1
        enddo
        if(mod(mang,2).eq.0) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          k=lml(isigm)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=-zlm(k)*coef
          nterm=nterm+1
        endif
  270 enddo
C---
CSS      ixy = 0
CSS      iz = 0
CSS      do 300 lang=1,lproju
CSS        do 290 mang=0,lang-1
CSS          nsigma=min(1,mang)
CSS          ndelta=max(0,1-mang)
CSS          anum = ((lang-mang)*(lang+mang+1))
CSS          aden = (2*(2-ndelta))
CSS          coef=sqrt(anum/aden)
CSS          do 280 isigma=nsigma,0,-1
CSS            isign=2*isigma-1
CSS            ixy = ixy+1
CSS            flmtx(1,ixy) = (isign)*coef
CSS            flmtx(2,ixy) = coef
CSS            if(mang.ne.0) then
CSS              iz=iz+1
CSS              flmtx(3,iz) = -(mang*isigma)
CSS            endif
CSS  280     enddo
CSS  290   enddo
CSS        iz=iz+1
CSS        flmtx(3,iz) = -(lang)
CSS  300 enddo
c  column and row indices for angular momentum matrix elements.
CSS      iadd = 1
CSS      do i=1,2*lproju-1
CSS        mc(1,i) = i
CSS        mc(2,i) = i
CSS        mc(3,i) = i+1
CSS        mr(1,i) = i+iadd
CSS        mr(2,i) = i+2
CSS        mr(3,i) = i+2
CSS        iadd = 4-iadd
CSS      enddo
C--
      Return
      End
