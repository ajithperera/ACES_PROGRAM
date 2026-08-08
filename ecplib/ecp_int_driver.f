










      subroutine Ecp_int_driver(Natoms, Ntotatoms, IGenby, Coord, Cint,
     &                          Grads)
c----------------------------------------------------------------------
c     Ecp_int_driver drives the calculation of 1e-matrices in the basis of
c     cartesian atomic orbitals using subroutine ecp_int.
c----------------------------------------------------------------------
c     tol    = integral cutoff (ints to be neglected, if value expected
c              to be < exp(-tol)(i)
c     Cint   = 1e-integral array ('batch') of dimension nfij
c       
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C-----------------------------------------------------------------------

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
  


CSSS      parameter (nftmax = (ndilmx+1)*(ndilmx+2)*(ndilmx*3)/6)
c     symmetry parameter:
C-----------------------------------------------------------------------
C     nftmax  total number of reducible gaussians up to {ndilmx+1}
C     SPECIAL parameter needed in gradient calculations
C-----------------------------------------------------------------------

      logical iandj,forget,Zero_int, Grads
C 
C Ecpints need to be (2*Maxang+1)**2 times the square of the 
C maximum number of contracted functions per shell (currently
C set at 20 for ECP integrals). The current Mxcbf is 1000. 

      dimension Cint(Maxmem), Ecpint_4Shell(Maxints_4shell),
     &          Ecpint(Mxcbf*(Mxcbf+1)/2), IGenby(Ntotatoms),
     &          Coord(3,Ntotatoms)

CSSS      common /infoa / xyz(3,Mxatms),charg(Mxatms),wmass(Mxatms),natoms
C
CSSS      common /symshe/ mulsh(Mxtnsh)
CSSS      common /modez / zetm(Mxtnsh)
CSSS      common /pairij/ ipq(ndi4+1)
CSSS      common /typinf/ ltmax,ltdrv,nftnft,nftdrv
CSSS      common /powers/ jx(nftmax),jy(nftmax),jz(nftmax),
CSSS     &                ix(nftmax),iy(nftmax),iz(nftmax)
CSSS      common/xa/xand(3,8,4*(Maxang-1)-3)
CSSS      common/symsh2/mij,kij,hkij,kb
C
      integer and,eor,or,dstrt
C
c-----------------------------------------------------------------------
c     common from AcesII (created in Readin)
c-----------------------------------------------------------------------
      COMMON /DAT/ EXPA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     &             CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA,TLC

      common/indx/pc(512),dstrt(8,MXCBF),ntap,lu2,nrss,nucz,itag,
     &            maxlop,maxlot,kmax,nmax,khkt(7),mult(8),isytyp(3),
     &            itype(7,28),and(8,8),or(8,8),eor(8,8),nparsu(8),
     &            npar(8),mulnuc(Mxatms),nhkt(MXTNSH),mul(MXTNSH),
     &            nuco(MXTNSH),nrco(MXTNSH),jstrt(MXTNSH),
     &            nstrt(MXTNSH),mst(MXTNSH),jrs(MXTNSH)
c-----------------------------------------------------------------------
CSS      Nshell = 2
      Write(6,"(a,I3)") "The Debug value for Nshell: ", Nshell
CSS      nprims(2) = 1

      Iloc   = 0
      Indecp = 0
      Do Icent = 1, Ntotatoms
      do ishell = 1, Nshell
CSSS      do ishell = 6, 6
        iatom  = katom(ishell)

        If (IGenby(Icent) .EQ. Iatom) Then

        xa = Coord(1,icent)
        ya = Coord(2,icent)
        za = Coord(3,icent)
        la = ktype(ishell) - 1

        kprimi=kprim(ishell)-1
        idegen=ndegen(ishell)
        numcoi=numcon(ishell)
        imin  =kmini(ishell)
        imax  =kmaxi(ishell)

C begin loop jshell 

        Jloc = 0
        Do Jcent = 1, Icent
        do jshell=1, Nshell
CSSS        do jshell=5, 5

          Jdegen=ndegen(Jshell)
          jatom=katom(jshell)
          lb = ktype(jshell) - 1
          jmin = kmini(Jshell)
          jmax = kmaxi(Jshell)
C
          If (IGenby(Jcent) .EQ. Jatom) Then
          If (.NOT. (ICent .EQ. Jcent .AND. JShell .GT. IShell)) Then
C
          xb=Coord(1,jcent)
          yb=Coord(2,jcent)
          zb=Coord(3,jcent)

          kprimj=kprim(jshell)-1
          numcoj=numcon(jshell)

          Call Dzero(Ecpint_4shell, Maxints_4shell)
C
          Do iprim = 1, nprims(ishell)
CSSS                  Do iprim = 3, 3

              exp1 = expnt(kprim(ishell)+iprim-1)
              Indnpri = Kprimi + Iprim
                     
              Do jprim = 1, nprims(jshell)
CSSS                     Do jprim = 3, 3
                 exp2 = expnt(kprim(jshell)+jprim-1)
                 Indnprj = Kprimj + Jprim
                 Call Dzero(Cint, Maxmem)
C
                 Call Ecp_int_4prim(Xa, Ya, Za, Xb, Yb, Zb, 
     &                              La, Lb, Coord, Exp1, Exp2,
     &                              Natoms, Ntotatoms, IGenby,
     &                              Cint, Int, Zero_int, 
     &                              Grads)

                 If (.NOT. Zero_int)
     &              Call ecp_int_4shell(Cint, Ecpint_4shell, 
     &                                  La, Lb, Numcoi, Numcoj, 
     &                                  Indnpri, Indnprj, Iprim,
     &                                  Jprim, Jnt)
             Enddo
          Enddo        
C 
      Write(6,*)
      Write(6,"(a,4(1x,I2))")
     &"The Contracted integral for shell pair", Ishell, JShell, Idegen,
     & Jdegen

      If (Ishell .NE. Jshell) Then
          Do I = 1, Idegen
             Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,Jdegen)
          Enddo
      Else
          Do I = 1, Idegen
             Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,I)
          Enddo
      Endif
                  Call Ecp_int(Ecpint, Ecpint_4shell, La, Lb, Iloc, 
     &                         Jloc, Numcoi, Numcoj, Imin, Imax, 
     &                         Jmin, Jmax, Indecp)
C
          Jloc = Jloc + (Lb+1)*(Lb+2)*Numcoj/2
C
        Endif
        Endif
C
        Enddo
        Enddo

        Iloc = Iloc + (La+1)*(La+2)*Numcoi/2
C
      Endif 
      Enddo 
      Enddo
C 
      Write(6,*)  Iloc*(Iloc+1)/2
      w=0.0d0
      Write(6,"(a)") "@-Ecpint_Driver, The ECP integrals"
      Write(6, "(6(1x,F10.7))") (Ecpint(I), I=1, Iloc*(Iloc+1)/2)
      do i =1, iloc*(iloc+1)/2
      W = w+ecpint(i)*ecpint(i)
      enddo
      Write(6, "(a,(1x,F20.13))") "The int. check sum = ", w
C
      Call Putrec(20, "JOBARC", "ECP1INTS", iloc*(iloc+1)/2, Ecpint)
C
      return
      end
