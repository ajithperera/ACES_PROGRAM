      Subroutine Ecp_grdint_driver(Ntotatoms, IGenby, Coord, 
     &                             Cint, Dens_fao, Naobfns, Grad_xyz,
     &                             Grads)

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
     &          Ecpint(Mxcbf*(Mxcbf+1)/2),
     &          Ecpgrd_x(Maxints_4shell), Ecpgrd_y(Maxints_4shell), 
     &          Ecpgrd_z(Maxints_4shell), Dens_fao(Naobfns,*),
     &          Grad_xyz(3,Ntotatoms), Dens_4shell(Maxints_4shell),
     &          IGenby(Ntotatoms), Coord(3,Ntotatoms)

      common /infoa /natoms
C
      common /symshe/ mulsh(Mxtnsh)
      common /modez / zetm(Mxtnsh)
CSSS      common /pairij/ ipq(ndi4+1)
CSSS      common /typinf/ ltmax,ltdrv,nftnft,nftdrv
CSSS      common /powers/ jx(nftmax),jy(nftmax),jz(nftmax),
CSSS     &                ix(nftmax),iy(nftmax),iz(nftmax)

      common/xa/xand(3,8,4*(Maxang-1)-3)
CSSS      common/symsh2/mij,kij,hkij,kb
C
      integer and,eor,or,dstrt
C
c-----------------------------------------------------------------------
c     common from AcesII (created in Readin)
c-----------------------------------------------------------------------
      COMMON /NUCLEI/ CHARGE(MXATMS), CORD(MXATMS,3),
     &                DCORD(MXATMS,3),DCORGD(MXATMS,3),
     &                DOPERT(0:3*MXATMS)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(4,10), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXATMS,0:7), NCRREP(0:7),
     &                IPTCOR(MXATMS*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
c-----------------------------------------------------------------------
CSS      Nshell = 2
      Write(6,"(a,2I3)") "The Debug value for Nshell: ", Nshell,
     &                    Ntotatoms
      open(unit=99,status="new",form="formatted",file="ECPGRD")
CSS      nprims(2) = 1
      Iloc   = 0
      Indecp = 0

      Do Icent = 1, Ntotatoms 
      do ishell = 1, Nshell
CSSS      do ishell = 1, 1
C
        iatom  = katom(ishell)
        If (IGenby(Icent) .EQ. Iatom) Then

        xa = Coord(1,Icent)
        ya = Coord(2,Icent)
        za = Coord(3,Icent)
        la = ktype(ishell) - 1

        kprimi=kprim(ishell)-1
        idegen=ndegen(ishell)
        idegen_grd=Ideg_grd(la)
        numcoi=numcon(ishell)

C begin loop jshell 
 
        Jloc = 0
        Do Jcent = 1, Ntotatoms
        do jshell=1, Nshell
CSSS        do jshell=1, 1

C
C get symmetry infor regarding shell j
C
          jatom=katom(jshell)
          lb = ktype(jshell) - 1
          Jdegen=ndegen(jshell)
C
          If (IGenby(Jcent) .EQ. Jatom) Then
 
          xb=Coord(1,Jcent)
          yb=Coord(2,JCent)
          zb=Coord(3,Jcent)

          kprimj=kprim(jshell)-1
          numcoj=numcon(jshell)

C Get the density for this shell pair.
     
          Call Dzero(Dens_4shell, Maxints_4shell)
          Call Get_shell_den(Dens_fao, Dens_4shell, Naobfns, La, Lb,
     &                       Iloc, Jloc, Numcoi, Numcoj)

          Do iprim = 1, nprims(ishell)
CSSS          Do iprim =1, 1

             exp1 = expnt(kprim(ishell)+iprim-1)
             Indnpri = Kprimi + Iprim
                     
            Do jprim = 1, nprims(jshell)
CSSS            Do jprim = 1, 1 

CSS                Call Dzero(Ecpint_4shell, Maxints_4shell)
CSS                Call Dzero(Cint, Maxmem)

                exp2 = expnt(kprim(jshell)+jprim-1)
                Indnprj = Kprimj + Jprim
C
      Write(6, "(a, 4(1x,I2))") "The Ishell, Jshell, La, Lb: ", 
     &                           Ishell, Jshell, La, Lb
      Write(6, "(a, 4(1x,I2))") "The prim. pair; Iprim, Jprim",
     &                           Iprim, Jprim
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of A : ", Xa, ya, Za
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of B : ", Xb, yb, Zb
      Write(6, "(a, 2(1x,F10.6))") "The Primtive pair: ", Exp1, Exp2
      Write(6,*)
      Write(6,*) "Initilize debug data"
                Do Icnt = 1,  Ntotatoms
                   Isym_unq_cnt = IGenby(Icnt)
                   If (Ipseux(Isym_unq_cnt) .NE. 0) Then
                      
                      If (.NOT. (Grads .AND. Icent .EQ. Icnt)) Then
                         Iecp_cnt = Isym_unq_cnt
                         Call Dzero(Ecpint_4shell, Maxints_4shell)
                         Call Dzero(Cint, Maxmem)
 
                         Call Ecp_grdint_4prim(Xa, Ya, Za, Xb, Yb, 
     &                                         Zb, La, Lb, Coord, Exp1,
     &                                         Exp2, Natoms, Ntotatoms, 
     &                                         Cint, Int, Zero_int, 
     &                                         Icnt, Iecp_cnt, Grads)

      Write(6,"(a)") "Shell Densities"
      Write(6, "(6(1x,F10.7))") (dens_4shell(i), i=1,Numcoi*Numcoj*
     &                           Idegen*Jdegen)

                         If (.NOT. Zero_int)
     &                   Call ecp_int_4shell(Cint, Ecpint_4shell, 
     &                                       La, Lb, Numcoi, Numcoj, 
     &                                       Indnpri, Indnprj, Iprim,
     &                                       Jprim, Jnt, Grads)
      Write(6,*)
      Write(6,"(a,6(1x,I2))")
     &"The Contracted integral for shell pair", Ishell, JShell, Iprim,
     & Jprim, Idegen_grd, Jdegen
       Do I = 1, Idegen_grd
          Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,Jdegen)
      Enddo

                         Call Ecpgrd_int(Ecpint_4shell, La, Lb, 
     &                                   Idegen_grd, Idegen,
     &                                   Jdegen, Numcoi, Numcoj,
     &                                   Iprim, Jprim, Exp1, Exp2,
     &                                   Dens_4shell, Ecpgrd_x, 
     &                                   Ecpgrd_y, Ecpgrd_z, 
     &                                   Ntotatoms, Icent, ICnt, 
     &                                   Grad_xyz)
                      Endif
                   Endif 
C
                Enddo
C
             Enddo
          Enddo        

          Jloc = Jloc + (Lb+1)*(Lb+2)*Numcoj/2

        Endif 
        Enddo
        Enddo
        Iloc = Iloc + (La+1)*(La+2)*Numcoi/2

      Endif 
      Enddo 
      Enddo

      close(99)
C 
      return
      end
