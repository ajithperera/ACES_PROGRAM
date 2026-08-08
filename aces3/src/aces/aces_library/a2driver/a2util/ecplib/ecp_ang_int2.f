










      Subroutine Ecp_Ang_int2(Xph, Yph, Zph, Lamu, lab, Nx, Ny, Nz, 
     &                        Lprj, PQS_xyz, Ang)

      Implicit Double Precision (A-H, O-Z)
      Logical Skip


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
  


      Dimension PQS_Xyz(0:2*Maxang,3),TPQ_Xyz(0:2*Maxang,3),
     &          Ang(0:2*Maxang,-Maxproj:Maxproj,0:2*Maxang)


       Call make_BC_PQ_xyz(Nx, PQS_xyz(0,1), TPQ_xyz(0,1))
       Call make_BC_PQ_xyz(Ny, PQS_xyz(0,2), TPQ_xyz(0,2))
       Call make_BC_PQ_xyz(Nz, PQS_xyz(0,3), TPQ_xyz(0,3))
C
C loop over Mprj (Mprj= -Lprj:Lprj)
       Do Iang = 0, 2*Maxang
          Do iprj = -Maxproj, Maxproj
              Do Jang = 0, 2*Maxang
                 Ang(Jang, Iprj, Iang) = 0.0D0
              Enddo
          Enddo
      Enddo
C
      Do Mprj = -Lprj, Lprj
C
c Now we have loops over Nax, Nay, NaZ
c
         Do ix = 0, Nx
            Do iy = 0, Ny
               Do iz = 0, Nz
C
                  Ixyz = ix+iy+iz

                  l_high = Min(Lprj+Ixyz, Lamu)

                  Do Lam = 0, Lamu

c Loops over mu (-la to +la and -lb to +lb)

                     Ang_N_lam = 0.0D0

                     Do mu = -Lam, Lam
C
                        If (Lam .LE. Ixyz+Lprj .and. 
     &                      Mod(Lam+Lprj+Ixyz,2) .NE. 1) Then

CSSS                           Write(6,*) Lam, Ixyz+Lprj, Lam+Lprj+Ixyz
                            Call Screen_ang_typ2(Lprj, MPrj, Ix, Iy,
     &                                           Iz, Lam, Mu, Skip)
CSSS                            If (.NOT. Skip) Then

                               Call Make_res(Lam, Mu, Xph, Yph, Zph,
     &                                       Res)
C
                               Call Ang_int(Lprj, Mprj, Ix, Iy, Iz, 
     &                                      Lam, Mu, Ang_N_lam_mu)
                               Ang_N_lam = Ang_N_lam + Ang_N_lam_mu*
     &                                     Res*Sqrt_fpi
C
CSSS                           Endif
                        Endif
C
                     Enddo

                     If (Lam .LE. Ixyz+Lprj .and. Mod(Lam+Lprj+Ixyz,2)
     &                   .NE. 1) Then
CSSS                        Write(6,*) Lam, Ixyz+Lprj, Lam+Lprj+Ixyz

                         Ang(Ixyz, Mprj, Lam) = Ang(Ixyz, Mprj, Lam) +
     &                                          Ang_N_lam*
     &                                          TPQ_xyz(Ix,1)*
     &                                          TPQ_xyz(Iy,2)*
     &                                          TPQ_xyz(Iz,3)
C
                     Endif

                  Enddo

                Enddo
            Enddo
         Enddo
C
      Enddo

      Return
      End

