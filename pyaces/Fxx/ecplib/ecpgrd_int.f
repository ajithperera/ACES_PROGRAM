      Subroutine Ecpgrd_int(Ecpint_4shell, La, Lb, Idegen_grd, Idegen, 
     &                      Jdegen, Numcoi, Numcoj, Iprim, Jprim, 
     &                      Exp1, EXp2, Dens_4shell, Ecpgrd_x, Ecpgrd_y,
     &                      Ecpgrd_z, Natoms, Iatom, Jatom, Grad_xyz)

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
  

     
      Dimension Ecpint_4shell(*), Ecpgrd_x(*), Ecpgrd_y(*), 
     &          Ecpgrd_z(*), Dens_4shell(*), Grad_xyz(3, Natoms)
      
      Length = Idegen*Jdegen*Numcoi*Numcoj
      Call Dzero(Ecpgrd_x, Length)
      Call Dzero(Ecpgrd_y, Length)
      Call Dzero(Ecpgrd_z, Length)

      Texp1 = 2.0D0*Exp1
      Texp2 = 2.0D0*Exp2
      
      La_lo = Istart_grd(La)
      La_hi = Iend_grd(La)
      
      Lb_lo = Istart_grd(Lb)
      Lb_hi = Iend_grd(Lb) 

      Do Ndegen_grd = 1, Idegen_grd

         Ix = Lmnval_grd(4, Ndegen_grd + La_lo - 1)
         Iy = Lmnval_grd(5, Ndegen_grd + La_lo - 1)
         Iz = Lmnval_grd(6, Ndegen_grd + La_lo - 1)
         
         Tax = Texp1
         Tay = Texp1
         Taz = Texp1

         If (Lmnval_grd(7, Ndegen_grd + La_lo - 1) .Lt. 0) Then

            Tax = -(Lmnval_grd(1, Ndegen_grd + La_lo - 1) +1)
            Tay = -(Lmnval_grd(2, Ndegen_grd + La_lo - 1) +1)
            Taz = -(Lmnval_grd(3, Ndegen_grd + La_lo - 1) +1)

         Endif


         Ndim = Jdegen*Numcoi*Numcoj

         If (Ix .NE. 0) Then
            Ioff = (Ndegen_grd-1)*Ndim+1
            Joff = (Ix-1)*Ndim+1
            Call daxpy(Ndim, Tax, Ecpint_4shell(Ioff), 1, 
     &                 Ecpgrd_x(Joff), 1)
         Endif


         If (Iy .NE. 0) Then
            Ioff = (Ndegen_grd-1)*Ndim+1
            Joff = (Iy-1)*Ndim+1
            Call daxpy(Ndim, Tay, Ecpint_4shell(Ioff), 1, 
     &                 Ecpgrd_y(Joff), 1)
         Endif
             
         If (Iz .NE. 0) Then
            Ioff = (Ndegen_grd-1)*Ndim+1
            Joff = (Iz-1)*Ndim+1
            Call daxpy(Ndim, Taz, Ecpint_4shell(Ioff), 1, 
     &                 Ecpgrd_z(Joff), 1)
         Endif

      Enddo

      Write(6,*)
      Write(6,"(a)") "The differentiated integrals"
      Write(6,"(a,4(1x,I2))")
      Write(6,"(a)") "d/dx integrals"
      Write(6,"(4(1x,F20.13))")(Ecpgrd_x(j), J=1,Length)

      Write(6,"(a)") "d/dy integrals"
      Write(6,"(4(1x,F20.13))")(Ecpgrd_y(j), J=1,Length)

      Write(6,"(a)") "d/dz integrals"
      Write(6,"(4(1x,F20.13))")(Ecpgrd_z(j), J=1,Length)

      Write(6,*)
      Write(6,"(a)") "Shell Densities"
      Write(6, "(4(1x,F20.13))") (dens_4shell(i), i=1, length)

      Grad_x = Ddot(Length, Ecpgrd_x, 1, Dens_4shell, 1)
      Grad_y = Ddot(Length, Ecpgrd_y, 1, Dens_4shell, 1)
      Grad_z = Ddot(Length, Ecpgrd_z, 1, Dens_4shell, 1)

      Write(6,*) 
      Write(6,"(a)") "Grad_xyz"
      Write(6,"(3(1x,F20.13))") Grad_x, Grad_y, Grad_z
      Write(99, '(2I4)') Iatom, Jatom
      Write(99,"(3(1x,F20.13))") Grad_x, Grad_y, Grad_z
      
      Grad_xyz(1, Iatom) = Grad_x + Grad_xyz(1, Iatom)
      Grad_xyz(2, Iatom) = Grad_y + Grad_xyz(2, Iatom)
      Grad_xyz(3, Iatom) = Grad_z + Grad_xyz(3, Iatom) 

      Grad_xyz(1, Jatom) = -Grad_x + Grad_xyz(1, Jatom)
      Grad_xyz(2, Jatom) = -Grad_y + Grad_xyz(2, Jatom)
      Grad_xyz(3, Jatom) = -Grad_z + Grad_xyz(3, Jatom) 

      Return
      End
