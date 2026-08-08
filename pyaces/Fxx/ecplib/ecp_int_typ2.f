      Subroutine Ecp_int_typ2(Iecp_cnt, Lamau, Lambu, La, Lb, Exp1,
     &                        Exp2, Exp12, Xahat, Yahat, Zahat, Xbhat,
     &                        Ybhat, Zbhat, Xp, Yp, Zp, Xc, Yc, Zc, 
     &                        Acs_xyz, Bcs_xyz, Ca, Cb, Pc2, Fact_ab,
     &                        CA_Zero, CB_zero, Cint, Rad2_zero, 
     &                        Grads)

      Implicit Double Precision (A-H, O-Z)
      Logical Rad2_zero, CA_zero, CB_zero, Grads


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
  

   
      Dimension ACS_xyz(0:2*Maxang,3),BCS_xyz(0:2*Maxang,3)
      Dimension Rad2(0:2*Maxang,0:2*Maxang, 0:2*Maxang),
     &          Ang2_A(0:2*Maxang,-Maxproj:Maxproj,0:2*Maxang),
     &          Ang2_B(0:2*Maxang,-Maxproj:Maxproj,0:2*Maxang),
     &          Cint(Maxmem)

      Lmax  = LLMax(Iecp_cnt)
      If (Grads) Then
          Lap = La + 1
         Ltot = La + Lb + 1
      Else
          Lap = La
         Ltot = La + Lb 
      Endif
C     
C
C This handles U(l)-U(lmax), l = 0,...,lmax-1 
C
      Test = ((exp1*exp2)/exp12)*(Ca-Cb)**2
      If (Test .Gt. Tol) Then
         Rad2_zero = .True.
         Return
      Endif
           
      If (.Not. CA_zero) Then
          Lamau = Lamau + Lmax - 1
      Else
CSSS          Lamau = 1
          Lamau = 0
      Endif
      If (.Not. CB_zero) Then
          Lambu = Lambu + Lmax - 1
      Else
CSSS          Lambu = 1
          Lambu = 0
      Endif

      If (CA_zero .And. CB_zero) then
          Lprjhi = Min(Lmax-1, Lap+1, Lb+1)
          Lprjlo = Mod(Lap, 2)
          Inc    = 2
C
          If (Lprjlo .Ne. Mod(Lb,2) .Or. Lprjlo .Gt. Lprjhi) Then
             Rad2_zero = .True.
             Return
          Endif

      Elseif (CA_zero) Then
          Lprjhi = Min(Lmax-1, Lap+1)
          Lprjlo = Mod(Lap, 2)
          Inc    = 2

          If (Lprjlo .Gt. Lprjhi) Then
             Rad2_zero = .True.
             Return
          Endif
C
      Elseif (CB_zero) Then
          Lprjhi = Min(Lmax-1, Lb+1)
          Lprjlo = Mod(Lb, 2)
          Inc    = 2

          If (Lprjlo .Gt. Lprjhi) Then
             Rad2_zero = .True.
             Return
          Endif
C
      Else     
          Lprjhi = Lmax-1 
          Lprjlo = 0
          Inc    = 1
      Endif
C          
C

      Do Lprj = Lprjlo, Lprjhi, Inc

         Lamalo = Max(Lprj-Lap, 0)
         Lamahi = Min(Lprj+Lap, Lamau)
         Lamblo = Max(Lprj-Lb, 0)
         Lambhi = Min(Lprj+Lb, Lambu)

        Call Dzero(Rad2, (2*Maxang+1)**3)
	Kcrl = Kfirst(Lprj+2, Iecp_cnt)
        Kcru = Klast(Lprj+2, Iecp_cnt)

         Do Kcr = Kcrl, kcru
        
            Zeta = Zlp(Kcr)
            Np   = Nlp(Kcr)
            Dlj  = Clp(Kcr)

C
            Call Ecp_rad_int_typ2(Lamalo, Lamahi, Lamblo, Lambhi,
     &                            Lamau, Lambu, Lap, Lb, Ltot, Lprj, 
     &                            Np, Exp1, Exp2, Exp12, Zeta, Dlj, 
     &                            Xahat, Yahat, Zahat, Xbhat, Ybhat, 
     &                            Zbhat, Xp, Yp, Zp, Xc, Yc, Zc, CA,
     &                            CB, Pc2, Fact_ab, Rad2, Rad2_zero)
         Enddo
C
           If (.Not. Rad2_zero) 
     &     call Ecp_ang_int_typ2(Xahat, Yahat, Zahat, Xbhat, Ybhat, 
     &                           Zbhat, La, Lb, Ltot, Lamalo, Lamahi, 
     &                           Lamblo, Lambhi, Lamau, Lambu, Lprj, 
     &                           Acs_xyz, Bcs_xyz, Ang2_A, Ang2_B,
     &                           Rad2, Cint, Int, Grads)

      Enddo
      
CSSS      Call Dscal(Int, Fpi, Cint, 1) 

C 
      Return
      End
