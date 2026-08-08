










      subroutine Ecp_Grdint_4prim(Xa, Ya, Za, Xb, Yb, Zb, La, Lb, 
     &                            Coord, Exp1, Exp2, Natoms, Ntotatoms, 
     &                            Cint, Int, Zero_int, Icnt, Iecp_cnt, 
     &                            Grads)

C (XA,YA,ZA) and (XB,YB.ZB) are the coordintes of center A and B.
C LA and LB are the angular momentums of gaussians on A and B (0,...n)
C Cont1 and Cont1 are the contractions coefficients for gaussian A and B.
C Exp1 and Exp1 are the exponents of gaussian A and B
C
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
  


      Dimension Coord(3,Mxatms), ACS_xyz(0:2*Maxang,3),
     &          BCS_xyz(0:2*Maxang,3), Cint(Maxmem)
      Logical CA_zero, CB_zero, PC_zero, Zero_int, Rad1_zero,
     &        Rad2_zero, Grads, Zero_int_cnt(Mxatms)
C
C 
C Loop over all the ECP centers and compute the integral 
C for each pair of  primitives

      Ltot = La + Lb + 1
      Lap = La + 1

      Xc = Coord(1, Icnt)
      Yc = Coord(2, Icnt)
      Zc = Coord(3, Icnt)

C
C Build CA anb CB,(= C-A, C-B) vector components
C
      Xca = Xc - Xa
      Yca = Yc - Ya
      Zca = Zc - Za
      Ca2 = Xca**2 + Yca**2 + Zca**2

      if (Ca2 .Gt. 0.0D0) Then
         Ca = Dsqrt(Ca2)
         Lamau = Lap
         Xahat = -Xca/Ca
         Yahat = -Yca/Ca
         Zahat = -Zca/Ca
         CA_zero = .False.
      Else 
         Ca = 0.0D0
         Lamau = 0
         Xahat = 0
         Yahat = 0
         Zahat = 0
         CA_zero = .True.
      Endif
C
      Xcb = Xc - Xb
      Ycb = Yc - Yb
      Zcb = Zc - Zb
      Cb2 = Xcb**2 + Ycb**2 + Zcb**2

      if (Cb2 .Gt. 0.0D0) Then
          Cb = Dsqrt(Cb2)
          Lambu = Lb
          Xbhat = -Xcb/Cb
          Ybhat = -Ycb/Cb
          Zbhat = -Zcb/Cb
          CB_zero = .False.
      Else
          Cb = 0.0D0
          Lambu = 0
          Xbhat = 0
          Ybhat = 0
          Zbhat = 0
          CB_zero = .True.
      Endif
C
C Build PIn_i=n,,k_i=0...n_i}CA_i^^(n_i-k_i) i=x,y,z
C

      Call Make_xyz(Lap, Xca, ACS_xyz(0,1))
      Call Make_xyz(Lap, Yca, ACS_xyz(0,2))
      Call Make_xyz(Lap, Zca, ACS_xyz(0,3))

      Call Make_xyz(Lb, Xcb, BCS_xyz(0,1)) 
      Call Make_xyz(Lb, Ycb, BCS_xyz(0,2))
      Call Make_xyz(Lb, Zcb, BCS_xyz(0,3))


C
C Get the common center of gaussian on A and B and the new
C multification factor (Fact_AB).
C
      Call Get_center(Xa, Ya, Za, Xb, Yb, Zb, Exp1, Exp2,
     &                Xp, Yp, Zp, Exp12, Fact_Ab)

      Xpc = Xp - Xc
      Ypc = Yp - Yc
      Zpc = Zp - Zc
       
      Pc2 = Xpc**2 + Ypc**2 + Zpc**2

      If (Pc2 .Gt. 0.0D0) Then
         Pc    = Dsqrt(Pc2)
         Lamu  = Ltot
         Xhat  = Xpc/pc
         Yhat  = Ypc/pc
         Zhat  = Zpc/pc
         Beta  = 2.0D0*Exp12*PC
         PC_Zero = .False.
      Else 
          Pc    = 0.0D0
          Lamu  = 0
          Beta  = 0.0D0 
          Xhat  = 0.0D0
          Yhat  = 0.0D0
          ZhAt  = 0.0D0
          PC_Zero = .True.
      Endif

      Call Ecp_int_typ1(Iecp_cnt, Lamu, Ltot, La, Lb, Exp12, 
     &                  Xhat, Yhat, Zhat, Xp, Yp, Zp, Xc, Yc,
     &                  Zc, Acs_xyz, Bcs_xyz, Pc2, Fact_Ab,
     &                  PC_zero, Cint, Int, Rad1_zero, Grads)
      Call Ecp_int_typ2(Iecp_cnt, Lamau, Lambu, La, Lb, 
     &                  Exp1, Exp2, EXp12, Xahat, Yahat, Zahat, 
     &                  Xbhat, Ybhat, Zbhat, Xp, Yp, Zp, Xc, Yc,
     &                  Zc, Acs_xyz, Bcs_xyz, Ca, Cb, Pc2, 
     &                  Fact_Ab, CA_Zero, CB_Zero, Cint, 
     &                  Rad2_zero, Grads)

      Zero_int_cnt(Icnt) = Rad1_zero .AND. Rad2_zero 


      Zero_int = .True. 
      Do Kcnt = 1, Natoms
         If (.NOT. Zero_int_cnt(Kcnt)) Then
             Zero_int = .False.
             Return
         Endif
      Enddo 
C
C Do the contractions. 
C
      Return
      End
