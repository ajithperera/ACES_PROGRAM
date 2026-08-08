










      Subroutine Ecp_rad_int_typ2(Lamalo, Lamahi, Lamblo, Lambhi,
     &                            Lamau, Lambu, La, Lb, Ltot, Lprj,  
     &                            NP, Exp1, Exp2, Exp12, Zeta, Dlj, 
     &                            Xahat, Yahat, Zahat, Xbhat, Ybhat,
     &                            Zbhat, Xp, Yp, Zp, Xc,
     &                            Yc, Zc, CA, CB, Pc2, Fact_ab, Rad2,
     &                            Rad2_zero)

      Implicit Double Precision (A-H, O-Z)
      Logical Rad2_zero
       

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
      Dimension Rad2(0:2*Maxang,0:2*Maxang, 0:2*Maxang)

      Alpha = Exp12 + Zeta

      Beta1 = 2.0D0*Ca*Exp1
      Beta2 = 2.0D0*Cb*Exp2

      Talpha = 2.0D0*Alpha
      Falpha = 4.0D0*Alpha
      
      Rc  = (Beta1 + Beta2)/Talpha 
      Rc2 = (Beta1 + Beta2)**2/Falpha

C The commented block of code is how I would have computed the
C exponential factor. However, the one below is presented as
C simplifaction in IJQC, 40, 773., 1991. I do not see how. If
C I do not use what is on above paper, I can not match with 
C older codes. 
CSSS      Call Get_center(Xp, Yp, Zp, Xc, Yc, Zc, Exp12, Zeta,
CSSS     &                xq, Yq, Zq, Exp123, Fact_Pc)
CSSS
CSSS      Exp_fac = Fact_ab + Fact_pc

      Fact_ab = ((Exp1*Exp2)/(Exp12))*(Ca-Cb)**2
      Fact_pc = (Zeta*Alpha)*(Rc**2)/(Exp12)
      Exp_fac = Fact_ab + Fact_pc
      Rad2_zero = .True. 


      If (Exp_fac .Le. Tol) Then
      
         Rad2_zero = .False.

         Pre_fac = Dlj*Exp_fac

CSSS         If (Beta1 .EQ. 0 .AND. Beta2 .EQ. 0) 
CSSS            Beta = 0
CSSS            Xval = 0
CSSS            Rad2(Ltot,1,1) = Rad2(Ltot,1,1) + Qcomp(NP+Ltot,0)
CSSS         Elseif (Beta1 .EQ. 0) Then
CSSS            Beta = Beta2
CSSS            Xval = Rc2
CSSS            Do Lamb = Lprj, Lambhi
CSSS               Rad2(Ltot,Lamb,1) = Rad2(Ltot,1,1) + Qcomp(NP+Ltot,0)
CSSS            Enddo
         If (Rc2 .LT. 50.0D0) Then
C
C Use the single powers series with recrusions (Eqns. 51 and 52).

            Call Sps(Np, Lprj, La, Lb, Ltot, Lamalo, Lamahi, Lamblo, 
     &               Lambhi, Lamau, Lambu, Alpha, Beta1, Beta2, 
     &               Dlj, Exp_fac+Rc2, Rad2)
C 
         Else

C Use the Gaussian quadrature (points and weights)
 
             Prd = Dexp(-Exp_fac)*Dlj
            
             Trc1 = Max(2.0D0*Rc*Beta1,1.0D0)
             Trc2 = Max(2.0D0*Rc*Beta2,1.0D0)
             
             Fexp1 = (4.0D0*Exp1)**(La+1)
             Fexp2 = (4.0D0*Exp2)**(Lb+1)
             Taiaj = (4.0D0*Exp1*Exp2)
            
             Term1 = Dabs(Prd)/(Trc1*Trc2)
             Term2 = Dsqrt(4.0D0*Fexp1*Fexp2*Dsqrt(Taiaj)/Alpha)
         
             Ecp_lim = Term1*Term2
             If (Rc .Ge. Ca) Then
                Nplim = Np + La
                If (Rc .Ge. CB) Then
                    Nplim = Np + Lb
                Else
                    Eclim = Ecp_lim*Cb**Lb
                Endif
             Else
                Nplim = Np
                Eclim = Ecp_lim*Ca**la
                If (Rc .Ge. CB) Then
                    Nplim = Np + Lb
                Else
                    Eclim = Ecp_lim*Cb**Lb
                Endif
             Endif
             
             Aprx_Ecp = Eclim*Rc**Nplim

             If (Aprx_Ecp .Gt. Eps1) Call Gpwm(Np, Lprj, La, Lb, Ltot, 
     &                                         Lamalo, Lamahi, Lamblo, 
     &                                         Lambhi, Lamau, Lambu, 
     &                                         Alpha, Rc, Rc2, Prd, 
     &                                         Beta1, Beta2, Rad2)


         Endif
     
      Endif 
        
      Return
      End
