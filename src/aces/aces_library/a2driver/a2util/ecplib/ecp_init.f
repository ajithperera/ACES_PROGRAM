










      Subroutine Ecp_Init

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
  


C -Note that (PRMTESTING block is used only in the early stages of
C development and no longer supported). Ajith Perera.

      Pi  = Dacos(-1.0D0)
      Fpi = 4.0D0*Pi
      Sqrt_Fpi = Dsqrt(4.0D0*PI)
      Sqpi2 = Dsqrt(Pi/2.0D0)

C Set lmnval array:

      II = 0
      DO LVAL = 0,MAXANG
         DO L = LVAL,0,-1
            LEFT = LVAL - L
            DO M = LEFT,0,-1
               N = LEFT - M
               II = II + 1
               LMNVAL(1,II) = L
               LMNVAL(2,II) = M
               LMNVAL(3,II) = N
             ENDDO
          ENDDO
      ENDDO
C
C Prepare I*(I+1)/2 values

      Ideg(0)= 1
      II = 1
      Do I = 1, Maxang
         II = II + 1
         Ideg(I) = II*(II+1)/2
      Enddo

C Prepare Istart and Iend arrays 

      Istart(0) = 1
      Istart(1) = 2
      Do I = 2, Maxang
         Istart(I) = Istart(I-1)+ Ideg(I-1)
      Enddo
C
      Iend(0) = 1
      Do I = 1, Maxang
         Iend(I) = Iend(I-1) + Ideg(I)
      Enddo
C
C Prepare I(I+1)/2 + (I-1)I/2 array for gradient

      Ideg_grd(0) = 3
      II = 1
      Do I = 1, Maxang
         II = II + 1
         Ideg_grd(I) =  (II+1)*(II+2)/2  + (II-1)*(II)/2
      Enddo

C Prepare Istart and Iend arrays for gradients

      Istart_grd(0) = 1
      Istart_grd(1) = 4
      Do I = 2, Maxang
         Istart_grd(I) = Istart_grd(I-1)+ Ideg_grd(I-1)
      Enddo

      Iend_grd(0) = 3
      Do I = 1, Maxang
         Iend_grd(I) = Iend_grd(I-1) + Ideg_grd(I)
      Enddo
      Call Izero(Lmnval_grd, 7*Lmnmaxg)
      Ioff = 1
      Do Iang = 0, Maxang
         
         Do Nlm = Istart(Iang), Iend(Iang)

            If (Lmnval(1, Nlm) .Gt. 0) Then 
               Nn = Lmnval(1, Nlm) - 1
               Ll = Lmnval(2, Nlm)
               Mm = Lmnval(3, Nlm)
               Ifound = 0
               Do J = Istart_grd(Iang), Ioff -1
                  If (Lmnval_grd(1,J) .EQ. Nn .AND.
     &                Lmnval_grd(2,J) .EQ. Ll .AND.
     &                Lmnval_grd(3,J) .EQ. Mm) Ifound = J
               Enddo

               If (Ifound .NE. 0) Then
                   Lmnval_grd(4, Ifound) =  Nlm - Istart(Iang) + 1
               Else
                   Lmnval_grd(7, Ioff) = -1
                   Lmnval_grd(1, Ioff) = Nn
                   Lmnval_grd(2, Ioff) = LL
                   Lmnval_grd(3, Ioff) = Mm
                   Lmnval_grd(4, Ioff) = Nlm - Istart(Iang) + 1
                   Ioff = Ioff + 1
               Endif
            Endif

            If (Lmnval(2, Nlm) .Gt. 0) Then 
               Nn = Lmnval(1, Nlm) 
               Ll = Lmnval(2, Nlm) - 1
               Mm = Lmnval(3, Nlm)
               Ifound = 0
               Do J = Istart_grd(Iang), Ioff -1
                  If (Lmnval_grd(1,J) .EQ. Nn .AND.
     &                Lmnval_grd(2,J) .EQ. Ll .AND.
     &                Lmnval_grd(3,J) .EQ. Mm) Ifound = J
               Enddo

               If (Ifound .NE. 0) Then
                   Lmnval_grd(5, Ifound) =  Nlm - Istart(Iang) + 1
               Else
                   Lmnval_grd(7, Ioff) = -1
                   Lmnval_grd(1, Ioff) = Nn
                   Lmnval_grd(2, Ioff) = LL
                   Lmnval_grd(3, Ioff) = Mm
                   Lmnval_grd(5, Ioff) =  Nlm - Istart(Iang) + 1
                   Ioff = Ioff + 1
               Endif
           Endif

           If (Lmnval(3, Nlm) .Gt. 0) Then 
               Nn = Lmnval(1, Nlm) 
               Ll = Lmnval(2, Nlm) 
               Mm = Lmnval(3, Nlm) - 1
               Ifound = 0
               Do J = Istart_grd(Iang), Ioff -1
                  If (Lmnval_grd(1,J) .EQ. Nn .AND.
     &                Lmnval_grd(2,J) .EQ. Ll .AND.
     &                Lmnval_grd(3,J) .EQ. Mm) Ifound = J
               Enddo

               If (Ifound .NE. 0) Then
                   Lmnval_grd(6, Ifound) = Nlm - Istart(Iang) + 1
               Else
                   Lmnval_grd(7, Ioff) = -1
                   Lmnval_grd(1, Ioff) = Nn
                   Lmnval_grd(2, Ioff) = LL
                   Lmnval_grd(3, Ioff) = Mm
                   Lmnval_grd(6, Ioff) = Nlm - Istart(Iang) + 1
                   Ioff = Ioff + 1
               Endif
            Endif
C
            Nn = Lmnval(1, Nlm) + 1
            Ll = Lmnval(2, Nlm)
            Mm = Lmnval(3, Nlm)
            Ifound = 0
            Do J = Istart_grd(Iang), Ioff -1
               If (Lmnval_grd(1,J) .EQ. Nn .AND.
     &             Lmnval_grd(2,J) .EQ. Ll .AND.
     &             Lmnval_grd(3,J) .EQ. Mm) Ifound = J
            Enddo
            
            If (Ifound .NE. 0) Then
               Lmnval_grd(4, Ifound) = Nlm - Istart(Iang) + 1
            Else
               Lmnval_grd(7, Ioff) = 1
               Lmnval_grd(1, Ioff) = Nn
               Lmnval_grd(2, Ioff) = LL
               Lmnval_grd(3, Ioff) = Mm
               Lmnval_grd(4, Ioff) = Nlm - Istart(Iang) + 1
               Ioff = Ioff + 1
            Endif

            Nn = Lmnval(1, Nlm) 
            Ll = Lmnval(2, Nlm) + 1
            Mm = Lmnval(3, Nlm)
            Ifound = 0
            
            Do J = Istart_grd(Iang), Ioff -1
               If (Lmnval_grd(1,J) .EQ. Nn .AND.
     &             Lmnval_grd(2,J) .EQ. Ll .AND.
     &             Lmnval_grd(3,J) .EQ. Mm) Ifound = J
            Enddo
            If (Ifound .NE. 0) Then
               Lmnval_grd(5, Ifound) = Nlm - Istart(Iang) + 1
            Else
               Lmnval_grd(7, Ioff) = 1
               Lmnval_grd(1, Ioff) = Nn
               Lmnval_grd(2, Ioff) = LL
               Lmnval_grd(3, Ioff) = Mm
               Lmnval_grd(5, Ioff) = Nlm - Istart(Iang) + 1
               Ioff = Ioff + 1
            Endif

            Nn = Lmnval(1, Nlm) 
            Ll = Lmnval(2, Nlm)
            Mm = Lmnval(3, Nlm) + 1
            Ifound = 0
            Do J = Istart_grd(Iang), Ioff -1
               If (Lmnval_grd(1,J) .EQ. Nn .AND.
     &             Lmnval_grd(2,J) .EQ. Ll .AND.
     &             Lmnval_grd(3,J) .EQ. Mm) Ifound = J
            Enddo
   
            If (Ifound .NE. 0) Then
               Lmnval_grd(6, Ifound) = Nlm - Istart(Iang) + 1
            Else
               Lmnval_grd(7, Ioff) = 1
               Lmnval_grd(1, Ioff) = Nn
               Lmnval_grd(2, Ioff) = LL
               Lmnval_grd(3, Ioff) = Mm    
               Lmnval_grd(6, Ioff) = Nlm - Istart(Iang) + 1
               Ioff = Ioff + 1
             Endif
C
         Enddo
      Enddo
C
      call Factorial
      call Dfactorial
      call Factorialo
      call Binomial_coefs
C
C Prepare Zlm, real spherical harmonic table and table 
C required to do the angular integral (Ftab). See the comments
C in the header of the individual subroutines for further details.
      
      Call Make_ztab
      Call Make_ftab

C -Note that (PRMTESTING block is used only in the early stages of
C development and no longer supported). Ajith Perera.


      call Factorial
      call Dfactorial
      call Factorialo
      call Binomial_coefs
C
C Prepare Zlm, real spherical harmonic table and table 
C required to do the angular integral (Ftab). See the comments
C in the header of the individual subroutines for further details.
      
      Call Make_ztab
      Call Make_ftab

C -Note that (PRMTESTING block is used only in the early stages of
C development and no longer supported). Ajith Perera.


      Return
      End












