










      Subroutine Radint_recur(Nlq, Ltot, Qnl)
  
      Implicit Double Precision (A-H, O-Z)
C

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
      Dimension Qnl(0:2*Maxang, 0:2*Maxang)
      Logical Eve, Odd
      Data A0, Ahalf, A1, A2, A3 /0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0/

      write(6,*) "@-Radint_recur the Qstore common block values"
      Write(6,"(a,3(1x,F18.13),2(1x,i2))")"Alpha,Beta,Xval,Nlq,Ltot: ", 
     &                      Alpha,Beta,Xval, Nlq, Ltot
      Write(6,*) 
C
      Odd = ((Mod(Nlq,2)+1) .EQ. 2)
      Eve = ((Mod(Nlq,2)+1) .EQ. 1)

C
      Talpha = 2.0D0*Alpha
C
C Compute the radial integral, Int {0,inf} exp(-alpha r^n) M_l(beta r)
C employing recursion relations given in Mcmurchie and Davidson,
c Journal of Computational Physics, 44, 289 (1981).
C
      If (Beta .Eq. 0.0D0) then

C Using standard integral for Bessel functions, Page 302, Text book
C of Mathematical functions (Abramowitz and Stegun, Eds. 1964, US
C BUreau of Standards).

         If (Odd) Then

            If (Ltot .LE. 1) Then

                If (Nlq .EQ. 1) Then

                    Qnl(0,0) = Ahalf*(A1/Alpha)

                Else

                    Qnl(0,0) =  Ahalf*(A1/Alpha)
                       Term1 = (A1/Alpha)
                    
                    Do NN = 2, Nlq, 2
                        
                        Qnl(0,0) = Qnl(0,0)*Term1
                           Term1 = Term1 + (A1/Alpha)
                    
                    Enddo
                
                Endif
            
            Else
                
                Do L = 0, Ltot, 2
                   
                   NN = (Nlq + L)/2
                   Qnl(L,0) = Fact(NN)/(A2*(Alpha**(NN+1)))
                
                Enddo
            
            Endif
         
         Else
             
             If (Ltot .LE. 1) Then
                    
                 If (Nlq .EQ.0) Then
                    
                    Qnl(0,0) =  Ahalf*Dsqrt(Pi/Alpha)
                 
                 Else
                     
                     Qnl(0,0) =  Ahalf*Dsqrt(Pi/Alpha)
                        Term1 = Ahalf*(A1/Alpha)

                     Do NN = 2, Nlq, 2

                        Qnl(0,0) = Qnl(0,0)*Term1
                           Term1 = Term1 + (A1/Alpha)

                     Enddo

                 Endif

             Else

                 Do L = 0, Ltot, 2

                    NN = (Nlq + L)/2
                    Qnl(L,0) = Faco(NN)/((A2**(NN+1))*(Alpha**NN))*
     &                        (Dsqrt(Pi/Alpha)) 

                 Enddo

             Endif

         Endif

      Else

         If (Ltot .LE. 1) Then
            Do L = 0, Ltot
               Do N = 0, L
C
C Only (N-L) positive and even is needed.
C
                   If ((N-L) .GE. 0) Qnl(N,L) = Qcomp(N+Nlq, L)

               Enddo
            Enddo
         Endif

         If (Ltot .LE. 1) Return
C
         If (Nlq .GT. 0) Then

            If (Eve) Then

C Npl > 0 and even
C
C N = 0+2
               Qnl(0,0) = Qcomp(2, 0)

               Do LL = 0, Ltot-1
C N = LL+3-2
                  N = LL+1
                  L = LL+1
                  Qnl(N,L) = Beta*(Qnl(N-1,L-1))/Talpha

               Enddo

               If ((Nlq-2) .EQ. 0) Then
C
C Use recursion 38D of Davidson paper.
C
                   Lmax = Ltot
                   Do NN = 0, Ltot-2, 2

                      Lmax = Lmax - 2
                      Do L= 0, Lmax

C NLM = N+2-1; N=0+2, 0+4,...

                         N = L + 2 + NN
                         NLM = N + L + 1
                         Qnl(N,L) = ((NLM)*Qnl(N-2,L) +
     &                                Beta*Qnl(N-1,L+1))/Talpha

                      Enddo
                   Enddo

               Else
C
C Use recursion 38D in MD paper (multiple times).
C
                  Call Qcomp_38D(Qnl, Nlq, Ltot, Talpha)

               Endif

            Else

               If (Xval .GT. 3.0D0) Then
C
C Npl > 0, and odd  Xval > 3.0
C N = 0+1, 1+1...
                   Qnl(0,0) = Qcomp(1,0)
                   Qnl(1,1) = Qcomp(2,1)
C
C Use recursion 38E in MD paper.
C
                   Do LL = 3, Ltot+1

                      L = LL - 1
                      N = LL - 1
                      Qnl(N,L) =  ((N+L-2)*Qnl(N-2,L-2) +
     &                             (Beta-(2*L-1)*(Talpha/Beta))*
     &                              Qnl(N-1,L-1))/Talpha
                   Enddo
C
C Use recursion 38D in MD paper (multiple times).
C
                   Call Qcomp_38D_odd(Qnl, Nlq, Ltot, Talpha)

               Else
C
C Npl > 0, and odd  Xval < 3.0
C N = Ltot+1, Ltot+1+....

                  Qnl(Ltot-1,Ltot-1)   = Qcomp(Ltot, Ltot-1)
                  Qnl(Ltot,Ltot) = Qcomp(Ltot+1, Ltot)
C
C Using 38F from MD paper.
C
                  Do LL = (Ltot-1), 1, -1
                     L = LL - 1
                     N = LL - 1

                     Qnl(N,L) = (Talpha*Qnl(N+2,L+2) - (Beta-(2*L+3)*
     &                           (Talpha/Beta))*Qnl(N+1,L+1))/(N+L+2)
                  Enddo
C
C Use recursion 38D in MD paper (multiple times).
C
                  Call Qcomp_38D_Odd(Qnl, Nlq, Ltot, Talpha)

               Endif

            Endif

         Else
C
C Nlq=0, Ltot >=2 bblock
C
            If (Ltot .EQ. 2) Then
C
C Using prescription A in Fig. 1 of MD paper
c
               Qnl(2,0) = Qcomp(2,0)

               If (Xval .Gt. 3.0D0) Then

                  Qnl(0,0) = Qcomp(0,0)
C
C Using 38B of MD paper
C
                  Qnl(1,1) = (Talpha*Qnl(2,0) - Qnl(0,0))/Beta
C
C Combining 38A and 38B of MD paper,note the cancellation of constants.
C
CSSSS Qnl(3,1) = Beta*Q(2,0)/Talpha
C
                  Qnl(2,2) = Qnl(2,0) - 3.0D0*Qnl(1,1)/Beta

               Else
C
C Using prescription C in Fig. 1 and combining with 38C and 38A.
c
                  Qnl(2,2) = Qcomp(2,2)

CSSS Qnl(1,1) = Beta*Qnl(2,0)/Talpha
CSSS Qnl(3,1) = Beta*Qnl(2,0)/Talpha

                  Qnl(1,1) = Beta*(Qnl(2,0) - Qnl(2,2))/3.0D0
                  Qnl(0,0) = (Talpha*Qnl(2,0) - Beta*Qnl(1,1))

               Endif

            Else
C
C Nlq=0, Ltot>2 block
C Using recursion 38A of MD paper.
C
               Qnl(2,0) = Qcomp(2,0)

               Do L = 0, Ltot - 3

                  Qnl(3+L,L+1) = Beta*(Qnl(L+2,L))/Talpha

               Enddo

               If (Xval .Gt. 3.0D0) Then

                   Qnl(0,0) = Qcomp(0,0)

                   Do L = 1, Ltot - 1
C
C Using 38B of MD paper.
C
                      Qnl(L,L) = (Talpha*Qnl(L+1,L-1) - (2*L-1)*
     &                            Qnl(L-1,L-1))/Beta
                   Enddo
C
C Combining 38A and 38B of MD paper,note the cancellation of constants.
C
                   Qnl(Ltot,Ltot) = Qnl(Ltot,Ltot-2) - (2*Ltot-1)*
     &                              Qnl(Ltot-1, Ltot-1)/Beta

C Using 38D of MD paper (Note the error in MD paper)
C
                   Lmax = Ltot - 2
                   Do NN = 0, Ltot-4, 2
                      Lmax = Lmax - 2
                      Do L= 0, Lmax

                         N = L + 4 + NN
                         Qnl(N,L) = ((N+L-1)*Qnl(N-2,L) +
     &                                Beta*Qnl(N-1,L+1))/Talpha

                      Enddo
                   Enddo

               Else
C
C Using prescription C in Fig. 1 and combining with 38C and 38A.
c
                  Qnl(2,2) = Qcomp(2,2)

CSSS         Qnl(1,1) = Beta*Qnl(2,0)/Talpha
CSSS         Qnl(3,1) = Beta*Qnl(2,0)/Talpha

                  Qnl(1,1) = Beta*(Qnl(2,0) - Qnl(2,2))/3.0D0

                  Do L = Ltot-3, 0, -1

                     Qnl(L,L) = (Talpha*Qnl(L+2,L) - Beta*Qnl(L+1,L+1)/
     &                          (2*L+1))
                  Enddo
C
C Using 38D of MD paper (Note the error in MD paper),
C
                   Lmax = Ltot - 2
                   Do NN = 0, Ltot-4, 2
                      Lmax = Lmax - 2
                      Do L= 0, Lmax

                         N = L + 4 + NN
                         Qnl(N,L) = ((N+L-1)*Qnl(N-2,L) +
     &                                Beta*Qnl(N-1,L+1))/Talpha

                      Enddo
                   Enddo

C
               Endif

            Endif
C
         Endif
C
      Endif


      Return
      End


                              
