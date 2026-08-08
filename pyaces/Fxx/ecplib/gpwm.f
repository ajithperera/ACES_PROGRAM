      Subroutine Gpwm(Nlq, Lprj, La, Lb, Ltot, Lamalo, Lamahi, Lamblo, 
     &                Lambhi, Lamau, Lambu, Alpha, Rc, Rc2, Prd, Beta1, 
     &                Beta2, Rad2)
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

      
      Dimension Rad2(0:2*Maxang,0:2*Maxang,0:2*maxang), Z(35), W(35),
     &          Bessa(20,0:Maxang), Bessb(20,0:Maxang), 
     &          Pntw(20,0:Maxang)

C
      data              z /-.20201828704561d+01,-.95857246461382d+00,
     & .00000000000000d+00, .95857246461382d+00, .20201828704561d+01,
     &                     -.34361591188377d+01,-.25327316742328d+01,
     &-.17566836492999d+01,-.10366108297895d+01,-.34290132722370d+00,
     & .34290132722370d+00, .10366108297895d+01, .17566836492999d+01,
     & .25327316742328d+01, .34361591188377d+01,
     &                     -.53874808900112d+01,-.46036824495507d+01,
     &-.39447640401156d+01,-.33478545673832d+01,-.27888060584281d+01,
     &-.22549740020893d+01,-.17385377121166d+01,-.12340762153953d+01,
     &-.73747372854539d+00,-.24534070830090d+00, .24534070830090d+00,
     & .73747372854539d+00, .12340762153953d+01, .17385377121166d+01,
     & .22549740020893d+01, .27888060584281d+01, .33478545673832d+01,
     & .39447640401156d+01, .46036824495507d+01, .53874808900112d+01/
      data              w / .19953242059046d-01, .39361932315224d+00,
     & .94530872048294d+00, .39361932315224d+00, .19953242059046d-01,
     &                      .76404328552326d-05, .13436457467812d-02,
     & .33874394455481d-01, .24013861108231d+00, .61086263373533d+00,
     & .61086263373533d+00, .24013861108231d+00, .33874394455481d-01,
     & .13436457467812d-02, .76404328552326d-05, 
     &                      .22293936455342d-12, .43993409922732d-09,
     & .10860693707693d-06, .78025564785321d-05, .22833863601635d-03,
     & .32437733422379d-02, .24810520887464d-01, .10901720602002d+00,
     & .28667550536283d+00, .46224366960061d+00, .46224366960061d+00,
     & .28667550536283d+00, .10901720602002d+00, .24810520887464d-01,
     & .32437733422379d-02, .22833863601635d-03, .78025564785321d-05,
     & .10860693707693d-06, .43993409922732d-09, .22293936455342d-12/

      Write(6,"(a,9(1x,I3))") "Nlp,lprj,La,Lb,ltot,lamalo,lamahi,
     &lamblo,lambhi: ",Nlq,lprj,La,Lb,ltot,lamalo,lamahi,lamblo,
     &                lambhi
      Write(6,*)
      TRc2  = Rc2


      If (Trc2 .Gt. 5.0D4) Then
         Npnts = 5
         Idifs = 0 
      Else  if (Trc2 .Gt. 5.0D2) Then
         Npnts = 10
         Idifs = 5
      Else
         Npnts = 20
         Idifs = 15
      Endif

      Sqalpha = Dsqrt(Alpha)
      Rcalpha = Rc
      Prd     = Prd/Sqalpha

      Do Ipnts = 1,  Npnts
        
         Pnt = Z(Ipnts + Idifs)/Sqalpha + Rcalpha 
         Do Lama = Lamalo, Lamahi
            Bessa(Ipnts, Lama) = Bess(Beta1*Pnt, Lama)

         Enddo
        
         
         Do Lamb = Lamblo, Lambhi
            Bessb(Ipnts, Lamb) = Bess(Beta2*Pnt, Lamb)
         Enddo

         If (Nlq .Gt. 0) Then
   
            Pntw(Ipnts, 0) = Prd*pnt**Nlq

            Do L = 1, Ltot
        
               Pntw(Ipnts, L) = pnt*Pntw(Ipnts, L-1)

            Enddo

         Else

            Pntw(Ipnts, 0) = Prd

            Do L = 1, Ltot
        
               Pntw(Ipnts, L) = pnt*Pntw(Ipnts, L-1)

            Enddo

         Endif
       
         Do Lama = Lamalo, Lamahi

            Na_min = Iabs(Lprj-Lama) 

            Do Lamb = Lamblo, Lambhi

               Nb_min = Iabs(Lprj-Lamb)

               Nlo = Na_min + Nb_min
               Nhi = (Ltot - Mod(La-Na_min,2)) - Mod(Lb-Nb_min,2)

               Do N = Nlo, Nhi, 2


                   Rad2(N, Lamb, Lama) = Rad2(N, Lamb, Lama) + 
     &                                   ((W(Ipnts+Idifs)*
     &                                    Bessa(Ipnts, Lama))*
     &                                    Bessb(Ipnts, Lamb))*
     &                                    Pntw(Ipnts, N)
               Enddo
            Enddo
         Enddo

      Enddo

      Return
      End

