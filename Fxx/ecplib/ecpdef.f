      SUBROUTINE ECPDEF(NATOMS,MTYPE,Ecpnam,CNICK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c     NOTE : THE ATOMIC SYMBOL must PRECEED THE ECP NICK NAME
c
c     ---------------------------------------------------------------
c     MEANING OF VARIABLES :
c     ---------------------------------------------------------------
c
c     --------- INPUT -----------------------------------------------
c     MTYPE   : ATOMIC IDENTIFIER (STARTING WITH ELEMENT SYMBOL)
c     NATOMS  : NUMBER OF ATOMS
c     --------- OUTPUT ----------------------------------------------
c     IPSEUX  : PSEUDO POTENTIAL TYPE OF ATOM
c               0  MEANS : NO PSEUDO POTENTIAL OR NOT YET ASSIGNED
c     CNICK   : PSEUDO POTENTIAL NICKNAMES INDEXED BY PSEUDO POTENTIAL
c               TYPES
c
      CHARACTER*16 BLNKBN
      PARAMETER (BLNKBN = '                ')

      CHARACTER ATOSYM*2,MTYPE(NATOMS)*(*)
      CHARACTER*4 DUMPF
      CHARACTER*80 ECPNAM(NATOMS)
      CHARACTER*80 CNICK(NATOMS)
      LOGICAL bDOIT



c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end




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
  

c
      CALL GETREC(20,'JOBARC','NATOLD',1,NATOLD)
      CALL GETREC(20,'JOBARC','ECPNAM',80*NATOLD,ECPNAM)
 
      do i=1,natold
      Write(6,*) ECPNAM(i)
      enddo
c  
c first detect different types of ecp names
C
 
      DO I=1,NATOMS
         IPSEUX(I)=0
      END DO

      DO I=1,NATOMS
         CNICK(I)=BLNKBN
      END DO
C
c ECPNAM is the full list of ECP names per atom including dummy atoms
C
      NOFECP=0

      DO IAT=1, NATOLD

c loop over the atoms and remove all redundant and NONE lines
         IF ((ECPNAM(IAT).NE.BLNKBN).AND.
     &       (INDEX(ECPNAM(IAT),'NONE').EQ.0)) THEN
            bDOIT=.TRUE.
            DO II=1,NOFECP
               IF (ECPNAM(IAT).EQ.CNICK(II)) bDOIT=.FALSE.
            END DO
            IF (bDOIT) THEN
               NOFECP=NOFECP+1
               CNICK(NOFECP)=ECPNAM(IAT)
            END IF
         END IF
      END DO
      ICHECK=0
      DO IPOT=1,NOFECP
         NPOS=INDEX(CNICK(IPOT),':')

         IF (NPOS.EQ.0) THEN
            WRITE(*,*) '@ECPDEF: Error - missing seperator (:) ',
     &                 'between atom symbol and ECP nickname.'
            CALL ERREX
         END IF
c
c extract the full atomic symbol ( C != CL )
c
         DUMPF = '  '
         DUMPF = CNICK(IPOT)(1:NPOS-1)

c loop over the symmetry unique atoms and assign matches to the 
c current CNICK line number

         DO IAT=1,NATOMS

            IF (DUMPF(1:2).EQ.MTYPE(IAT)(1:2)) THEN

c bomb if two orbits have the same atom (CNICK only contains 
C unique ECP definitions)

               if (IPSEUX(IAT).ne.0) then
                  print *, '@ECPDEF: ERROR - Two orbits have the same ',
     &                     'atom and different ECP definitions.'
                  print *, '         The current implementation cannot',
     &                     ' distinguish between the two.'
                  call errex
               end if
               IPSEUX(IAT)=IPOT
               ICHECK=ICHECK+1
            END IF
         END DO
      END DO

      IF (ICHECK.EQ.0) THEN
         IF (NOFECP.EQ.0) THEN
            WRITE(*,*) '@ECPDEF: Error - no ECP nicknames found in ',
     &                 'JOBARC.'
            CALL ERREX
         END IF
         WRITE(*,*) '@ECPDEF: Error - ECP nickname and atom symbols ',
     &              'in ZMAT are incompatible.'
         CALL ERREX
      END IF

      RETURN 
      END

