





































































































































































































      Subroutine ECP_MAIN(Iecp,Natoms,Namat,Iqmstr,Jcostr,Nucstr,
     &                    Nrcstr)
C      
      Implicit Double Precision(A-H, O-Z)


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
  

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end

      Character*2 mtype(Mxatms)
      Character*80 Cnick(Mxatms)
      Character*80 Ecpnam(Mxatms)
      Character*4 namat(Mxatms)
      Logical Nharm, Grads
C
      Dimension iqmstr(Mxatms),jcostr(Mxatms,Maxang), 
     &          nucstr(Mxatms,Maxang,Mxprim),
     &          nrcstr(Mxatms,Maxang,Mxprim),IGenby(Mxatms),
     &          Coord(3,Mxatms)
C
      Dimension  Cint(Maxmem)
C
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     & MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     & AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     & NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     & NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
C
      COMMON /DAT/  EXPA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     &              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
C
C ECP integrals. Except for simplifications, I have kept the original
C form of all the input/output ECP processing routines. This was 
C originaly written by Christien Huber, Ajith Perera, 12/2001.
C
C - Natoms is the number of symmetry unique atoms
C - Ipseux is 0 when no ECP definitions; > 0 otherwise.
C - Cnick is the name of then ECP 
C - iqmstr highest ang. mom. qun. num.(AMQN) for each atom.
C - jcostr has the AQMN for each shell (each atom). Normally the # of shells
C   is the same as the Max AMQN but some cases there may be mutiple shells
C   with the same AQMN, so jcostr tells you how many.
C - nucstr is the number of primitives per each shell (per atom,per AQMN).
C - nrcstr is the number of contracted functions per each shell
C   (per atom,per AQMN)
C - nhram is logical that tells you whether this is sperical or 
C   cartesian basis calculation (key-word driven)
C - I must add that the names choosen for those 4 variables are dumb!!
C   
      nharm = (Iflags(62) .EQ. 1)
      Itol  = Iflags(108)
      Thres = Dble(10.0D+00*(-Itol))
 
      if (iecp.ne.0) then
c-----------------------------------------------------------------------
c get atomsymbols and put them to mtype-array
c-----------------------------------------------------------------------
         ierr=0
         call chrges(CHARGE,mtype,natoms,ierr)
c-----------------------------------------------------------------------
c read ecp data and prepare data for evulation of the c integrals
c-----------------------------------------------------------------------

         if (ierr.eq.0) then
             ipr=0
             call ecpdef(natoms,mtype,Ecpnam,cnick)

             if (ierr.eq.0) then
                ipr=0
                call ecppar(ipr,cnick,CHARGE,natoms,.FALSE.,0)
C
c-----------------------------------------------------------------------
c  Summation over the charges and storage of this value
c-----------------------------------------------------------------------

                chgsum=0
                do i=1,natoms
C
C account for symmetry 
C
                    chgsum=chgsum+charge(i)*fmult(mulnuc(i))
                Enddo

                nproton=idint(chgsum)
c
                call putrec(20,'JOBARC','NMPROTON',1,nproton)

c-----------------------------------------------------------------------
c prepare data for ecp integral evaluation
c-----------------------------------------------------------------------
                Call ecp_init
                Call ecp_basis_init(namat,iqmstr,jcostr,nucstr,
     &                              nrcstr, NHARM)
          
c-----------------------------------------------------------------------
             else
                write(6,"(a)") 'ERROR while reading ECP data!'
                write(6,"(a)") 'Please check yourECPDATA file.'
                Call Errex
             endif

         else
             write(6,"(a)") 'E R R O R in subroutine chrges!'
             Call Errex
         endif
c-----------------------------------------------------------------------

      Endif

      Call Dzero(Cint, Maxmem)
      Call Getrec(20, 'JOBARC', 'NATOMS  ', 1, Ntotatoms) 
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Ntotatoms, Coord)

      Call Symeqv(Ntotatoms, IGenby)
      Call ECP_int_driver(Natoms, Ntotatoms, IGenby, Coord, Cint, 
     &                    .False.)
      
      Return
      End
