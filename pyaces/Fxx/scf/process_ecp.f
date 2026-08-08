










      Subroutine process_ecp(Ecpint, Onehamil, Ldim, Nbas, Nirrep, 
     &                       Nbfirr)
C
      Implicit Double Precision (A-H, O-Z)



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
      Dimension Ecpint(Ldim), Onehamil(Ldim), Ao2so(Mxcbf*Mxcbf),
     &          Tmpmat1(Mxcbf*Mxcbf), Tmpmat2(Mxcbf*Mxcbf),
     &          Nbfirr(Nirrep)
 
      Call Getrec(20, 'JOBARC','NAOBASFN', 1, Naobfns)
      Length = Naobfns*(Naobfns + 1)/2
      
      Call Getrec(20, "JOBARC", "ECP1INTS", Length*Iintfp, Tmpmat2)
      Call Expnd2(Tmpmat2, Tmpmat1, Naobfns)


      Call Getrec (20, "JOBARC", "CMP2ZMAT", Nbas*Naobfns*Iintfp,
     &             Ao2So)

      
      Call Xgemm("N", "N", Naobfns, Nbas, Naobfns, 1.0D0, Tmpmat1,
     &            Naobfns, Ao2So, Naobfns, 0.0D0, Tmpmat2, Naobfns)

      Call Xgemm("T", "N", Nbas, Nbas, Naobfns, 1.0D0, Ao2so,
     &            Naobfns, Tmpmat2, Naobfns, 0.0D0, Tmpmat1, Naobfns)

      Call sympack(Tmpmat1, Ecpint, Naobfns, Nirrep, Nbfirr)


      Call Daxpy(Ldim, 1.0D0, Ecpint, 1, Onehamil, 1)

    
      Return
      End
     
             
    
