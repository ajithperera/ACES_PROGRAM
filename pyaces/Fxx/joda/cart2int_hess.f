
C  *******************************************************
C  Conversion of the hessian from cartesian to red. internal
C  coordinates.
C  Luis Galiano, 07/07/03
C  *******************************************************

      SUBROUTINE CART2INT_HESS(FI,HC,TOTREDNCO,NRATMS,HI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER TOTREDNCO
c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end
c     Maximum string length of absolute file names
      INTEGER FNAMELEN
      PARAMETER (FNAMELEN=80)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/ XYZIN,NWFINDIF
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
c ric_heap.com : begin

c This common block contains the heap address and array indices for
c processing RICs. New additions must be initialized to 1 in
c bd_ric_heap.F and set in init_ric_heap.F.

      external bd_ric_heap

      double precision dRICHeap(1)
      integer          z_RICHeap, z_DerBMat, z_BMat, z_GMat, z_BTGInv

      common /ric_heap_com/ dRICHeap,
     &                      z_RICHeap, z_DerBMat, z_BMat, z_GMat,
     &                      z_BTGInv
      save   /ric_heap_com/

c ric_heap.com : end

C The following arrays are need to be dynamically allocated
C and deallocated.
C
C DIFTEMP: The is used in COMPKMAT as a temporary array to keep
C          an intermediate. It is of length 9*MXATMS*MXATMS.
C TEMP1  : Antoher temporary array to keep intermediates during
C          the transformations of the Cartesian Hessian. It is
C          of length 3*MAXREDUNCO*MXATMS.
C GMINBT : Keep the G-mat for Hessian transformations. It is of
C          length 9*MXATMS*MXATMS.
C GMINB  : To keep the  transpose of the G-mat, and it is of
C          the same lengthe as GMINBT.
C PFI    : Keep the gradients during the gradient projection.
C          It is of length MAXREDUNCO.
      DIMENSION DIFTEMP(9*MXATMS*MXATMS),FI(TOTREDNCO),PFI(MAXREDUNCO),
     &          HC(3*NRATMS,3*NRATMS),
     &          HI(TOTREDNCO*TOTREDNCO),
     &          GMINBT(3*MXATMS*MAXREDUNCO),
     &          TEMP1(MAXREDUNCO*3*MXATMS),
     &          GMINB(MAXREDUNCO*MAXREDUNCO)

      CALL GETREC(-1,'JOBARC','REDNCORD',1,TOTREDNCO)

C Build the K-matrix and the (H-K) matrix as shown in JCP, 117, 9160, 2002.
      CALL COMPKMAT(dRICHeap(z_DerBMat),NRATMS,TOTREDNCO,FI,HC,DIFTEMP)

      LENGTH_BGMAT =3*NRATMS*TOTREDNCO
      CALL GETREC(20,'JOBARC','BTGMIN',IINTFP*LENGTH_BGMAT,
     &            dRICHeap(z_BTGInv))
         print *, 'The BTGINV matrix'
         CALL OUTPUT(dRICHeap(z_BTGInv), 1, 3*NRATMS, 1, TOTREDNCO, 
     &               3*NRATMS, TOTREDNCO, 1)
      CALL GETREC(20,'JOBARC','GMATRIX',IINTFP*LENGTH_BGMAT,
     &            GMINBT)

         print *, 'The G matrix'
         CALL OUTPUT(GMINBT, 1, 3*NRATMS, 1, TOTREDNCO, 
     &               3*NRATMS, TOTREDNCO, 1)
C   Here we transpose the GMATRIX from JOBARC, because it was stored as
C   (G(-)*B)t and we need (G(-)*B)

      CALL TRANSP(GMINBT,GMINB,TOTREDNCO,3*NRATMS)

C   And here what we do is G(-)*B*(HC-K)*B(t)*G(-)

      CALL XGEMM('N','N',TOTREDNCO,3*NRATMS,3*NRATMS,1.0D0,
     &           GMINB,TOTREDNCO,DIFTEMP,3*NRATMS,0.0D0,
     &           TEMP1,TOTREDNCO)

      CALL XGEMM('N','N',TOTREDNCO,TOTREDNCO,3*NRATMS,1.0D0,
     &           TEMP1,TOTREDNCO,dRICHeap(z_BTGInv),3*NRATMS,0.0D0,
     &           HI,TOTREDNCO)

C Now we have the the Hessian in Redundant Internal Coordinates
         print *, 'Hessian in redundant internals'
         CALL OUTPUT(HI, 1, TOTREDNCO, 1, TOTREDNCO, TOTREDNCO,
     &               TOTREDNCO, 1)
      RETURN
      END

