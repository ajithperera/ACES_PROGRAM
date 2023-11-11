






































































































































































































      Subroutine Ecp_dgrad(Dens_tao, Naofns, Ncoord, Nbfirr, ECP_grd)
C
      Implicit Double Precision(A-h, O-Z)
C


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

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c
      Logical Spherical 

      Dimension Dens_tao(Naofns*(Naofns+1)/2), 
     &          Dens_fao(Mxcbf*Mxcbf), Grd_xyz(3,Mxatms),
     &          Mxpair(1), Map2zmat(Mxatms), Ecp_grd(3,Mxatms),
     &          Idexatom(Mxatms), Accord(3, Mxatms),
     &          Reorder_Mat(Mxcbf*Mxcbf), Tmpmat(Mxcbf*Mxcbf),
     &          Nbfirr(8)
    
      Common /Symminfo/ Ntotatoms, Idexatom, Acoord

      Call GetRec (20, "JOBARC", 'NIRREP  ', 1, Nirrep)
      Call Getrec(20, 'JOBARC','NATOMS  ', 1, Natoms) 
      call Getrec(20, 'JOBARC','MAP2ZMAT', Natoms, Map2zmat) 
      Call Getrec(20, 'JOBARC','NAOBASFN', 1, Naobfns)
      Call Getrec(20, "JOBARC", "NBASTOT ", 1, Nbfns)

      if (Naofns .Ne. Naobfns) then
         print *, '@REDUCE: Assertion failed.'
         print *, '   Nbas = ', Naofns
         print *, ' Naobf  = ', Naobfns
         call aces_exit(1)
      end if
      Call Dzero(Dens_fao, Naobfns*Naobfns)
      Call Sym_unpack(Dens_fao,Dens_tao,Naobfns,Nirrep,Nbfirr)


      Spherical = (iFlags(62).Eq.1)
      If (Spherical) Then
         Call Getrec(20, "JOBARC", "CART2CMP", Nbfns*Naobfns*Iintfp,
     &                Reorder_Mat)
         Call Xgemm("N", "N", Nbfns, Naobfns, Naobfns, 1.0D0, 
     &               Reorder_Mat, Nbfns, Dens_fao, Naobfns, 
     &               0.0D0, Tmpmat, Nbfns)
         Call Xgemm("N", "T", Nbfns, Nbfns, Naobfns, 1.0D0, Tmpmat, 
     &               Nbfns, Reorder_Mat, Nbfns, 0.0D0, Dens_fao, 
     &               Naobfns)
      Endif 

      Write(6,"(a)")"After transformation to spherical"
      call debugout(Dens_fao, nbfns)
  
      Call Getrec(20, "JOBARC", "CMP2ZMAT", Nbfns*Naobfns*Iintfp,
     &            Reorder_Mat)

      Call Xgemm("N", "N", Naobfns, Nbfns, Nbfns, 1.0D0, Reorder_Mat, 
     &            Naobfns, Dens_fao, Naobfns, 0.0D0, Tmpmat, Naobfns)
      Call Xgemm("N", "T", Naobfns, Naobfns, Nbfns, 1.0D0, Tmpmat, 
     &            Naobfns, Reorder_Mat, Naobfns, 0.0D0, Dens_fao, 
     &            Naobfns)

      Write(6,"(a)")"The reordered AO basis density matrix"
      call debugout(Dens_fao, nbfns)
      Ntotatoms = Natoms
      call Getrec(20, 'JOBARC','COORD   ', 3*Ntotatoms*Iintfp, 
     &            Acoord) 
      Call Dzero(Grd_xyz, 3*Mxatms)
C
      Call ecp_grdint_main(Ntotatoms, Dens_fao, Naobfns, Grd_xyz)
C
      Write(6,*)
      Write(6,"(a)") "Gradients returned from ecp_grdint_main"
      Write(6,"(3(1x,F10.7))") ((Grd_xyz(i,j),i=1,3),j=1,natoms)

      Do Iatms = 1, Natoms

         If (Map2zmat(Iatms) .NE. 0) Then

            Do Ixyz = 1, 3
           
               Ecp_grd(Ixyz, Iatms) = Grd_xyz(Ixyz, Map2zmat(Iatms))*
     &                                2.0D0
            Enddo

         Endif
      Enddo

      Write(6,*)
      Write(6, "(a)") "The gradient due scalar-ECP part"
      Do iatms = 1, Natoms
         Write(6, "(3(1x, F15.12))") (Ecp_grd(Ixyz,iatms), ixyz=1,3)
      Enddo 

      Call Tran_c2sgrad(Ecp_grd, Grd_xyz, Ncoord)
      Call Dcopy(Ncoord, Grd_xyz, 1, Ecp_grd, 1)

      Write(6,*)
      Write(6, "(a)") "The sym. ada. gradient due scalar-ECP part"
      Do iatms = 1, Natoms
         Write(6, "(3(1x, F15.12))") (Ecp_grd(Ixyz,iatms), ixyz=1,3)
      Enddo 
      Return
      End
