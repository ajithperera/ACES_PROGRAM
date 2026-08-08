










      Subroutine Gen_curr_coords(React, Prodt, Currn, Ndim, XYZ)

      Implicit Double Precision (A-H,O-Z)
 
      Logical XYZ

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)


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



C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end



      Dimension React(Ndim), Prodt(Ndim), Currn(Ndim)

      Call Getrec(20, "JOBARC", "RXSTRUCT", Ndim*IINTFP, React)
      Call Getrec(20, "JOBARC", "PRSTRUCT", Ndim*IINTFP, Prodt)

      Write(6,"(a)")"The React/Product structues for LST optimizations"
      Write(6,"(3(1x,F12.6))") (React(i),i=1, Ndim)
      Write(6,*)
      Write(6,"(3(1x,F12.6))") (Prodt(i),i=1, Ndim)
CSSS      Call Dcopy(Ndim, React, 1, Currn, 1)

      If (XYZ) Then
      
         Do Iatm = 1, Ndim/3
            Do jatm = 1, Iatm
               Distij = Dist(React(3*(Iatm-1)+1), React(3*(Jatm-1)+1))
               Distkl = Dist(Prodt(3*(Iatm-1)+1), Prodt(3*(Jatm-1)+1))
               
               Distmn = Distij + (Distkl - Distij)/2.0D0
      Write(6,"(a)") "Distij, Distkl, Distmn"
      Write(6,"(3(1x,F12.6))") Distij, Distkl, Distmn*0.529177249D0
            Enddo
         Enddo
C
         Call Putrec(1,'JOBARC','CUSTRUCT',Ndim*iintfp,Currn)
      
      Else
  
         Do idim = 1, Ndim
 
CSSS            Currn(Idim) = React(Idim) + (Prodt(Idim) - 
CSSS     &                                   React(Idim))/2.0D0
         Enddo
 
CSSS         Call Dcopy(Ndim, Currn, 1, R, 1)
         Call Putrec(20, "JOBARC", "RXSTRUCT", Ndim*IINTFP, Currn)

      Endif
C
      Write(6,"(a)") "The current structue for LST optimizations"
      Write(6,"(3(1x,F12.6))") (currn(i),i=1, Ndim)
      Return
      End
