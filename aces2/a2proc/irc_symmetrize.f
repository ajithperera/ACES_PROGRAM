










      Subroutine Irc_symmetrize(Xbohr, Xmw, Nreals, Natoms)
C
C Keeps an IRC point in the point group of the reference geometry. The
C points reuse the MOL file of the reference (transition state), so
C they must keep its symmetry exactly, but round off (and its growth
C in the constrained search) moved atoms off the symmetry elements:
C the HCOH points left the Cs plane by 1E-5 bohr and xvmol2ja could
C not map them (SIGSEGV in ISCTR). For each generator in the MOL file
C (X, Y, Z, XY, XZ, YZ, XYZ: the axes it inverts), an atom whose
C inverted coordinates are all zero in REF_GEOM lies on that symmetry
C element, and those coordinates are set to zero in Xbohr and Xmw
C (Nreals atoms, dummies left out, Z-matrix order as REF_GEOM). 10/2026.
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



C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
      Dimension Xbohr(3*Nreals), Xmw(3*Nreals), Ref(3*Mxatms),
     &          Amass(Mxatms)
      Character*1 Kka(9)
      Character*80 Fname
      Logical Onelt
      Parameter (Tiny = 1.0D-10)
C
      Call Gfname('MOL', Fname, Length)
      Open(Unit=31, File=Fname(1:Length), Status='Old', Err=99)
      Read(31,*)
      Read(31,*)
      Read(31,*)
      Read(31,"(2I5,9A1)", Err=98) Nontyp, Nsymop, (Kka(I), I=1,9)
      Close(31)
      If (Nsymop .Eq. 0) Return
C
      Call Getrec(20, "JOBARC", "REF_GEOM", Natoms*3*IINTFP, Ref)
      Call Getrec(20, "JOBARC", "ATOMMASS", Natoms*IINTFP, Amass)
C
      Ireal = 0
      Do Iatm = 1, Natoms
         If (Amass(Iatm) .Ge. 0.5D0) Then
            Ireal = Ireal + 1
            Do Igen = 1, Nsymop
C              On this generator's symmetry element?
               Onelt = .True.
               Do K = 1, 3
                  Iaxis = Index('XYZ', Kka(3*(Igen-1)+K))
                  If (Iaxis .Gt. 0) Then
                     If (Dabs(Ref(3*(Iatm-1)+Iaxis)) .Gt. Tiny)
     &                   Onelt = .False.
                  Endif
               Enddo
               If (Onelt) Then
                  Do K = 1, 3
                     Iaxis = Index('XYZ', Kka(3*(Igen-1)+K))
                     If (Iaxis .Gt. 0) Then
                        Xbohr(3*(Ireal-1)+Iaxis) = 0.0D0
                        Xmw(3*(Ireal-1)+Iaxis)   = 0.0D0
                     Endif
                  Enddo
               Endif
            Enddo
         Endif
      Enddo
      Return
C
 98   Close(31)
 99   Write(6,"(a)") " @IRC_SYMMETRIZE: cannot read the MOL file."
      Call Errex
      End
