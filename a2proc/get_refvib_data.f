










      SUBROUTINE GET_REFVIB_DATA(Vcoords, Norm_coords, Vomega, Omega,
     &                           AtmMass, AtmLabel, SymLabel, Coords,   
     &                           Btmp, A2grad, Grad, Hess, Vhess,  
     &                           Imap, Vib_Type, Nreals, Natoms, 
     &                           Nvibs, B2ang, Au2Invcm, Trns_state,
     &                           Mass_Weigh_nm, Mass_Weigh_gr,
     &                           Get_hess, Get_Grad, Ivib_level,
     &                           Igrad_calc,Icol)

      Implicit Double Precision (A-H, O-Z)
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











































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




C
      Character*4 SymLabel(3*Natoms)
      Character*5 AtmLabel(Natoms)
      Character*11 Vib_Type(3*Natoms)
      Double Precision Norm_coords
C
      Logical Trns_state, Mass_weigh_nm, Get_Hess, Get_Grad, 
     &        Mass_Weigh_gr,Symmetry
      Logical HESSIANM_present
C
      Dimension VCoords(3*Natoms), Norm_coords(9*Nreals*Nreals), 
     &          Vomega(3*Nreals), Omega(3*Nreals), AtmMass(Natoms), 
     &          Coords(3*Nreals), Btmp(Natoms), Imap(Natoms),
     &          Hess(9*Natoms*Natoms), Vhess(9*Natoms*Natoms),
     &          Grad(3,Natoms), A2grad(3,Natoms)

      Data Thres /1.0D-09/

      Symmetry = (Iflags(60) .Gt. 0)
      Symmetry = .FALSE.
C The normal coordinates in NORMCORD are in M^(1/2)L dimensions.

      Call Getrec(20, "JOBARC", "REF_GEOM", Natoms*3*IINTFP,
     &            Vcoords)
      Call Getrec(20, "JOBARC", "NORMCORD", 3*Nreals*3*Nreals*
     &             IINTFP, Norm_coords)
      Write(6,"(2a)") "The reference (with dummy atoms geometry)",
     &                " (in Angstrom)"
      Write(6, "(3F10.5)") (B2ang*VCoords(i),i=1,3*Natoms)

C Remove the dummy atoms from the coordinates.

      Ioff = 1
      Joff = 1
      Ix   = 1
      Do Iatm = 1, Natoms
         Joff = Joff + (Iatm-1)*3
         If (AtmLabel(i)(1:5) .NE. 'X    ') Then
            Ioff = Ioff + (Ix-1)*3
            Call Dcopy(3,Vcoords(Joff),1,Vcoords(Ioff),1)
            Ix = Ix + 1
         Endif
      Enddo

C If projected frequencies are read, the VIB_TYPE record is 
C not in the correct order. Simply skip the first six frequencies.

      Call Getrec(20, "JOBARC", 'FORCECON', 3*Nreals*IINTFP,
     &            VOmega)
      Call Getcrec(20, "JOBARC", 'VIB_TYPE', 3*Nreals*11*IINTFP, 
     &             Vib_Type)

      write(6,*)
      Write(6,"(1x,2a)") " Vibration types (correspond to unprojected",
     +                   " Hessian)"
      Write(6,"(6(1x,a))") (Vib_type(i),i=1,3*Nreals)
C
      Ivibs = 0
      Do Imode = 1, 3*Nreals
C
         If (Vomega(Imode) .lt. 0.0D0) Then
            Vomega(Imode) = - DSQRT(DABS(Vomega(Imode)))*Au2Invcm
         Else
            Vomega(Imode) = DSQRT(Vomega(Imode))*Au2Invcm
         Endif
         If (Imode .Gt. Icol) Then
             Ivibs = Ivibs + 1
             Omega(Ivibs) = Vomega(Imode)
         Endif 
C
      Enddo
C
      If (Trns_state .AND. Omega(1) .GE. 0.0D0) Then
          Write(6, "(a)") "The IRC search must start from the TS."
CSSS          Call Errex
      Endif
C
      Call Getcrec(20, "JOBARC", "VIB_SYMS", 3*Nreals*4, SymLabel)
C
      Call Dcopy(Natoms, AtmMass, 1, Btmp, 1)
      Ireal = 0
      Icord = 1
      Jcord = 1
      Do Iatom = 1, Natoms
         Icord = Icord + (Iatom - 1)*3
         If (.Not. (Btmp(Iatom) .lt. 0.50D0)) Then
            Ireal = Ireal + 1
            Jcord = Jcord + (Ireal - 1)*3
            AtmMass(Ireal) = Btmp(Iatom)
            Call Dcopy(3, Vcoords(Icord), 1, Coords(Jcord), 1)
         Endif
         Icord = 1
         Jcord = 1
      Enddo
C
      Ivibs = 0
      Do Imode = 1, 3*Nreals
         If (Imode .Gt. Icol) Then
            Ivibs = Ivibs + 1
            SymLabel(Ivibs) = SymLabel(Imode)
         Endif
      Enddo
C
      Ivibs = 0
      Do Imodes = 1, 3*Nreals
         If (Imodes .Gt. Icol) Then
            Ivibs = Ivibs + 1
            Joff  = 3*Nreals*(Ivibs  - 1) + 1
            Ioff  = 3*Nreals*(Imodes - 1) + 1
            Call Dcopy(3*Nreals, Norm_coords(Ioff), 1,
     &                 Norm_coords(Joff), 1)
         Endif
      Enddo

         Write(6,*)  
         Write(6,"(2a)") "@get_ref_vibdata; Cartesian/normal mode",
     &                   " transformation from JOBARC (dimension L)" 
   
         Call output(Norm_coords,1,3*nreals,1,3*nreals,3*nreals,
     &               3*nreals,1)
C
      If (Mass_weigh_nm) Then 
          Call Get_massw_nrmlmodes(Norm_coords,AtmMass,Nvibs,Nreals)

      Endif 
C The Hessian that is read is exactly what enters into vib1.F in 
C joda. This is ordered as in ZMAT order and the dummy atoms have
C been removed.

      If (Get_Hess) then
          Length = 9*Nreals*Nreals 
          Call Getrec(20,'JOBARC','CARTHESC',Length*IINTFP,Hess)

         Write(6,*)  
         Write(6,*) "@get_ref_vibdata; Hessian from JOBARC"
         Call output(Hess,1,3*nreals,1,3*nreals,3*nreals,3*nreals,1)
      Endif 

C
C
      If (Get_Grad) Then

          Call Getrec(20,'JOBARC','GRADIENT',3*Nreals*IINTFP,
     &                A2Grad(1,1))
C
          Call Dzero(Grad, 3*Natoms)

C Remove dummy atoms 

          Do i=1, Natoms
             If (AtmLabel(i)(1:5) .NE. 'X    ') Then
                Do j=1,3
                   Grad(j,k)  = A2grad(j,i)
                Enddo
             Endif
          Enddo                  
C        
          Call Dcopy(3*Natoms, Grad, 1, A2grad, 1)

C Change to ZMAT order

          Do i=1, Natoms
             If (Imap(i) .NE. 0) Then
                 Call Dcopy(3, A2grad(1,Imap(i)), 1, Grad(1,i), 1)
             Endif
          Enddo
C
          If (Mass_Weigh_gr) Then 
             Do Iatom = 1, Nreals
                Do Ixyz = 1, 3
                   Grad(Ixyz, Iatom) = Grad(Ixyz, Iatom)/
     &                                 Dsqrt(AtmMass(Iatom))
                Enddo
             Enddo
          Endif
C
      Endif
C
      Return
      End

