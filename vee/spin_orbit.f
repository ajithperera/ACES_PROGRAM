










      Subroutine Spin_orbit(Tdens, Work, Maxcor, NBASIS)

      Implicit Double Precision (A-H, O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)

      Dimension Tdens(Nbasis, Nbasis), Work(Maxcor), Socc(3)
      Character*8 String(3)
      Character*80 Fname 
      Dimension Zeff(18)
      Logical Vpout_exist 
      Integer Charge(Mxatms), Atom_number 



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




       Data String /'   OPX  ','   OPY  ','   OPZ  '/

       Data Zeff /1.0D0,0.0D0,1.35D0,2.0D0,2.75D0,3.60D0,
     +            4.55D0,5.60D0,6.75D0,0.0D0,10.04D0,10.80D0,
     +            11.54D0,12.25D0,12.94D0,13.60D0,14.24D0,0.D0/

C Read the 1-el spin-orbit integral from VPOUT file

      CALL GFNAME('VPOUT', FNAME, ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=VPOUT_EXIST)

      If (Vpout_exist) Then
         OPEN (UNIT=30, FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     +        STATUS='OLD')
      Else
         Write(6,"(a)") " For 1el perturbative spin-orbit effects",
     +                  " calculations VPOUT file is needed."
         Call Errex
      Endif 
C
      Irwnd = 0
      Nsize = Nbasis * Nbasis 

      CALL GETREC (20, 'JOBARC', 'NREALATM', 1, NATOMS)
      CALL GETREC (20, 'JOBARC', 'ATOMCHRG', NATOMS,CHARGE)

      DO Iatom = 1, NATOMS 

         DO Ixyz = 1, 3
            I000 = 1
            I010 = I000 + Nsize * IINTFP

            CALL SEEKLB (String(Ixyz), IERR, IRWND, 30)
            IF (IERR .NE. 0) CALL ERREX

            CALL LOADINT (Work(I000), NATOM, NSIZE, NAO, IUHF)
            
            Atom_number = Charge(Iatom) 
            Z_effctive  = Zeff(Atom_number)
   
            Socc_bare  = Ddot(Nsize, WOrk(I000), 1, Tdens, 1)
            Socc(Ixyz) = Socc_bare * Z_effctive + Socc(Ixyz)
     +                 + Socc(Ixyz) 

            Irwnd = 1
         Enddo 
       Enddo

       Write(6,"(a)") " 1el perturbative effective spin-orbit effects"
       write(6,"(3(a,F12.7))") "SOX = ", Socc(1), "SOY = ", Socc(2), 
     +                      "SOZ = ", Socc(3)

      Return
      End

  
