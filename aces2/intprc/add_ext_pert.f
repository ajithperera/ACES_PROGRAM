














































































































































































































      Subroutine  Add_ext_pert(Work, Maxcor, X_field, Y_field, Z_field,
     &                         F, Nbas, Nbasx, Ispin)

      Implicit Double Precision (A-H, O-Z)
 
      Dimension Work(Maxcor), F(Nbas, Nbas)
      Logical X_field, Y_field, Z_field 
      Character*8 Label(3), DIRECTION



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





      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
C      
      Length = Nbas*(Nbas+1)/2
C   
      IXYZt_efld = 1
      IXYZf_efld = IXYZt_efld + Length
      IAO_efld   = IXYZf_efld + Length
      IMO_efld   = IAO_efld   + Nbas*Nbas
      IEVEC      = IMO_efld   + Nbas*Nbas
      IScr       = IEVEC      + Nbas*Nbasx
      Inext      = Iscr       + 2*Nbax*Nbasx
C
      Call ZERO(Work(IXYZf_efld), Length)

      IF (X_field) Then
         IFLD_STREN  = Iflags(23)
         DXFLD_STREN = DFLOAT(IFLD_STREN)*1.0D-06
         DIRECTION   = '     X  '
         WRITE(6,100) DIRECTION(6:6) , DXFLD_STREN

         Call Getrec(20, "JOBARC", Label(1), Length*Iintfp,
     &                         Work(IXYZt_efld))
         CALL DSCAL(Length, DXFLD_STREN, Work(IXYZt_efld), 1)
         Call Dcopy(Length, Work(IXYZt_efld), 1,  Work(IXYZf_efld), 1)
      Endif

      IF (Y_field) Then
         IFLD_STREN = Iflags(23)
         DXFLD_sTREN = DFLOAT(IFLD_STREN)*1.0D-06
         DIRECTION   = '     Y  '
         WRITE(6,100) DIRECTION(6:6), DYFLD_STREN

         Call Getrec(20, "JOBARC", Label(2), Length*Iintfp,
     &               Work(IXYZt_efld))
         Call Daxpy(Length, DYFLD_sTREN, Work(IXYZt_efld), 1, 
     &              Work(IXYZf_efld), 1)

      Endif 

      IF (Z_field) Then
          IFLD_STREN = Iflags(25)
          DZFLD_STREN = DFLOAT(IFLD_STREN)*1.0D-06
          DIRECTION   = '     Z  '
          WRITE(6,100) DIRECTION(6:6), DZFLD_STREN

         Call Getrec(20, "JOBARC", Label(3), Length*Iintfp,
     &               Work(IXYZt_efld))
         Call Daxpy(Length, DZFLD_sTREN, Work(IXYZt_efld), 1, 
     &              Work(IXYZf_efld), 1)

      Endif 

      Call Expnd2(Work(IXYZf_efld), Work(IAO_efld), Nbas)
      Call AO2MO2(Work(IAO_efld), Work(IMO_efld), Work(Ievec),
     &                 Work(Iscr), Nbas, Nbasx, Ispin)

      Call Daxpy(Nbas*Nbas, 1.0D0, Work(IMO_efld), 1, F, 1)


100   FORMAT(T3,'@FINFLD-I, Adding electric field perturbation to ',
     &          'Fock matrix.',/,T10,'Field direction',
     &          T30,': ',A1,/,T10,'Field strength',T30,': ',F10.7,
     &          ' a.u.')

       Return
       End 
