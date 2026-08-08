














































































































































































































      Program  pccd_drpmo
      Implicit none
      Integer iuhf
      Integer IGrad_calc
      Integer IHess_calc
      Integer Pccd_refr,Itest
      Integer OOmicroItr,OOmacroItr,OOtotalItr
      Integer Length, Junk
      Double Precision Scale 
      Character*80 Fname
      Character*6 Status
      Logical OOexist 
      DOUBLE PRECISION GLOBNORM
      INTEGER MACRONUM,MICRONUM,TOTALITR
      LOGICAL MACCONVG,MICCONVG
      Logical Opt_orbs, OO_constr, OV_constr 
      Logical Pccd,Ccd,Lccd,Symmetry 

c COMMON BLOCKS
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end


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




      Common /OO_info/Opt_orbs,OO_constr,OV_constr
      Common /CALC/PCCD,CCD,LCCD
      Common /Symm/Symmetry

c ----------------------------------------------------------------------

      Call Aces_init(icore, i0,icrsiz, iuhf, .true.)
      Pccd  = .FALSE.
      Ccd   = .FALSE.
      Lccd  = .FALSE.
      CALL Getrec(0,"JOBARC","PCCD_RUN",Length,Junk)
      If (Length .Gt. 0) pCCD = .True.
      If (IFLAGS(2).EQ.53) Lccd = .True.
      If (IFLAGS(2).EQ.54) Ccd  = .True. 

      IF(IFLAGS(60).EQ.2) Symmetry  = .TRUE.

      IGrad_calc = 1
      IHess_calc = 1
      Status     = "Begin "
      
      Call Gfname("ooinfo.dat",Fname,Length)
      Inquire(file=Fname(1:Length),exist=OOexist)
      If (OOexist) Then
         Open(Unit=1011,File=Fname(1:Length),Status="old")
         Read(1011,*) MACRONUM,MICRONUM,TOTALITR,MACCONVG,
     +                          MICCONVG,GLOBNORM
         call putrec(-20,"JOBARC","MACRONUM",1,MACRONUM)
         call putrec(-20,"JOBARC","MICRONUM",1,MICRONUM)
         call putrec(-20,"JOBARC","TOTALITR",1,TOTALITR)
         CALL putrec(-20,"JOBARC","MACCONVG",1,.False.)
         CALL putrec(-20,"JOBARC","MICCONVG",1,.False.)
         call putrec(-20,"JOBARC","GLOBNORM",1,0.0d0)
         close(1011, status='delete')
      Endif

      call getrec(-20,"JOBARC","MICRONUM",1,OOmicroItr)
      call getrec(-20,"JOBARC","MACRONUM",1,OOmacroItr)
      call getrec(-20,"JOBARC","TOTALITR",1,OOtotalItr)

      If (Iflags2(138) .Eq. 2) Then
          IGrad_calc = 0
          IHess_calc = 0
          Scale     = Iflags(57)
      Endif
      Opt_orbs  = (Iflags2(178) .Gt. 0)
      Oo_constr = (Iflags2(181) .Eq. 1)
      Ov_constr = (Iflags2(184) .Eq. 1)
C
      If (Igrad_calc .eq. 1 .and. IHess_calc .eq.1) Then

C If Bruekner optimization is sought iterate until Brueckner convergence
C is achieved before computing orbital gradients. 

         Call Pccd_drpmo_do_orbrots(Icore(I0),Icrsiz,Iuhf,IGrad_calc,
     +                              IHess_calc,Scale,OOmicroItr,
     +                              OOmacroItr,OOtotalItr)
      Endif

 10   Continue

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

