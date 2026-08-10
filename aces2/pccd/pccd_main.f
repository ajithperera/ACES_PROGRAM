














































































































































































































      SUBROUTINE PCCD(ICORE,ICRSIZ,IUHF)
      Implicit none
      Integer icore(icrsiz),icrsiz,iuhf,i0
      Integer IGrad_calc
      Integer IHess_calc
      Integer Pccd_refr,Itest
      Integer OOmicroItr,OOmacroItr,OOtotalItr
      Integer Length,Idrop
      Double Precision Scale 
      Character*80 Fname
      Character*6 Status
      Logical OOexist 
      DOUBLE PRECISION GLOBNORM
      INTEGER MACRONUM,MICRONUM,TOTALITR
      LOGICAL MACCONVG,MICCONVG
      Logical Opt_orbs, OO_constr, OV_constr, All_mos

c COMMON BLOCKS
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
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
c ----------------------------------------------------------------------
      I0 = 1
      CALL PUTREC(20,"JOBARC","PCCD_RUN",1,1)

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
      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)

      If (Iflags2(138) .Eq. 2) Then
          IGrad_calc = 0
          IHess_calc = 0
          Scale     = Iflags(57)
      Endif

      All_mos   = (Idrop .Eq. 0)
      Opt_orbs  = (Iflags2(178) .Gt. 0)
      Oo_constr = (Iflags2(181) .Eq. 1)
      Ov_constr = (Iflags2(184) .Eq. 1)
C
      If (Igrad_calc .eq. 1 .and. IHess_calc .eq.1) Then

C If Bruekner optimization is sought iterate until Brueckner convergence
C is achieved before computing orbital gradients. 

         If (Iuhf .Eq. 0 .AND. All_mos) Then 
            Call Pccd_banner(6,Status)
            Call Pccd_rdriver(Icore(i0),icrsiz/iintfp,Iuhf)

            Pccd_refr = Iflags2(179)
            Call Getrec(-20, 'JOBARC', 'BRUKTEST', 1, Itest)

            If ((Pccd_refr .Eq. 2) .And. 
     +           Itest .Ne. 1) Go to 10

            Call Pccd_ldriver(Icore(i0),icrsiz/iintfp,Iuhf)

            If (Opt_orbs) Then
               Call Pccd_init_2pdens_lists(Icore(i0),icrsiz/iintfp,Iuhf)
               Call Pccd_init_orbhes_lists(Icore(i0),icrsiz/iintfp,Iuhf)

               Call Pccd_dens(Icore(i0),icrsiz/iintfp,Iuhf)
               Call Pccd_do_orbrots(Icore(I0),icrsiz/iintfp,Iuhf,
     +                              Igrad_calc,IHess_calc,Scale,
     +                              OOmicroItr,OOmacroItr,OOtotalItr)
            Endif 
         Else
            Call Pccd_banner(6,Status)
            Call Pccd_urdriver(Icore(i0),icrsiz/iintfp,Iuhf)

            If ((Pccd_refr .Eq. 2) .And. 
     +           Itest .Ne. 1) Go to 10
            Pccd_refr = Iflags2(179)
            Call Getrec(-20, 'JOBARC', 'BRUKTEST', 1, Itest)

            Call Pccd_uldriver(Icore(i0),icrsiz/iintfp,Iuhf)

            If (Opt_orbs) Then
               Call Pccd_do_orbrots(Icore(I0),icrsiz/iintfp,Iuhf,
     +                              Igrad_calc,IHess_calc,Scale,
     +                              OOmicroItr,OOmacroItr,OOtotalItr)
            Endif 
         Endif 
      Else
            Call Pccd_form_nmrcl_grad_hess(Icore(I0),icrsiz/iintfp,Iuhf,
     +                                  Igrad_calc,IHess_calc,Scale)
      Endif

 10   Continue

C
c ----------------------------------------------------------------------
      RETURN
      End


C ----------------------------------------------------------------------
C Thin PROGRAM wrapper restoring the standalone xPCCD entry point
C (MAIN__), lost when this module was converted to a SUBROUTINE for
C pyaces (see project memory: ACES_PROGRAM unification, 2026-08). This
C calls the shared aces_init/aces_fin pair so the standalone binary
C manages its own memory exactly as the classic driver's own separate
C process used to, before the pyaces conversion removed that from
C inside the subroutine itself.
C ----------------------------------------------------------------------
      PROGRAM XPCCD
      IMPLICIT NONE


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
      INTEGER IUHF
      CALL ACES_INIT(ICORE, I0, ICRSIZ, IUHF, .TRUE.)
      CALL PCCD(ICORE(I0), ICRSIZ, IUHF)
      CALL ACES_FIN
      STOP
      END
