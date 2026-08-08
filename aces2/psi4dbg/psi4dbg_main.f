














































































































































































































      Program  psi4dbg
      Implicit none
      Integer iuhf
      Integer IGrad_calc
      Integer IHess_calc
      Integer Psi4dbg_refr,Itest
      Integer OOmicroItr,OOmacroItr,OOtotalItr
      Integer Length
      Integer Junk,Ncycle
      Double Precision Scale 
      Character*80 Fname
      Logical OOexist 
      Logical Macro_iter
      DOUBLE PRECISION GLOBNORM
      INTEGER MACRONUM,MICRONUM,TOTALITR
      LOGICAL MACCONVG,MICCONVG



c COMMON BLOCKS


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



       
c ----------------------------------------------------------------------
      Call Aces_init(icore, i0,icrsiz, iuhf, .true.)
      If (Iuhf .Eq. 1) Then
         Write(6,*)
         Write(6,"(2a)") " UHF reference pCCD calculations are not",
     +                   " supported changed the refrence to RHF"
         Write(6,"(a)")  " and rerun to proceed."
         Call Errex()
      Endif 
      CALL PUTREC(20,"JOBARC","PCCD_RUN",1,1)

      IGrad_calc = 1
      IHess_calc = 1
      
      Call Gfname("ooinfo.dat",Fname,Length)
      Inquire(file=Fname(1:Length),exist=OOexist)
      If (OOexist) Then
         Open(Unit=1011,File=Fname(1:Length),Status="old")
         Read(1011,*) MACRONUM,MICRONUM,TOTALITR,MACCONVG,
     +                         MICCONVG,GLOBNORM
         call putrec(20,"JOBARC","MACRONUM",1,MACRONUM)
         call putrec(20,"JOBARC","MICRONUM",1,MICRONUM)
         call putrec(20,"JOBARC","TOTALITR",1,TOTALITR)
         CALL putrec(20,"JOBARC","MACCONVG",1,.False.)
         CALL putrec(20,"JOBARC","MICCONVG",1,.False.)
         call putrec(20,"JOBARC","GLOBNORM",1,0.0d0)
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
C
      If (Igrad_calc .eq. 1 .and. IHess_calc .eq.1) Then

C If Bruekner optimization is sought iterate until Brueckner convergence
C is achieved before computing orbital gradients. 

        Call getrec(0,"JOBARC","ORBOPITR",Ncycle,junk)
        If (Ncycle .lt. 0) then
           Macro_iter = .True.
        Else
            call getrec(20,"JOBARC","ORBOPITR",1,ncycle)
            Macro_iter = .True. 
        Endif

         Call Psi4dbg_rdriver(Icore(i0),icrsiz/iintfp,Iuhf,
     +                        .Not.Macro_iter)
         Call Psi4dbg_ldriver(Icore(i0),icrsiz/iintfp,Iuhf,
     +                        .Not.Macro_iter)
         call psi4dbg_check_t2(icore(i0),icrsiz,iuhf)

         Psi4dbg_refr = Iflags2(179)
         Call Getrec(-20, 'JOBARC', 'BRUKTEST', 1, Itest)

         If (Psi4dbg_refr .Eq. 2 .And. Itest .Ne. 1) Go to 10

         Call Psi4dbg_init_2pdens_lists(Icore(i0),icrsiz/iintfp,Iuhf)

         Call Psi4dbg_dens(Icore(i0),icrsiz/iintfp,Iuhf)
         Call Psi4dbg_do_orbrots(Icore(I0),icrsiz/iintfp,Iuhf,
     +                           Igrad_calc,IHess_calc,Scale,
     +                           OOmicroItr,OOmacroItr,
     +                           OOtotalItr)
      Else
         Call Pccd_form_nmrcl_grad_hess(Icore(I0),icrsiz/iintfp,Iuhf,
     +                                  Igrad_calc,IHess_calc,Scale)
      Endif

 10   Continue

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

