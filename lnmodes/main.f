













































































































































































































      Program Main
C
C Description
C
C   Perform Foster-Boys type Localizations of vibrational
C   modes.Ajith Perera 06/2014
C
C Declarations
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Nocupy(2), Itemp(12)
      Character*5 Spntyp(2)
      Character*8 Dump(2), Pscf1(2), Pscf2(2)
      Logical Abort
      Parameter(Lirvpr=600, Oned10 = 1.0D-10, Oned8 = 1.0D-8, 
     &          MXNiter = 150, Mxtry = 3)
C
C Common Block informations


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

      Data Spntyp /'ALPHA', 'BETA'/
      Data Dump /'LOCALMOA', 'LOCALMOB'/
      Data Ione /1/

C All the initilizations 
      call aces_init_rte
      call aces_ja_init
      call getrec(1,'JOBARC','IFLAGS', 100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',500,iflags2)

      icrsiz = iflags(36)
      icore(1) = 0
      do while ((icore(1).eq.0).and.(icrsiz.gt.1000000))
         call aces_malloc(icrsiz,icore,i0)
         if (icore(1).eq.0) icrsiz = icrsiz - 1000000
      end do
      if (icore(1).eq.0) then
         print *, '@MAIN: unable to allocate at least ',
     &            1000000,' integers of memory'
         call aces_exit(1)
      end if
C
      Call Getrec(20, "JOBARC", "NREALATM", Ione, Nrel_atoms)
      Call Getrec(20, "JOBARC", "NATOMS  ", Ione, Ntot_atoms)

      Call get_lnmodes(Icore(i0),icrsiz, Nrel_atoms, Ntot_atoms)

      Call Aces_fin
C
      Stop
      End





