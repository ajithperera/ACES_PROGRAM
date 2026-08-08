













































































































































































































      Subroutine Run1props(Iuhf)

      Implicit None
      Integer Iuhf 
      Logical FileExist

c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end
c istart.com : begin
      integer*8       i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end
c flags.com : begin
      integer*8      iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer*8       iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer*8       iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end
c jobarc.com : begin

c This data tracks the contents of the JOBARC file. 'physical' records refer
c to direct I/O while 'logical' records refer to the ACES archive elements.

c marker(i) : the name of logical record i
c rloc(i)   : the integer index in JOBARC that starts logical record i
c rsize(i)  : the integer-length of logical record i
c nrecs  : the number of physical records in the JOBARC file
c irecwd : the integer-length of a physical record
c irecln : the    recl-length of a physical record

      integer*8 MAX_JA_RECS 
      integer*8 JA_RECWD
      integer*8 JA_UNIT 
      
      parameter(MAX_JA_RECS = 1000)
      parameter(JA_RECWD    = 128)
      parameter(JA_UNIT     = 75)

      character*8     marker(MAX_JA_RECS)
      integer*8       rloc  (MAX_JA_RECS),
     &                rsize (MAX_JA_RECS),
     &                nrecs, irecwd, irecln
      common /jobarc/ marker,
     &                rloc,
     &                rsize,
     &                nrecs, irecwd, irecln
      save   /jobarc/

c bJAUp  : a flag for bombing in get/putrec if aces_ja_init has not been called
c bJAMod : a flag for updating JAINDX in aces_ja_fin

      logical           bJAUp, bJAMod
      common /ja_flags/ bJAUp, bJAMod
      save   /ja_flags/

c jobarc.com : end

      FileExist=.false.
      inquire(file='IIII',exist=FileExist)
 
      If (FileExist) Then 
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
         Call Propints(icore(i0),icrsiz/iintfp,icore(i0),icrsiz)
      Else 
         Call Ints_(Icore(i0),Icrsiz,Iuhf)
         Call V2ja(icore(i0),Icrsiz,Iuhf)
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
         Call Propints(icore(i0),icrsiz/iintfp,icore(i0),icrsiz)
      Endif 

      If (Iflags(2) .Gt. 0) Then
         Call aces_init_chemsys()
         Call Vtran(icore(i0),Icrsiz,Iuhf)
         Call aces_init_chemsys()
         Call Init_legacies()
         Call Intproc(icore(i0),Icrsiz,Iuhf)
         Call Vcc(icore(i0),Icrsiz,Iuhf)
         Call Vlambda(icore(i0),Icrsiz,Iuhf)
         Call Vdens(icore(i0),Icrsiz,Iuhf)
         Call Cmp_props(Icore(I0),Icrsiz,Iuhf)
      Else
         Call Cmp_props(Icore(I0),Icrsiz,Iuhf)
      Endif 

      Call aces_fin
      Call C_free(Icore)

      Return
      End

