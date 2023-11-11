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
