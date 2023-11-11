c lists.com : begin
c These common blocks contain global information about the arrays in storage.
c Elements prepended with "bw" are for storing the file metadata while working
c on multiple references.

c moio  (iGrp,iFam) : the physical record that contains the first element
c                     of the array (iGrp,iFam)
c moiowd(iGrp,iFam) : the integer-word index of the first element
c moiods(iGrp,iFam) : the number of columns in the array
c moiosz(iGrp,iFam) : the number of rows    in the array
c moiofl(iGrp,iFam) : the external file unit that contains the array

      integer*8 MAX_IO_GRPS
      integer*8 MAX_IO_FAMS
      integer*8 FAMS_PER_FILE
      integer*8 FIRST_IO_LUN
      integer*8 MAX_IO_LUNS
      integer*8 maxref 

      parameter(MAX_IO_GRPS=10)
      parameter(MAX_IO_FAMS=500)
      parameter(FAMS_PER_FILE=100)
      parameter(FIRST_IO_LUN=50)
      parameter(MAX_IO_LUNS=5)
      parameter(MAXREF=1)

      integer*8        moio  (MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiowd(MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiosz(MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiods(MAX_IO_GRPS,MAX_IO_FAMS),
     &                 moiofl(MAX_IO_GRPS,MAX_IO_FAMS),
     &               bwmoio  (MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiowd(MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiosz(MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiods(MAX_IO_GRPS,MAX_IO_FAMS,maxref),
     &               bwmoiofl(MAX_IO_GRPS,MAX_IO_FAMS,maxref)
      common /lists/   moio,   moiowd,   moiosz,   moiods,   moiofl,
     &               bwmoio, bwmoiowd, bwmoiosz, bwmoiods, bwmoiofl
      save   /lists/

c moiomxsz(iGrp,iFam) : the original length of a one-dimensional array
c                       (This is shameful. Arrays should not be re-dimensioned
c                        at will during a job.)

      integer*8           moiomxsz(MAX_IO_GRPS,MAX_IO_FAMS),
     &                    bwmoiomxsz(MAX_IO_GRPS,MAX_IO_FAMS,maxref)
      common /lists_mxsz/   moiomxsz,
     &                    bwmoiomxsz
      save   /lists_mxsz/

c pRec(i)    : the index of the physical record in file i containing free space
c              (i is the internal unit number of the storage file.)
c iIntOff(i) : the integer offset from the beginning of the physical record
c              needed to address the free space

      integer*8          pRec   (MAX_IO_LUNS),
     &                   iIntOff(MAX_IO_LUNS),
     &                 bwpRec   (MAX_IO_LUNS,maxref),
     &                 bwiIntOff(MAX_IO_LUNS,maxref)
      common /io_ptrs/   pRec,   iIntOff,
     &                 bwpRec, bwiIntOff
      save   /io_ptrs/

c bIOUp  : a flag for bombing in get/putlst if aces_io_init has not been called
c bIOMod : a flag for updating the records in aces_io_fin

      logical           bIOUp, bIOMod
      common /io_flags/ bIOUp, bIOMod
      save   /io_flags/

c lists.com : end
