











      subroutine nl_init(nltitle,err,printdef0)
c
c The nl_xxx series of routines implements a fortran namelist.
c
c A namelist starts with an "*" in the first column followed immediately
c by the name of the namelist.  This is followed by any number of fields
c of the form KEY=VAL separated by spaces or newlines.  Any number of
c fields may appear on a single line provided the line is no more than
c 80 characters long.  The end of the namelist is signaled by the end of
c the file or a line with an "*" in the first column.
c
c Keys do not have to be a rigidly defined word.  Rather, they can be
c abbreviated.  A sample keyword might be "CALC*ULATION".  This means
c that the first four characters are required, and any remaining
c characters may be omitted, but if they are given, must match the
c keyword appropriately.  In other words, "CALC", "CALCUL" and
c "CALCULATION" all work but "CALCCC" won't.


c This opens the ZMAT file and searches for a namelist.  If found, it is
c read into the array namelist for parsing.
c
c If printdef is true, the namelist will default to printing (though it
c can be overridden by the PRINT_NL keyword in the namelist itself).
c Otherwise it will default to not printing.








c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






c This contains the global string for identifying the current subroutine
c or function (provided the programmer set it).  cf. tools/callstack.F
c BE GOOD AND RESET CURR ON EXIT!

      character*64                callstack_curr,callstack_prev
      common /callstack_curr_com/ callstack_curr,callstack_prev
      save   /callstack_curr_com/




c This common block contains information about a namelist which is being
c parsed.

c   nlmaxline   : the maximum number of lines which can be in the namelist
c   nllinelen   : the maximum length of each line

c   nlnumline   : the number of lines in the namelist (blank lines are
c                 removed)
c   nltext      : the text in the namelist
c   prt_nl      : a logical which is read in from the namelist to see
c                 if the values read in are printed or not
c
      integer nlmaxline, nllinelen
      parameter(nlmaxline=64)
      parameter(nllinelen=132)

      character*(nllinelen) nltext(nlmaxline)
      integer nlnumline
      logical prt_nl

      common /namelistc/ nltext
      common /namelist/  nlnumline
      common /namelistl/ prt_nl
      save /namelist/
      save /namelistc/
      save /namelistl/


c This common block contains values for all of the Aces3 keywords as
c specified by default or in the ZMAT file.  For historical reasons,
c the flags are broken up into two blocks (originally only 100 were
c used and later an additional 500 were added).  It is hoped that a
c more flexible way of doing this will be implemented.

c nflags,nflags2  : the number of flags in each block
c iflags,iflags2  : the values for all flags
c iprint          : a special flag containing print information
c iuhf            : 1 if a UHF calculation is done
c needbas         : set to 1 if basis set info is needed
c debug           : a logical flag which may be useful to some
c multipoint      : A flag that says this more than one single point calculation.
c getall          : logical flag that is used in mrcc
c get_nonsense    : if .not. get_nonsense aces3 stops calculation if it thinks it is wrong.
c development_version: indiacates that this is a development version rather than production.

      integer nflags,nflags2
      parameter (nflags=100)
      parameter (nflags2=500)
      
      integer iflags(nflags),iflags2(nflags2),iprint,iuhf,needbas
      logical debug, multipoint, get_nonsense, development_version

      common /flags/  iflags,iflags2,iprint,iuhf,needbas
      common /flagsl/ debug, multipoint, get_nonsense,
     &    development_version

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
c These are parameters containing various file numbers for IO purposes.

c Standard input/output.
      integer stdout,stdin
      parameter (stdin=5)
      parameter (stdout=6)

c MOL file, vmol input created by joda
      integer molio
      character *(*) molfil
      parameter (molio=3)
      parameter (molfil='MOL')

c ZMAT file
      integer zio
      character *(*) zfil
      parameter (zio=4)
      parameter (zfil='ZMAT')

c POLYRATE file
      integer polyio
      character *(*) polyfil
      parameter (polyio=7)
      parameter (polyfil='POLYRATE')

c IIII integral file, created by vmol
      integer iiiiio
      character *(*) iiiifil
      parameter (iiiiio=10)
      parameter (iiiifil='IIII')

c IJIJ integral file, created by vmol
      integer ijijio
      character *(*) ijijfil
      parameter (ijijio=21)
      parameter (ijijfil='IJIJ')

c IIJJ integral file, created by vmol
      integer iijjio
      character *(*) iijjfil
      parameter (iijjio=22)
      parameter (iijjfil='IIJJ')

c IJKL integral file, created by vmol
      integer ijklio
      character *(*) ijklfil
      parameter (ijklio=23)
      parameter (ijklfil='IJKL')

c VPOUT integral file, created by vprops
      integer vpoutio
      character *(*) vpoutfil
      parameter (vpoutio=30)
      parameter (vpoutfil='VPOUT')

c---------------------------------------
c Files for storing lists
c MOINTS file
      integer mointsio
      character *(*) mointsfil0
      parameter (mointsio=50)
      parameter (mointsfil0='MOINTS')

c GAMLAM file
      integer gamlamio
      character *(*) gamlamfil0
      parameter (gamlamio=51)
      parameter (gamlamfil0='GAMLAM')

c MOABCD file
      integer moabcdio
      character *(*) moabcdfil0
      parameter (moabcdio=52)
      parameter (moabcdfil0='MOABCD')

c DERINT file
      integer derintio
      character *(*) derintfil0
      parameter (derintio=53)
      parameter (derintfil0='DERINT')

c DERGAM file
      integer dergamio
      character *(*) dergamfil0
      parameter (dergamio=54)
      parameter (dergamfil0='DERGAM')

c LIST6 file
      integer list6io
      character *(*) list6fil0
      parameter (list6io=55)
      parameter (list6fil0='LIST6')

c LIST7 file
      integer list7io
      character *(*) list7fil0
      parameter (list7io=56)
      parameter (list7fil0='LIST7')

c LIST8 file
      integer list8io
      character *(*) list8fil0
      parameter (list8io=57)
      parameter (list8fil0='LIST8')

c LIST9 file
      integer list9io
      character *(*) list9fil0
      parameter (list9io=58)
      parameter (list9fil0='LIST9')

c LIST10 file
      integer list10io
      character *(*) list10fil0
      parameter (list10io=59)
      parameter (list10fil0='LIST10')

c---------------------------------------

c NEWMOS file
      integer newmosio
      character *(*) newmosfil0
      parameter (newmosio=71)
      parameter (newmosfil0='NEWMOS')

c JOBARC file
      integer jobarcio
      character *(*) jobarcfil0
      parameter (jobarcio=75)
      parameter (jobarcfil0='JOBARC')

c JAINDX file
      integer jaindxio
      character *(*) jaindxfil0
      parameter (jaindxio=75)
      parameter (jaindxfil0='JAINDX')

c FRQARC file
      integer frqio
      character *(*) frqfil
      parameter (frqio=78)
      parameter (frqfil='FRQARC')

c Joda finished file
      integer doneio
      character *(*) donefil
      parameter (doneio=80)
      parameter (donefil='JODADONE')

c NUCDIP file
      integer nucdio
      character *(*) nucdfil
      parameter (nucdio=81)
      parameter (nucdfil='NUCDIP')

c------------------------------------------------------------------
c Module:  intgrt
c------------------------------------------------------------------

c The list of grids to use, one for each radial point in a numerical
c integration.
      integer radgio
      character *(*) radgfil
      parameter (radgio=40)
      parameter (radgfil='RADGRD')

c File for plot points
      integer pltgio
      character *(*) pltgfil
      parameter (pltgio=41)
      parameter (pltgfil='PLTGRD')

c Cartesian coordinates for plot grid points
      integer pgrdptio
      character *(*) pgrdptfil
      parameter (pgrdptio=44)
      parameter (pgrdptfil='PGRDPTS')

c HF electron density
      integer pdensio
      character *(*) pdensfil
      parameter (pdensio=44)
      parameter (pdensfil='PDENSITY')

c HF kinetic energy
      integer pekinio
      character *(*) pekinfil
      parameter (pekinio=44)
      parameter (pekinfil='PKINETIC')

c HF nuclear-electron attraction energy
      integer penatrio
      character *(*) penatrfil
      parameter (penatrio=44)
      parameter (penatrfil='PNATR')

c LDA exchange energy
      integer pdldaxio
      character *(*) pdldaxfil
      parameter (pdldaxio=44)
      parameter (pdldaxfil='PXLDAX')

c Becke exchange energy
      integer pdbeckeio
      character *(*) pdbeckefil
      parameter (pdbeckeio=44)
      parameter (pdbeckefil='PBECKE')

c LDA energy
      integer pdldaio
      character *(*) pdldafil
      parameter (pdldaio=44)
      parameter (pdldafil='PLDA')

c LYP correlation energy
      integer pdlypio
      character *(*) pdlypfil
      parameter (pdlypio=44)
      parameter (pdlypfil='PLYP')

c Thomas Fermi kinetic energy 
      integer pdtfkeio
      character *(*) pdtfkefil
      parameter (pdtfkeio=44)
      parameter (pdtfkefil='PTFKE')

c Weizsacker kinetic energy 
      integer pdwkeio
      character *(*) pdwkefil
      parameter (pdwkeio=44)
      parameter (pdwkefil='PWKE')

c Projected Fock potential
      integer pdfockio
      character *(*) pdfockfil
      parameter (pdfockio=44)
      parameter (pdfockfil='PFOCKP')

c Projected Brueckner potential
      integer pdbrknrio
      character *(*) pdbrknrfil
      parameter (pdbrknrio=44)
      parameter (pdbrknrfil='PBRKNRP')

c Coulomb energy 
      integer pcoulio
      character *(*) pcoulfil
      parameter (pcoulio=44)
      parameter (pcoulfil='PCOUL')

c Coulomb potential
      integer pcpotio
      character *(*) pcpotfil
      parameter (pcpotio=44)
      parameter (pcpotfil='PCPOT')

c Exchange energy 
      integer pexchio
      character *(*) pexchfil
      parameter (pexchio=44)
      parameter (pexchfil='PEXCH')

c Exchange potential
      integer pepotio
      character *(*) pepotfil
      parameter (pepotio=44)
      parameter (pepotfil='PEPOT')

c Exchange DFT potential
      integer pexdftio
      character *(*) pexdftfil
      parameter (pexdftio=44)
      parameter (pexdftfil='PEXDFT')

c Correlated electron density
      integer cdensio
      character *(*) cdensfil
      parameter (cdensio=45)
      parameter (cdensfil='CDENSITY')

c Correlated kinetic energy
      integer cekinio
      character *(*) cekinfil
      parameter (cekinio=45)
      parameter (cekinfil='CKINETIC')

c Correlated nuclear-electron attraction energy
      integer cenatrio
      character *(*) cenatrfil
      parameter (cenatrio=45)
      parameter (cenatrfil='CNATR')

c Thomas Fermi kinetic energy (correlated density)
      integer cdtfkeio
      character *(*) cdtfkefil
      parameter (cdtfkeio=45)
      parameter (cdtfkefil='CTFKE')

c Weizsacker kinetic energy (correlated density)
      integer cdwkeio
      character *(*) cdwkefil
      parameter (cdwkeio=45)
      parameter (cdwkefil='CWKE')

c Coulomb energy (correlated density)
      integer ccoulio
      character *(*) ccoulfil
      parameter (ccoulio=45)
      parameter (ccoulfil='CCOUL')

c Coulomb potential(correlated density)
      integer ccpotio
      character *(*) ccpotfil
      parameter (ccpotio=45)
      parameter (ccpotfil='CCPOT')

c Exchange energy (correlated density)
      integer cexchio
      character *(*) cexchfil
      parameter (cexchio=45)
      parameter (cexchfil='CEXCH')

c Exchange potential(correlated density)
      integer cepotio
      character *(*) cepotfil
      parameter (cepotio=45)
      parameter (cepotfil='CEPOT')

c Exchange DFT potential
      integer cexdftio
      character *(*) cexdftfil
      parameter (cexdftio=45)
      parameter (cexdftfil='CEXDFT')

c LDA exchange energy (correlated density)
      integer cdldaxio
      character *(*) cdldaxfil
      parameter (cdldaxio=45)
      parameter (cdldaxfil='CXLDAX')

c Becke exchange energy (correlated density)
      integer cdbeckeio
      character *(*) cdbeckefil
      parameter (cdbeckeio=45)
      parameter (cdbeckefil='CBECKE')

c LDA energy (correlated density)
      integer cdldaio
      character *(*) cdldafil
      parameter (cdldaio=45)
      parameter (cdldafil='CLDA')

c LYP correlation energy (correlated density)
      integer cdlypio
      character *(*) cdlypfil
      parameter (cdlypio=45)
      parameter (cdlypfil='CLYP')

c Coupled-cluster energy
      integer eccradio
      character *(*) eccradfil
      parameter (eccradio=45)
      parameter (eccradfil='ECCRAD')

      character *(80) mointsfil, gamlamfil, moabcdfil,
     &    derintfil, dergamfil, list6fil, list7fil,
     &    list8fil, list9fil, list10fil,
     &    newmosfil, jobarcfil, jaindxfil
      common /filenames/  mointsfil, gamlamfil, moabcdfil,
     &    derintfil, dergamfil, list6fil, list7fil,
     &    list8fil, list9fil, list10fil,
     &    newmosfil, jobarcfil, jaindxfil

c------------------------------------------------------------------

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
      character *(*) nltitle
      integer err
      logical printdef, printdef0

      character*(nllinelen) line,title
      character*(nllinelen) readline
      integer linenum,l,strlen
      logical eof,streq

      callstack_curr='NL_INIT'
c
      if (.not. multipoint .or. iprint .ge. 1) then 
        printdef = printdef0
      else
        printdef = .false.
      endif
c
c Open ZMAT file and search for the namelist
      open(unit=zio,file=zfil,form='formatted')
      linenum=0
      nlnumline=0
      title='*'//nltitle
      l=strlen(title)
      eof=.false.
   10 line=readline(linenum,zio,eof)
      if (eof) goto 999
      if (.not.streq(title,line(1:l),.true.)) goto 10

c Once the namelist is found, remove the '*TITLE' from it and if
c anything is left, put it as the first line of the namelist text
      line=line(l+1:strlen(line))
      if (strlen(line).gt.0) then
        call upcase(line)
        nltext(1)=line
        nlnumline=1
      endif

c Read until we are done with the file or find a line starting with '*'
   20 line=readline(linenum,zio,eof)
      if (eof) goto 999
      if (line(1:1).eq.'*') goto 999
      nlnumline=nlnumline+1
      if (nlnumline.gt.nlmaxline) then
        write(stdout,'(a)')
     &      '@NL_INIT-F:  maximum length of namelist exceeded'
        call errex
      end if
      call upcase(line)
      nltext(nlnumline)=line
      goto 20

  999 close(zio)
      call nl_log('print_nl',printdef,prt_nl)
      if (prt_nl) call nl_prttop(nltitle)
      if (nlnumline.eq.0) err=1
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
