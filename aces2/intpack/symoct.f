










      subroutine symoct
     &    (atomchrg,atmvc,rij,aij,cdnt,rsqrd,rrtmp,wtintr,
     &    totwt,radgrid,radpt,rwt,integaxis,gridxyz,gridwt,grdangpts,
     &    iradpt,iangpt,grid,evalpt)

c This routine sets up the cartesian coordinates for the numerical 
c integration grid








c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>






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
     &    newmosfil, jobarcfil, jaindxfil
      common /filenames/  mointsfil, gamlamfil, moabcdfil,
     &    derintfil, dergamfil, list6fil, list7fil,
     &    newmosfil, jobarcfil, jaindxfil

c------------------------------------------------------------------




c This contains flags that are set in the INTGRT namelist.  See the
c file initintgrt.F for a description of each of them.
c
c The following are exceptions:
c    int_ks           : .true. if we are doing Kohn-Sham
c    int_ks_finaliter : .true. if this is the final iteration of KS
c    int_ks_exch      : which potential to use to calculate exchange
c    int_ks_corr      : which potential to use to calculate correlation
c    int_kspot        : which hybrid functional to use to calculate potential
c                       (if equal to fun_special, use int_ks_exch and
c                       int_ks_corr)
c    int_dft_fun      : which functional to use with any SCF density
c                       Added and modified by Stan Ivanov
c                       (if fun_special the functional is user-defined
c                        if fun_hyb_name then use hybrid functional) 
c    int_printlev     : 0 if we are doing a dft calculation, 1 if we
c                       are doing the final iteration of a KS calculation,
c                       2 if we are doing a KS iteration.
c These are set in the calling routines, NOT in the namelist.
c                     : Additions by S. Ivanov
c     num_acc_ks      : .true.  if numerical accelerator is used for KS 
c                        Default is .true.
c     ks_exact_ex     : .true. if exact LOCAL exchange is used for KS
c                        Deafult is .false.
c     int_tdks        : .true. if time-dependent KS calculation is
c                        requested
c                        Default is .false.
c     int_ks_scf      : .true. if the actual KS SCF energy is being
c                        calculated and printed out. Default is .false.
c
      integer int_numradpts,int_radtyp,int_partpoly,int_radscal,
     &    int_parttyp,int_fuzzyiter,int_defenegrid,int_defenetype,
     &    int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp

      logical int_ks,num_acc_ks,ks_exact_ex,int_tdks,int_ks_scf,
     &        int_ks_finaliter

      double precision
     &    int_radlimit,coef_pot_nonlocal

      common /intgrtflags/  int_numradpts,int_radtyp,int_partpoly,
     &    int_radscal,int_parttyp,int_fuzzyiter,int_defenegrid,
     &    int_defenetype,int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_ks_finaliter,int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem,int_overlp
c

c  prakash added int_ksmem to the common block
      common /intgrtflagsd/ int_radlimit,coef_pot_nonlocal
      common /intgrtflagsl/ int_ks,num_acc_ks,ks_exact_ex,int_tdks,
     &                      int_ks_scf

      save /intgrtflags/
      save /intgrtflagsl/
      save /intgrtflagsd/

c The following are parameters used in the namelist

      integer int_prt_never,int_prt_dft,int_prt_ks,int_prt_always
      parameter (int_prt_never   =1)
      parameter (int_prt_dft     =2)
      parameter (int_prt_ks      =3)
      parameter (int_prt_always  =4)

      integer int_radtyp_handy,int_radtyp_gl
      parameter (int_radtyp_handy=1)
      parameter (int_radtyp_gl   =2)

      integer int_partpoly_equal,int_partpoly_bsrad,
     &    int_partpoly_dynamic
      parameter (int_partpoly_equal  =1)
      parameter (int_partpoly_bsrad  =2)
      parameter (int_partpoly_dynamic=3)

      integer int_radscal_none,int_radscal_slater
      parameter (int_radscal_none  =1)
      parameter (int_radscal_slater=2)

      integer int_parttyp_rigid,int_parttyp_fuzzy
      parameter (int_parttyp_rigid=1)
      parameter (int_parttyp_fuzzy=2)

      integer int_gridtype_leb
      parameter (int_gridtype_leb=1)






c This file contains a number of physical constants and conversion
c factors.  Physical constants are all stored in a variable named
c CONST_name where 'name' is the name of the physical constant.
c 
c Conversion factors are stored in variables named CONV_unit1_unit2
c where 'unit1' and 'unit2' are the names of two physical units.  To
c get from unit1 to unit2, multiply by this factor.
c 
c Example: To convert 20.4 meters/s to miles/hour in a perl script,
c          include the following line at the top of the script:
c              require "Constants.pl";
c          and then, at the point in the script where the conversion
c          is required, the following returns the desired value:
c              20.4 * $CONV_m_mile / $CONV_s_hr
c 
c          To do the same thing in either a C or Fortran program, include
c          the file "Constants.h" or "Constants.f" as appropriate and
c          at the point where the conversion is required, the following
c          returns the desired value:
c              20.4 * CONV_m_mile / CONV_s_hr
c 
c This file was generated automatically on 6/4/96.
c 
c Do not edit this file.  If you wish to add or change a conversion
c factor or constant, edit the file GenerateConstants.pl and rerun it.





c Physical Constants
c ====================
c 
c e       = e
c pi      = pi
c c       = speed of light in a vacuum     m s-1
c g       = gravitational acceleration     m s-2
c G       = gravitational constant         N m+2 kg-2
c me      = mass of an electron            kg
c mn      = mass of a neutron              kg
c mp      = mass of a proton               kg
c mmu     = mass of a mu particle          kg
c u       = atomic mass unit               kg
c ec      = elementary charge              C
c h       = planck's' constant              J s
c hbar    = h/2 pi                         J s
c k       = boltmann's' constant            J K-1
c u0      = permeability of vacuum         N A-2
c e0      = permittivity of vacuum         C+2 N-1 m-2
c re      = classical electron radius      m
c alpha   = fine structure constant   
c a0      = bohr radius                    m
c RH      = quantum hole resistance        ohm
c Rh      = Rydberg constant               m-1
c phi0    = magnetic flux quantum (h/2 ec) m+2 kg s-2 A-1
c uB      = Bohr magneton                  m+2 A
c ue      = electron magnetic moment       m+2 A
c un      = neutron magnetic moment        m+2 A
c uN      = nuclear magneton               m+2 A
c up      = proton magnetic moment         m+2 A
c umu     = mu particle magnetic moment    m+2 A
c lambdac = compton electron wavelength    m
c lambdacp= compton proton wavelength      m
c sigma   = Stefan-Boltmann constant       W m-2 K-4
c NA      = avogadro's' number              mole-1
c Vm      = ideal gas volume at STP        m+3 mole-1
c R       = gas constant                   J K-1 mole-1
c F       = faradays constant              C mole-1
c ea      = atomic unit of energy          J








c TIME (s)
c ====================
c 
c s      = second
c min    = minute
c hr     = hour
c day    = day
c week   = week
c yr     = calendar year (365 days)
c yrleap = calendar year (leap year)
c 
c yrave  = average year (calendar year averaged over 4 years)
c yrside = year (sidereal)
c yrtrop = year (tropical)
c monave = month (averaged over 4 calendar years)
c daysid = sidereal day
c 
c shake  = shake






c LENGTH (m)
c ====================
c 
c METRIC
c km     = kilometer
c m      = meter (SI)
c cm     = cm
c mm     = mm
c um     = micrometer
c nm     = nanometer
c pm     = picometer
c micron = micron
c mmicr  = millimicron
c fermi  = fermi
c 
c ATOMIC
c a      = angstrom
c a0     = bohr radius
c 
c AMERICAN/BRITISH
c hand   = hand
c ell    = ell
c in     = inch
c ft     = foot
c yd     = yard
c mile   = mile
c mileu  = mile (US survey)
c mil    = mil
c rod    = rod
c fur    = furlong
c chaing = chain (Gunter's')
c chainr = chain (Ramsden's')
c leag   = league
c cable  = cable length (U.S.)
c calib  = caliber
c cubit  = cubit
c ftu    = foot (US survey)
c barley = barleycorn (Brit)
c x      = x-unit
c span   = span
c nail   = nail (Brit)
c 
c NAUTICAL
c milen  = mile (nautical)
c leagn  = league (nautical)
c fathom = fathom
c degn   = nautical degree
c circn  = nautical circle
c cablen = cable length (international)
c 
c ASTRONOMICAL
c ly     = light year
c au     = astronomical unit
c pc     = parsec
c ls     = light second
c lm     = light minute
c 
c MISC
c bolt   = bolt (of cloth)






c MASS (kg)
c ====================
c 
c METRIC
c g      = gram
c kg     = kilogram (SI unit of mass)
c ktonm  = kiloton (metric)
c tonne  = tonne
c tonm   = ton (metric)
c 100wtm = hundredweight (metric)
c caratm = carat (metric)
c ng     = nanogram
c pg     = picogram
c mg     = milligram
c ug     = microgram
c 
c ATOMIC
c amu    = atomic mass units
c 
c AVOIRDUPOIS (US)
c 100wt  = hundredweight (short)
c 100wtl = hundredweight (long)
c cental = cental
c dram   = dram (solid)
c geelb  = geepound
c lb     = pound
c oz     = ounce (avoirdupois)
c slug   = slug
c stone  = stone
c ton    = ton (short)
c tonl   = ton (long)
c 
c TROY
c ozt    = ounce (Troy or apothecary)
c penny  = pennyweight
c scruple= scruple
c grain  = grain
c lbt    = pound (Troy)
c dramt  = dram (Troy or apothecary)






c ELECTRIC CURRENT (A)
c ====================
c 
c abA    = abampere
c A      = ampere
c Aint   = ampere (international)
c Aus    = ampere (U.S.)
c biot   = biot
c gilb   = gilbert






c TEMPERATURE INTERVAL (K)
c ====================
c 
c degc   = celcius degree
c degf   = farenheit degree
c degr   = rankine degree
c K      = kelvin degree






c LUMINOUS INTENSITY
c ====================
c 
c cd     = candela
c hef    = hefner unit
c lumPsr = lumen per steradian






c ANGLES (rad)
c ====================
c 
c amin   = minutes of an angle
c as     = seconds of an angle
c circum = circumference
c deg    = degrees
c gon    = gon (grade)
c quad   = quadrant
c rad    = radians
c rev    = revolution






c SOLID ANGLES (sr)
c ====================
c 
c sphere = sphere
c sr     = steradians
c sqdeg  = square degree
c sphra  = spherical right angle






c AREA (m2)
c ====================
c 
c METRIC
c cm2    = square cm
c km2    = square kilometer
c m2     = square meter
c are    = are
c circmm = circular millimeter
c hect   = hectare
c mm2    = square millimeter
c barn   = barn
c 
c US
c acre   = acre
c ft2    = square ft
c in2    = square inch
c mile2  = square mile
c yd2    = square yd
c acreus = acre (US survey)
c chaig2 = square chain (Gunter's')
c chair2 = square chain (Ramsden's')
c chaiu2 = square chain (US survey)
c circin = circular inch
c circmil= circular mil
c ftu2   = square foot (US survey)
c linkg2 = square link (Gunter's')
c linkr2 = square link (Ramsden's')
c mil2   = square mil
c mileu2 = square mile (US survey)
c rod2   = square rod
c town   = township (US)
c darcy  = darcy






c VOLUME (l)
c ====================
c 
c METRIC
c cc     = cubic centimeter
c cl     = centiliter
c cm3    = cubic centimeter
c cupm   = cup (metric)
c dm3    = cubic decimeter
c km3    = cubic kilometer
c l      = liter
c litero = old (1901-1964) value of liter
c m3     = cubic meter
c ml     = milliliter
c mm3    = cubic millimeter
c 
c US/British
c acreft = acre-foot
c acrein = acre-inch
c bag    = bag (Brit)
c bbl    = barrel (petroleum)
c bblbb  = barrel (Brit, beer)
c bblbw  = barrel (Brit, wine)
c bbluc  = barrel (US, cranb)
c bblud  = barrel (US, dry)
c bblul  = barrel (US, liquid)
c board  = board foot
c boardf = board foot
c bu     = bushel (US)
c bub    = bushel (Brit)
c buck   = bucket (Brit)
c cord   = cord
c cordft = cord-foot
c cup    = cup
c dra    = dram (US, liquid)
c drach  = drachm (Brit, liquid)
c firkb  = firkin (Brit)
c firku  = firkin (US)
c floz   = ounce (US, liquid)
c ft3    = cubic foot
c gal    = gallon (US, liquid)
c galb   = gallon (Brit)
c gald   = gallon (US, dry)
c gill   = gill (US)
c gillb  = gill (Brit)
c in3    = cubic inch
c mile3  = cubic mile
c minim  = minim (US)
c minimb = minim (Brit)
c ozbf   = ounce (Brit, liquid)
c peck   = peck (US)
c peckb  = peck (Brit)
c pt     = pint (US, liquid)
c ptb    = pint (Brit)
c ptd    = pint (US, dry)
c qt     = quart (US, liquid)
c qtb    = quart (Brit)
c qtd    = quart (US, dry)
c regton = register ton
c scrup  = scruple (Brit, liquid)
c seam   = seam (Brit)
c yd3    = cubic yard






c FREQUENCY (Hz)
c ====================
c 
c Hz     = hertz
c pS     = per seconds (s-1)






c FORCE (N)
c ====================
c 
c dyn    = dyne
c kgf    = kg-force
c N      = newton
c lbal   = poundal
c lbf    = pound-force
c mgf    = mg-force






c PRESSURE (Pa)
c ====================
c 
c atm      = atmosphere
c atmt     = atmosphere (tech)
c bar      = bar
c barye    = barye
c cmh2o    = cm of water
c cmhg     = cm of mercury
c dynPcm2  = dyne per square cm
c fth2o    = foot of water
c gPcm2    = gram-force per square cm
c inh2o    = inch of water
c inhg     = inch of mercury
c kgPcm2   = kg-force per square cm
c kgPm2    = kg-force per square m
c kgPmm2   = kg-force per square mm
c kPa      = kilopascal
c lbalPft2 = poundal per square foot
c lbPft2   = pound-force per square foot
c lbPin2   = pound-force per square inch
c mbar     = millibar
c megaPa   = megapascal
c mh2o     = meter of water
c mmh2o    = mm of water
c mmhg     = mm of mercury
c NPcm2    = newton per square cm
c NPm2     = newton per square m
c NPmm2    = newton per square mm
c Pa       = pascal (N m-2) (SI)
c psf      = pound-force per square foot
c psi      = pound-force per square inche
c t        = torr
c tonlPft2 = ton-force (long) per square foot
c tonlPin2 = ton-force (long) per square inch
c tonmPin2 = ton-force (metric) per square inch
c tonmPm2  = ton-force (metric) per square meter
c tonsPft2 = ton-force (short) per square foot
c tonsPin2 = ton-force (short) per square inch






c ENERGY, WORK, QUANTITY OF HEAT (J)
c ====================
c 
c cal      = calorie (thermochemical)
c kcal     = kilocalorie (thermochemical)
c ea       = hartree
c erg      = erg (g cm+2 s-2)
c ev       = electron volt
c hz       = energy measured as hertz
c J        = joule (N m) (SI)
c Nm       = newton-meter
c megaJ    = megajoule
c kcalPmol = kcal/mole
c kJPmol   = kjoule/mole
c wn       = cm-1
c ryd      = rydberg
c btu      = British thermal unit
c ft3atm   = cubic foot-atmosphere
c ft3lbPin2= cubic foot-pount-force/square inche
c ftlbal   = foot-poundal
c ftlb     = foot-pound-force
c hphr     = horsepower-hour
c hphrm    = horsepower-hour (metric)
c kgm      = kilogram-force-meter
c gcm      = gram-force-cm
c kWhr     = kilowatt-hour
c gWhr     = gigawatt-hour
c latm     = liter-atmosphere
c lbar     = liter-bar
c Whr      = watt-hour
c Ws       = watt-sec
c dyncm    = dyne-cm
c chu      = centigrade heat unit
c btu39    = btu (39 deg F, 4 deg C)
c but60    = btu (60 deg F, 15.6 deg C)
c btuave   = btu (mean)
c btuth    = btu (thermochemical)
c calor    = calorie
c kcalor   = kilocalorie
c cal15    = calorie (15 deg C)
c cal20    = calorie (20 deg C)
c calave   = calorie (mean)
c calth    = calorie (thermochemical)
c ccatm    = cubic centimeter-atmosphere






c POWER, RADIANT FLUX (W)
c ====================
c 
c btuPhr  = btu per hour
c btuPmin = btu per minute
c btuPs   = btu per sec
c calPmin = calorie per minute
c calPs   = calorie per second
c ergPs   = erg per second
c ftlbPhr = foot-pound-force per hour
c ftlbPmin= foot-pound-force per minute
c ftlbPs  = foot-pound-force per second
c hp      = horsepower
c hpb     = horsepower (boiler)
c hpe     = horsepower (electric)
c hpm     = horsepower (metric)
c hpw     = horsepower (water)
c JPhr    = joule per hour
c JPmin   = joule per minute
c JPs     = joule per second
c kcalPhr = kilocalories per hour
c kcalPmin= kilocalories per minute
c kcalPs  = kilocalories per second
c kgmPhr  = kilogram-force-meter per hr
c kgmPmin = kilogram-force-meter per minute
c kgmPs   = kilogram-force-meter per second
c W       = watt
c kW      = kilowatt
c Wave    = watt (int. mean)
c Wu      = watt (int. US)






c QUANTITY OF ELECTRICITY, ELECTRIC CHARGE (C)
c ====================
c 
c abC    = abcoulomb
c As     = ampere-second
c Ahr    = ampere-hour
c C      = coulomb
c frank  = franklin






c ELECTRIC POTENTIAL, POTENTIAL DIFFERENCE, ELECTROMOTIVE FORCE (V)
c ====================
c 
c abV    = abvolt
c kV     = kilovolt
c V      = volt
c Vave   = volt (int. mean)
c Vus    = volt (int. US)






c CAPACITANCE (F)
c ====================
c 
c abF    = abfarad
c F      = farad
c Fave   = farad (int. mean)
c Fus    = farad (int. US)
c uF     = microfarad






c ELECTRICAL RESISTANCE
c ====================
c 
c abO    = abohm
c O      = ohm
c Oave   = ohm (int. mean)
c Ous    = ohm (int. US)






c CAPACITANCE (S)
c ====================
c 
c S      = siemens
c mhO    = mho (omh-1)
c abmhO  = abmho






c MAGNETIC FLUX (Wb)
c ====================
c 
c max   = maxwell
c pole  = unit pole
c Vs    = volt-sond
c Wb    = weber






c MAGNETIC FLUX DENSITY (T)
c ====================
c 
c T        = tesla
c WbPm2    = weber per square meter






c INDUCTANCE (H)
c ====================
c 
c abH    = abhenry
c H      = henry
c Have   = henry (int. mean)
c Hus    = henry (int. US)






c ILLUMINANCE (lx)
c ====================
c 
c lx     = lux
c lmPcm2 = lumen per square cm
c lmPm2  = lumen per square meter
c lmPft2 = lumen per square foot
c ph     = phot






c ACTIVITY (Bq)
c ====================
c 
c Bq     = becquerel
c cur    = curie






c ABSORBED DOSE (Gy)
c ====================
c 
c Gy     = gray
c JPkg   = joule per kg






c SOUND INTENSITY
c ====================
c 
c db     = decibel
c neper  = neper






      double precision
     &    zero,one,two,three,four,five,six,seven,eight,nine,ten

      parameter (zero =0.0d0)
      parameter (one  =1.0d0)
      parameter (two  =2.0d0)
      parameter (three=3.0d0)
      parameter (four =4.0d0)
      parameter (five =5.0d0)
      parameter (six  =6.0d0)
      parameter (seven=7.0d0)
      parameter (eight=8.0d0)
      parameter (nine =9.0d0)
      parameter (ten  =1.0d1)

      double precision
     &    half,third,fourth,fifth,sixth

      parameter (half  =one/two)
      parameter (third =one/three)
      parameter (fourth=one/four)
      parameter (fifth =one/five)
      parameter (sixth =one/six)

      double precision
     &    pi

      parameter (pi    =3.14159265358979d0)

c Hopefully we can get away from storing a bunch of strings

c Abelian groups

c The following constants are generated in a3const.F.

      double precision
     &    CONST_nan

      common /const/ CONST_nan




c This common block contains the molecule and basis set information
c read in from the JOBARC file.  Since much of this information is
c used in a large number of modules, and since most of the information
c is relatively small compared to the other things held in memory,
c a large percentage of the data stored in the JOBARC file is stored
c here, even though some modules will not use all of it.

c   maxangshell - The maximum number of angular momentum shells.  Since this
c                 is used VERY infrequently, set it high enough to never
c                 cause a problem.
c   spinc(2)    - The characters 'A' and 'B' (useful for alpha/beta labels)
c   natoms      - The number of atoms in the Z-matrix (including X/GH).  After
c                 remove is called, natoms becomes equivalent to nrealatm.
c   natomsx     - The number of atoms in the Z-matrix (including X/GH).  This
c                 does not change.
c   nrealatm    - The number of atoms in the Z-matrix (including GH).
c   naobasfn    - The number of AOs in the basis
c   nbastot     - The number of symmetry adapted orbitals in the basis (the AO
c                 basis may be larger than the SO basis if spherical orbitals
c                 are used since harmonic contaminants are deleted)
c   linear      - 1 if the molecule is linear
c   orientmt    - 3x3 matrix which relates the computational and canonical
c                 orientations
c   nucrep      - Nuclear repulsion energy in a.u.
c   nmproton    - Number of protons in the molecule.
c
c   compptgp    - Point group
c   fullptgp    -
c   compordr    - Order of the point group
c   fullordr    -
c   compnirr    - Number of irreps in the point group
c   fullnirr    -
c   compnorb    - Number of unique atoms (orbits) in the point group
c   fullnorb    -
c   c1symmet    - 1 if the molecule is C1 symmetry
c   nirrep      - The same as compnirr (since nirrep is used so commonly,
c                 this is included for conveniance)
c                 ***NOTE*** nirrep is read in twice and is stored in /sym/
c                            so it is not actually included here
c   totprim     - Total number of primitive functions in the molecule
c   maxshlprim  - Largest number of primitives in a single shell
c   maxshlao    - Largest number of AOs in a single shell
c   maxshlorb   - Largest number of primitive orbitals (primitive functions
c                 times the number of AOs) in a single shell
c   maxangmom   - Largest angular momentum for any atom
c   maxshell    - Larges number of angular momentum shells for any atom
c   noccorb(2)  - The number of alpha and beta occupied orbitals
c   nvrtorb(2)  - The number of alpha and beta virtual orbitals

c The parameter maxorbit is needed because of how dynamic memory is used.
c Two runs of the program are needed.  The first to calculate memory usage,
c the second to use it.  In order to calculate totprim, we have to know the
c orbit population vector (the number of each type of atom).  BUT, this is
c stored in dynamic memory since we do not know how long this vector is.
c In the future, joda or vmol will write this information to JOBARC, and
c this problem will disappear.  In the meantime, we have to introduce a
c genuine limit on the size of the molecule.  It may have no more than
c maxorbit sets of unique atoms.  This limit is ONLY used in the subroutine
c basis, so it probably will disappear when the information in the MOL file
c is put in JOBARC.
c    maxorbit   - the number of symmetry unique atoms

c The following are pointers to real arrays
c
c   zatommass(natoms)  - Atomic mass of all atoms (X=0.0, GH=100.0)
c   zcoord(3,natoms)   - Coordinates of all atoms (computational orientation)
c   zalpha(totprim)    - The alpha for each primitive function
c   zprimcoef(totprim,naobasfn)
c                      - The primitive to AO coefficients
c
c The following are pointers to integer arrays
c
c   patomchrg(natoms)  - Atomic number of all atoms (X=0, GH=110)
c   pfullclss(fullordr)- Class type vector
c   pcompclss(compordr)-
c   pfullpopv(natoms)  - Number of atoms in each orbit
c   pcomppopv(natoms)  -
c   pfullmemb(natoms)  - Atoms sorted by point group orbits
c   pcompmemb(natoms)  -
c   pnprimatom(natoms) - Number of primitive functions for each atom
c   pnshellatom(natoms)- Number of different angular momentum shells for each
c                        atom (takes on values of 1,4,9,16, etc.)
c   pnangatom(natoms)  - The number of different angular momentum for each
c                        atom (takes on values of 1,2,3,4, etc.)
c   pnaoatom(natoms)   - Number of AOs for each atom
c   pnshellprim(maxshell,natoms)
c                      - The number of primitive functions in each shell
c                        of each atom
c   pnshellao(maxshell,natoms)
c                      - The number of AOs in each shell of each atom
c   pprimoff(maxshell,natoms)
c   paooff(maxshell,natoms)
c                      - The primcoef matrix is a block diagonal matrix.
c                        Each shell of each atom has a block.  If you have
c                        a list of all primitive functions, pprimoff(ishell,
c                        iatom) tells the location of the first primitive
c                        function in the block (ishell,iatom) and paooff
c                        contains similar information for the AOs.
c
c ***NOTE***  Because joda stores pfullpopv/pcomppopv as size natoms, we
c             do to, but they should be of size fullnorb/compnorb.  The
c             first ones have real values.  The remaining ones are 0.

      double precision orientmt(3,3),nucrep
      integer natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb(2),
     &    nvrtorb(2),maxshlao,maxangmom,natomsx
      integer patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,
     &    zalpha,zprimcoef,pprimoff,paooff,pnangatom

      common /mol_com/ orientmt,nucrep,
     &    natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb,
     &    nvrtorb,maxshlao,maxangmom,natomsx,
     &    patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,zalpha,
     &    zprimcoef,pprimoff,paooff,pnangatom
      save   /mol_com/

      character*4 compptgp,fullptgp
      character*1 spinc(2)
      common /molc_com/ compptgp,fullptgp,spinc
      save   /molc_com/




c This contains information about each of the possible grids for
c performing the numerical integration.

c###########################################################################
c MISC
c###########################################################################
c maxgrdatm : The largest atomic number for which the Slater and Bragg-Slater
c             atomic size is known.
c atmrad    : Atomic size using Slater's' rules for the radial integration
c xbsl      : The Bragg-Slater radii (one for each atom)
c             ***NOTE*** This is a genuine constraint.  Only atoms smaller
c             then this (currently 86) may be calculated.
c numangfct : number of angular momentum functions
c minpt     : The number of points used in the integration over
c             interatomic paths to find the minimum density point
c             between two atoms
c pangfct   : a pointer to the array of x,y,z angular momentum for each
c             angular momentum function

      integer maxgrdatm
      parameter (maxgrdatm=86)
      integer minpt
      parameter (minpt=100)

      integer numangfct

      double precision
     &    atmrad(maxgrdatm),xbsl(maxgrdatm),
     &     TA(maxgrdatm),multiEX(maxgrdatm)
      common /grid/  numangfct
      save /grid/

      common /gridd/ atmrad,xbsl,TA,multiEX
      save /gridd/

      integer pangfct
      common /gridp/ pangfct
      save /gridp/

c###########################################################################
c RADGRD file
c###########################################################################
c maxanggrid: The maximum number of different angular grids which can be
c             used in any given calculation.
c             ***NOTE*** This is a genuine constraint, but it must be used
c             since we must be able to keep a record of which angular grids
c             are used (before we have any allocated memory) since we have
c             to know how many grids are used in order to determine how much
c             memory to allocate.  This is set high enough it should never
c             be a problem.
c gridlist  : A list of all grids used (see the comment on maxanggrid).  It
c             is of dimension (maxanggrid,3) to keep track of the type and
c             subtype, and the number of times each grid is used.  The type
c             refers to how the grid is arrived at.
c                1 : Lebedev
c                2 : ?
c             The subtype refers to the degree of the grid of this type.
c numgrid   : The number of different angular grids used in the calculation.
c maxangpts : The maximum number of points in any of the angular grids used.
c maxanggrd : The grid with the maximum number of angular points.
c numradpts : The number of different radial points.
c ntotrad   : The total number of points in all angular grids at all radial
c             points (i.e. the entire integration grid)
c
c iradint   : determines if the Handy method (1) or Gauss-Legendre (2)
c             radial integration is used
c autosiz   : A flag which sets whether the polyhedra are (1) equally
c             sized, (2) sized according to Bragg-Slater radii or
c             (3) automatically sized according to the minimums in
c             density.
c slater    : A flag which determines whether (0) Slater's' rules are
c             used to determine the atomic size and scale the radial
c             integration or (1) no scaling is used.
c rigid     : A flag which determines whether rigid (0) or fuzzy
c             partitioning is used.
c nitr      : The number of iterations of the equations which create the
c             'fuzzy' boundary.

      integer maxanggrid
      parameter (maxanggrid=1000)

      integer gridlist(maxanggrid,3),numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad

      common /radgrd/  gridlist,numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad
      save /radgrd/

c Memory pointers
c
c pradgrid(numradpts) : The angular grid to use at each radial point.
c pgrdangpts(numgrid) : The number of angular points in each grid.
c zgridxyz(3,maxangpts,numgrid)
c                     : The x,y,z coordinate of each angular point in each grid
c zgridwt(maxangpts,numgrid)
c                     : The weight at each point.
c pintegaxis(natoms,3): Contains information about how much of each axis to
c                       integrate over.  If integaxis(iatom,i) is set to i,
c                       integrate only over the positive half of the i^th
c                       axis.  Otherwise, integrate over the entire axis.

      integer pgrdangpts,pradgrid,zgridxyz,zgridwt,pintegaxis

      common /radgrdp/ pgrdangpts,pradgrid,zgridxyz,zgridwt,
     &    pintegaxis
      save /radgrdp/

c###########################################################################
c Old stuff
c###########################################################################

c polist    : Contains an ordered list of unique atoms
c zatmvc    : The x, y, and z distance between each pair of atoms.
c zrij      : The distance between each pair of atoms.
c zatmpth   : The cartesian coordinates for the path integration between
c              all atom pairs
c zptdis    : The distance from atom i to a point along the path between
c              atoms i and j
c zprsqrd   : The distance squared from each atom to a point along all the
c              paths between all the atoms
c zpthpt    : The cartesian coordinates with respect to each atom for
c              the points along all the paths between all the atoms
c zbslrd    : The Bragg-Slater radii.
c zaij      : Surface shifting parameter dependent on the distance between
c               pairs of atoms.

      integer polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,zptdis,
     &    zprsqrd,zpthpt,zbslrd,zaij
      common /gridold/ polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,
     &    zptdis,zprsqrd,zpthpt,zbslrd,zaij
      save /gridold/




c This commonblock contains values for the functionals involved in the
c numerical integration either for DFT, plotting or fitting

c totele    : The total number of electrons
c etotekin  : The total kinetic energy
c etotenatr : The total nuclear-electron attraction energy
c etottf    : The total Thomas Fermi kinetic energy
c etotw     : The total Weizacker kinetic energy
c xldax     : The total LDA exchange energy
c becke     : The total Becke exchange energy
c lda       : The total LDA correlation energy
c xlyp      : The total LYP correlation energy
c icntr     : The integration center
c idns      : a flag =0 for SCF orbitals and =1 for natural orbitals

      integer icntr,idns
      double precision
     &    totele,etotkin,etotnatr,etottf,etotw,xldax,becke,
     &    lda,xlyp

      common /int_com/  icntr,idns
      save   /int_com/
      common /intr_com/ totele,etotkin,etotnatr,etottf,etotw,xldax,
     &                  becke,lda,xlyp
      save   /intr_com/

c array pointers

c zpcoeff(2): alpha/beta MO to primitive function transformation matrix
c zxocc     : alpha/beta orbital occupation

      integer
     &    zpcoeff(2),zxocc
      common /molecp/
     &    zpcoeff,zxocc
      save /molecp/



      integer
     &    atomchrg(natoms),radgrid(int_numradpts),
     &    integaxis(3,natoms),grdangpts(numgrid),
     &    iradpt,iangpt,grid
      double precision
     &    cdnt(natoms,3),atmvc(natoms,natoms,3),
     &    rsqrd(natoms),rij(natoms,natoms),
     &    aij(natoms,natoms),wtintr(natoms),totwt,
     &    rrtmp(natoms),radpt(int_numradpts),
     &    rwt(int_numradpts,maxangpts),
     &    gridxyz(3,maxangpts,numgrid),gridwt(maxangpts,numgrid)
      logical evalpt

      integer 
     &    i,n

      double precision
     &    atmwt,fctr,xx(0:3)

c ********************************************************************
c ********************************************************************

      call callstack_push('SYMOCT')

c Determine the angular weights and the cartesian coordiantes of each point
      xx(0)=zero

c Determine the cartesian coordinates of the point with respect to the
c integration center
      do 185 i=1,3
        xx(i)=radpt(iradpt)*gridxyz(i,iangpt,grid)
  185 continue

c Only keep points which fall in the xyz octant
      if (xx(integaxis(1,icntr)).ge.zero.and.
     &    xx(integaxis(2,icntr)).ge.zero.and.
     &    xx(integaxis(3,icntr)).ge.zero) then
        
        call crtgrid(icntr,atomchrg,atmvc,rij,aij,cdnt,rsqrd,
     &      rrtmp,wtintr,atmwt,xx(1))

c Determine total weight
        totwt=pi*four*atmwt*rwt(iradpt,iangpt)*
     &      gridwt(iangpt,grid)

c Multiply weight by number of symmetry equivalent points
        n=0
        if (xx(integaxis(1,icntr)).ne.zero) n=n+1
        if (xx(integaxis(2,icntr)).ne.zero) n=n+1
        if (xx(integaxis(3,icntr)).ne.zero) n=n+1
        fctr=two**n

        totwt=fctr*totwt

      else
        evalpt=.false.
      endif

      call callstack_pop
      return
      end
