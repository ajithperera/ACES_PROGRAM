






































































































































































































      Subroutine Drive_lhbar(Work, Length, Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Length)
      Dimension Ecorr(3)
      Logical Sing 
      Logical MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      Integer AAAA_LENGTH_IJKA



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c files.com : begin
      integer        luout, moints
      common /files/ luout, moints
c files.com : end

       Common /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
       Common /Energy/Energy(500,2),Ixtrle(500)

       Sing   = .False.
       Icontl = Iflags(4) 
       Iconvg = 1
       Ncycle = 0
       Sing   = (Iflags(2) .gt. 9)
       Luout  = 6
       Moints = 50

       Coulomb = .False.
       Call Parread(iuhf)
           write(6,2010) paralpha
           write(6,2011) parbeta
           write(6,2012) pargamma
           write(6,2013) pardelta
           write(6,2014) parepsilon
 2010      format(' PCCSD   alpha parameter : ', F14.6)
 2011      format(' PCCSD    beta parameter : ', F14.6)
 2012      format(' PCCSD   gamma parameter : ', F14.6)
 2013      format(' PCCSD   delta parameter : ', F14.6)
 2014      format(' PCCSD epsilon parameter : ', F14.6)
       If (ispar) Then
          write(6,"(a,a)") ' Perform a parameterized CC HBAR',
     &                    ' calculations'
          Write(6,*)
          if (coulomb) write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
          write(6,*)
          Fae_scale    = (Paralpha - 1.0D0)
          Fmi_scale    = (Parbeta  - 1.0D0)
          Wmnij_scale  = Pargamma
          Wmbej_scale  = Pardelta
          Gae_scale    = Paralpha  
          Gmi_scale    = Parbeta  
      Else
          write(6,*) '   Perform a regular CC HBAR calculations'
          write(6,*)
          Fae_scale    = 0.0D0
          Fmi_scale    = 0.0D0
          Wmnij_scale  = 1.0D0
          Wmbej_scale  = 1.0D0
          Gae_scale    = 1.0D0
          Gmi_scale    = 1.0D0
      Endif
C
C Unfortunately I have to use the METH common block. At the moment 
C We will set to CCSD uing calc key-word. 
C
       CCSD = (Iflags(2) .eq. 10)
C
C Copy the converged T1 and T2 amplitudes to L1 and L2 lists.
C These are our starting L1 and L2 amplitudes.
C
       Call init2_leom(Iuhf, Sing)
       Call Drmove(Work, Length, Iuhf, 100, Sing)
       Call Inilam(Iuhf)
       Call Rnabij(Work,Length,iuhf,"L")
       Call Diislst_leom(1, Iuhf, Sing)

       Irrepx = 1
      Write(6,*) "The Hbar lists"
      Call checkhbar(Work,length,Iuhf)
       Do While (Iconvg .eq. 1)

          Ncycle = Ncycle + 1
     
          Call Lguess(Work, Length, Iuhf, Sing)

      Write(6,*) "The lambda residuals before multiplication"
      Call check_leom(Work,length,Iuhf)
          Call Hbar_mult_l(Work, Length, Iuhf)

          Call Newl1_leom(Work,Length,Iuhf)
          Call Newt2(Work,Length,Iuhf)

      Write(6,*) "The lambda residuals after  multiplication"
      Call check_leom(Work,length,Iuhf)
          Call Drtsts(Work, Length, Ncycle, Iuhf, Iconvg, Icontl,
     &                Sing, 100, "L") 
          Call Dodiis0_l(Work, Length, Iuhf, 1, Ncycle, Iconvg,
     &                   Icontl, Sing, 144, 61, 190, 0, 90, 
     &                   2, 70, "       ")
          Call Drmove(Work, Length, Iuhf, 100, Sing)
          Call Rnabij(Work,Length,iuhf,"L")

          Call Cmpeng2(Work,Length,60,2,Ecorr,Energy(Ncycle+1,1),
     &                 Energy(Ncycle+1,2),Iuhf,1)

       Enddo 

       Return
       End
