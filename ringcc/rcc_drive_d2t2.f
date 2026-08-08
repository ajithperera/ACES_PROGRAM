










      Subroutine Rcc_drive_d2t2(Work,Maxcor,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Maxcor)



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

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Write(6,*)
      Write(6,"(a)") "--------Entering rcc_drive_d2t2---------"
      Write(6,*)
      
C Obtain the maximum length of the Denominator AA and BB array

      Maxaabb = 1
      Do Ispin = 1, Iuhf+1
         Do Irrep = 1, Nirrep
            Nrow = Irpdpd(Irrep,4+Ispin)
            NCol = Irpdpd(Irrep,6+Ispin)
            Length = Nrow * Ncol
            If (Length .LT. 0) Call Trap_intovf("Rcc_make_d2t2",1)
            Maxaabb=Max(Maxaabb,Length)
         Enddo 
      Enddo 

      Nbasis = Nocco(1) + Nvrto(1)
      Itmp1  = Nocco(1) * Nocco(1)
      Itmp1  = Nvrto(2) * Nvrto(2)
      Isize  = Itmp1 + Itmp2 + 2*Nirrep 
C  
      I000 = 1
      I010 = I000 + 2*Nbasis
      I020 = I010 + Maxaabb
      I030 = I020 + Isize
      I030 = Iend 
      If (Iend .LT. 0) Call Trap_intovf("Rcc_make_d2t2",2)
      If (Iend .GT. Maxcor) Call Insmem("Rcc_make_d2t2",Iend,Maxcor)
      Call Rcc_make_d2(Work(I000),Work(I010),Work(I020),Nbasis,
     +                 Isize,Maxaabb,Iuhf)

      Call Rcc_make_t2(Work,Maxcor,Iuhf)
     
      Return
      End 
