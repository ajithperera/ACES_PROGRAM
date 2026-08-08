






































































































































































































      Subroutine post_dccl_mods(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Length)



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

      Imode  = 0
      Irrepx = 1
C
C Needs the G-intermediates scales by Gae_ and G_mi scales.
C
      IRREPX      = 1
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)
      IPPA_LENGTH = IRPDPD(IRREPX,19)
      IPPB_LENGTH = IRPDPD(IRREPX,20)

      Call Updmoi(1,IHHA_LENGTH,9,191, 0, 0)
      Call Updmoi(1,IPPA_LENGTH,9,192, 0, 0)
      If (Iuhf .ne.0) Then
         Call Updmoi(1,IHHB_LENGTH,10,191, 0, 0)
         Call Updmoi(1,IPPB_LENGTH,10,192, 0, 0)
      Endif

      Call zerlst(Work,IHHA_LENGTH,1,1,9,191)
      Call zerlst(Work,IPPA_LENGTH,1,1,9,192)
      If (iuhf .ne. 0) then
            call zerlst(Work,IHHB_LENGTH,1,1,10,191)
            call zerlst(Work,IPPB_LENGTH,1,1,10,192)
      Endif
CSSS      If (Coulomb) Then
CSSS         Call Pdcc_formg1(Work,Length,Iuhf,1.0D0,8)
CSSS      Else 
CSSS        Call Pdcc_formg1(Work,Length,Iuhf,Gae_scale,8)
CSSS      Endif 
      Call Pdcc_formg1(Work,Length,Iuhf,Gae_scale,8)
      Call Pdcc_formg2(Work,Length,Iuhf,Gmi_scale,8)

C Store a flag to indicate that lambda has been done. This flag
C is used in EOM code to test whether lambad is done. 

      Call Putrec(20,"JOBARC","LAMBDA  ",1,Length)


      Return
      End
 
