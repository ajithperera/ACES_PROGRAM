













































































































































































































      Subroutine construct_t2wmnij(Work, Length, Iuhf)

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

      Imode = 0
      Write(6,"(a,a)") " Construction of T2 piece of the", 
     &                 " W(mn,ij) = 1/4t2(ij,ef)*<mn||ef>."
      Write(6,*)

      Call Inipck(1,14,14,253,Imode,0,1)

      If (Iuhf .ne. 0) Then
         Call Inipck(1,3,3,251,Imode,0,1)
         Call Inipck(1,4,4,252,Imode,0,1)
      Endif
C
C No need to do anything here for methods that do not have
C singles, but needs the 251-253 lists.
C
      If (Iflags(2) .lt. 9) Return
C
C Construct the W(mn,ij) intermediate as it would be for
C standard Cc calc.
C
      Intpck = 5
      Fact   = -(Pargamma - 1.0D0)

      Call Pdcc_quad1(Work,Length,Intpck,Iuhf,Fact)
      Call pdcc_drlad(Work,Length,Iuhf,1)

      Return 
      End 

