








































































































































































































      SUBROUTINE FORM_PERCENT_EXCITE(SCR, MAXCOR, IRREPX, IUHF,
     &                               ROOT, ECC, IROOT) 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION SCR(MAXCOR)

c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end


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




      Nbas = Nocco(1) + Nvrto(1)
  
      Iunit  = 556

C The unit 556 is opened in  vee.F and correspond to "EOM_excite"

      Itop   = 1
      Imap_a = Itop 
      Imap_b = Imap_a + Nbas
      Iend   = Imap_b + Nbas
      If (Iend .GE. Maxcor) CALL INSMEM("@-FORM_PERCENT_EXCITE",Iend,
     &                                   Maxcor)
      
      Call Getrec(20, "JOBARC", "ORBMAP_A", Nbas, Scr(Imap_a))
      Call Getrec(20, "JOBARC", "ORBMAP_A", Nbas, Scr(Imap_b))
      If (Iuhf .Ne. 0) Call Getrec(20, "JOBARC", "ORBMAP_B", Nbas, 
     &                             Scr(Imap_b))

      Write(Iunit,*)
      Write(Iunit,881) iroot, irrepx
 881  Format(' Root # ', i4, 3x, ' [', i1,']')
      Write(iunit,882) Root*27.2113957d0, Root, Root+Ecc
 882  Format( ' Excitation energy (eV): ', F12.6, ' in a.u. ', F14.8,
     &      /,  ' Total energy : ', F16.8)

      Call get_topamps(Scr(Iend),Maxcor-Iend,Scr(Imap_a),Scr(Imap_b),
     &                 Iuhf,Nbas,Irrepx,Iunit)

      Return
      End
     
