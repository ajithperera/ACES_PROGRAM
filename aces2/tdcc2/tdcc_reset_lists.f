













































































































































































































      Subroutine Tdcc_reset_lists(Work,Memleft,Iuhf,Irrepx)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_IJAB,BBBB_LENGTH_IJAB,AABB_LENGTH_IJAB
      Integer PHA_Length,PHB_Length,HHA_length,HHB_Length
      Integer PPA_Length,PPB_Length

      Dimension Work(Memleft)



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

C There are three type of amplitudes of singles and doubles dimension.
C Let us call them and resized the pointers as follows
C 
C 390,392,394 are singles lists of mu_0, mu_t, mu_d (dot) and mu_td
C amplitudes (i,a).

      PHA_Length = Irpdpd(Irrepx,9)
      PHB_Length = Irpdpd(Irrepx,10)
      Write(6,*) PHA_length,phb_length
      
      Call Aces_list_resize(1,390,PHA_Length)
      Call Aces_list_resize(1,392,PHA_Length)
      Call Aces_list_resize(1,394,PHA_Length)
      Call Aces_list_resize(1,396,PHA_Length)

      If (Iuhf .NE. 0) Then
         Call Aces_list_resize(2,390,PHB_Length)
         Call Aces_list_resize(2,392,PHB_Length)
         Call Aces_list_resize(2,394,PHB_Length)
         Call Aces_list_resize(2,396,PHB_Length)
      Endif

C The doubles lists for mubar, mubar_tilde and mubar_dot and
c mubar_tilde_dot. Originally these were initialized to irrep1,
C but can accomodate any irrep.  

      If (Iuhf .EQ. 0) Then
         Call Newtyp2(Irrepx,316,13,14,.True.)
         Call Newtyp2(Irrepx,326,13,14,.True.)
         Call Newtyp2(Irrepx,336,13,14,.True.)
         Call Newtyp2(Irrepx,346,13,14,.True.)
      Else
C Mubar
         Call Newtyp2(Irrepx,314,1,3,.True.)
         Call Newtyp2(Irrepx,315,2,4,.True.)
         Call Newtyp2(Irrepx,316,13,14,.True.)
C Mubar_tilde
         Call Newtyp2(Irrepx,324,1,3,.True.)
         Call Newtyp2(Irrepx,325,2,4,.True.)
         Call Newtyp2(Irrepx,326,13,14,.True.)
C Mubar_dot
         Call Newtyp2(Irrepx,334,1,3,.True.)
         Call Newtyp2(Irrepx,335,2,4,.True.)
         Call Newtyp2(Irrepx,336,13,14,.True.)
C Mubar_tilde_dot
         Call Newtyp2(Irrepx,344,1,3,.True.)
         Call Newtyp2(Irrepx,345,2,4,.True.)
         Call Newtyp2(Irrepx,346,13,14,.True.)
      Endif

C Temporary storage of (ai,bi) quantities 
    
      If (Iuhf .EQ. 0) Then
         Call Newtyp2(Irrepx,302,9,10,.True.)
      Else
         Call Newtyp2(Irrepx,300,9,9,.True.)
         Call Newtyp2(Irrepx,301,10,10,.True.)
         Call Newtyp2(Irrepx,302,9,10,.True.)
      Endif

C Lists for MO dipole moment integrals (OO,OV and VV blocks).

      HHA_Length = Irpdpd(Irrepx,21)
      HHB_Length = Irpdpd(Irrepx,22)
      PPA_Length = Irpdpd(Irrepx,19)
      PPB_Length = Irpdpd(Irrepx,20)

      Call Aces_list_resize(1,380,HHA_Length)
      Call Aces_list_resize(1,382,PPA_Length)
      Call Aces_list_resize(1,384,PHA_Length)

      If (Iuhf .NE. 0) Then
         Call Aces_list_resize(1,380,HHB_Length)
         Call Aces_list_resize(1,382,PPB_Length)
         Call Aces_list_resize(1,384,PHB_Length)
      Endif

      Return
      End
