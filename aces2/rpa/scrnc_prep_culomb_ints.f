










      Subroutine Scrnc_prep_culomb_ints(Coulmb_oo,Coulmb_vv,Coulmb_vo,
     +                                  Iuhf,Len_oo,Len_vv,Len_vo,
     +                                  Ipq_pair_4irep,Irrepx,Irrepr)

      Implicit Double Precision (A-H,O-Z)
      Integer Ipq_pair_4irep

      Dimension Coulmb_oo(Len_oo),Coulmb_vv(Len_vv),Coulmb_vo(Len_vo)



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

      List_aixx   = 205
      List_abxx   = 207
      List_ijxx   = 209

      Ioff_oo = 1
      Ioff_vv = 1
      Ioff_vo = 1
      Write(6,"(a,a,i5,i1)") " Integrals are obtained for pair",
     +                      " and irrep:",
     +                     Ipq_pair_4irep,Irrepr 
     
      Do Ispin = 1, Iuhf+1
         Call Getlst(Coulmb_oo(Ioff_oo),Ipq_pair_4irep,1,1,Irrepr,
     +               List_ijxx+Ispin)
         Call Getlst(Coulmb_vv(Ioff_vv),Ipq_pair_4irep,1,1,Irrepr,
     +               List_abxx+Ispin)
         Call Getlst(Coulmb_vo(Ioff_vo),Ipq_pair_4irep,1,1,Irrepr,
     +               List_aixx+Ispin)

         Ioff_oo = Ioff_oo + Irpdpd(Irrepx,21)
         Ioff_vv = Ioff_vv + Irpdpd(Irrepx,19)
         Ioff_vo = Ioff_vo + Irpdpd(Irrepx,9)
      Enddo 

      Write(6,*)
      write(6,"(a)") " <oo|J|XX>,<vv|J|XX> and <vo|J|XX> integrals"
      Ioff_oo = 1
      Ioff_vv = 1
      Ioff_vo = 1
      Lenoo  = Irpdpd(Irrepx,21)
      Lenvv  = Irpdpd(Irrepx,19)
      Lenvo  = Irpdpd(Irrepx,9)
      Write(6,"(a,3(1x,i3))")" Lenoo,Lenvv,Lenvo:", Lenoo,Lenvv,Lenvo
      call checksum("<oo|j|xx> Alpha: ",Coulmb_oo(ioff_oo),Lenoo)
      call checksum("<vv|j|xx> Alpha: ",Coulmb_vv(ioff_vv),Lenvv)
      call checksum("<vo|j|xx>:Alpha: ",Coulmb_vo(ioff_vo),Lenvo)
      Ioffoo = Ioff_oo + Irpdpd(Irrepx,21)
      Ioffvv = Ioff_vv + Irpdpd(Irrepx,19)
      Ioffvo = Ioff_vo + Irpdpd(Irrepx,9)
      call checksum("<oo|j|xx> Beta : ",Coulmb_oo(ioff_oo),Lenoo)
      call checksum("<vv|j|xx> Beta : ",Coulmb_vv(ioff_vv),Lenvv)
      call checksum("<vo|j|xx>:Beta : ",Coulmb_vo(ioff_vo),Lenvo)
      Return
      End
