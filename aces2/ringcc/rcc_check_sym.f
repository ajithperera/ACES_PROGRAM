










      Subroutine rcc_check_sym(Work,Maxcor,Ispin,List,Type,flag)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Character*6 Flag



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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
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

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Character*4 type

      Write(6,"(a,a,i4,a,1x,6a)") " Antisymmetry/symmetry check",
     +                            " of list",List,
     +                            " and type",type
      Irrepx = 1
      Ioff   = 1

      If (Flag .Eq. "noread") Then

      Do Irrep_r = 1, Nirrep

         Irrep_l = Dirprd(Irrep_r, Irrepx)

         Nrow  = Irpdpd(Irrep_l,20+ispin)
         Ncol  = Irpdpd(Irrep_r,20+Ispin)

         Do icol = 1, Ncol  
         Write(6,*) 
         Write(6,"(a,i4)") " Printing <AB|IJ> for IJ pair: ", icol

         do irrep_sr =1, Nirrep
            Irrep_sl = Dirprd(Irrep_sr, Irrep_l)

            If (Type .eq. "pphh") Then
               Nrow_s = Vrt(irrep_sl,Ispin)
               Nrow_l = Vrt(irreP_sr,Ispin)

               call output(Work(ioff),1,Nrow_s,1,Nrow_l,
     +                     Nrow_s,Nrow_l,1)
           Endif
           Ioff = Ioff + Nrow_s*Nrow_l
         Enddo
         Enddo 
      Enddo

      Elseif (Flag .Eq. "doread") Then

      Do Irrep_r = 1, Nirrep

         Irrep_l = Dirprd(Irrep_r, Irrepx)

         Nrow  = Irpdpd(Irrep_l,Isytyp(1,List))
         Ncol  = Irpdpd(Irrep_r,Isytyp(2,List))
      
         Call Getlst(Work,1,Ncol,1,Irrep_r,List)

         Do irrep_sr =1, Nirrep
            Irrep_sl = Dirprd(Irrep_sr, Irrep_l)
 
            If (Type .eq. "pphh") Then
               Nrow_s = Vrt(irrep_sl,Ispin)
               Nrow_l = Vrt(irreP_sr,Ispin)
    
               call output(Work(ioff),1,Nrow_s,1,Nrow_l,
     +                     Nrow_s,Nrow_l,1)
           Endif
           Ioff = Ioff + Nrow_s*Nrow_l

         Enddo
      Enddo 

      Endif 

      Return
      End
