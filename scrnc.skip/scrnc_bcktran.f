










      Subroutine Scrnc_bcktran(Work,Maxmem,Iuhf,JJ_respn_ijkl,Npairs,
     +                         Irrepx)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxmem)
      Double Precision JJ_respn_ijkl(Npairs,Npairs)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)


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

      Write(6,"(1x,a)") "-----Entering scrnc_bcktran------"

      Write(6,"(a)") " The input JJ response matrx"
      call output(JJ_respn_ijkl,1,Npairs,1,Npairs,Npairs,Npairs,1)
      Write(6,*)
      Write(6,"(a)") " The screened <xx|xx> integrals" 

      Nmo  = Nocco(1)+Nvrto(1)
      Ione = 1
      Call Getrec(20,"JOBARC","NBASTOT",Ione,Nao)

      Do irrep_r = 1, Nirrep 
         Irrep_l = Dirprd(Irrep_r,Irrepx)

         Nij_cols = Irpdpd(Irrep_r,21)
         Nij_rows = Irpdpd(Irrep_l,21)
         Nxj_rows = Irpdpdaomo_occbk(Irrep_l,1)
         Nxj_cols = Irpdpdaomo_occbk(Irrep_r,1)
         Nxx_rows = Irpdpdao(Irrep_l)

         I000 = 1
         I010 = I000 + Nxx_rows * Nij_cols
         I020 = I010 + Nao * Nmo 
         I030 = I020 + Nxj_rows
         I040 = I030 + Nxx_rows * Nij_cols
         I040 = Iend
         if (Iend .GE. Maxmem) Call Insmem("Scrnc_backtrn",Iend,
     +                                      Maxmem)
         Ioff_xx_l = I000
         Do Icol = 1, Nij_cols

            Call Scrnc_bcktran_left(JJ_respn_ijkl(1,Icol),
     +                              Work(Ioff_xx_l),Work(I010),
     +                              Work(I020),Nij_rows,Nxj_rows,
     +                              Nxx_rows,Nao,Nmo,Irrep_l)
            Ioff_xx_l = Ioff_xx_l + Nxx_rows
         Enddo 

C      Write(6,"(a)") " The two index transformed JJ response matrx"
C      call output(Work(I000),1,Nxx_rows,1,Nij_cols,Nxx_rows,
C     +            Nij_cols,1)
      call checksum("(xx|IJ): ",Work(I000),Nxx_rows*Nij_cols)

C After completion of the above loop, we should have (xx|pq) batches of 
C integrals. Transpose it form (pq|xx) follwed by the loop to form (xx|xx).

         Ncols = Nij_cols
         Nrows = Nxx_rows
         Call Transp(Work(I000),Work(I030),Ncols,Nrows)
         Call Dcopy(Ncols*Nrows,Work(I030),1,Work(I000),1)

         Nxx_cols = Nxx_rows
         Nxx_rows = Irpdpdao(Irrep_r)
         Nxj_rows = Nxj_cols
         
         I020 = I010 + Nao * Nmo
         I030 = I020 + Nxx_rows * Nxx_cols 
         I040 = I030 + Nxj_rows 
         Iend = I040 
         if (Iend .GE. Maxmem) Call Insmem("Scrnc_backtrn",Iend,
     +                                      Maxmem)
         Ioff_xx_r = I000
         Ioff_xx_l = I020

         Do Icol = 1, Nxx_cols

            Call Scrnc_bcktran_left(Work(Ioff_xx_r),Work(Ioff_xx_l),
     +                              Work(I010),Work(I030),Nij_cols,
     +                              Nxj_rows,Nxx_rows,Nao,Nmo,Irrep_r)
            Ioff_xx_r = Ioff_xx_r + Nij_cols
            Ioff_xx_l = Ioff_xx_l + Nxx_rows

         Enddo 
C      Write(6,"(a)") " The two index transformed JJ response matrix"
C      call output(Work(I020),1,Nxx_rows,1,Nxx_cols,Nxx_rows,
C     +            Nxx_cols,1)
      call checksum("(xx|xx): ",Work(I020),Nxx_rows*Nxx_cols)

      Enddo
      
      Return
      End 
