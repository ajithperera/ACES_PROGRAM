











      Subroutine Scrnc_bcktran_left(JJ_respn_oo,JJ_respn_xx,
     +                              Eveca,JJ_respn_xo,
     +                              Nij_rows,Nxj_rows,Nxx_rows,
     +                              Nao,Nmo,Irrep_l)

      Implicit Double Precision (A-H,O-Z)

      Double Precision JJ_respn_oo(Nij_rows)
      Double Precision JJ_respn_xo(Nxj_rows)
      Double Precision JJ_respn_xx(Nxx_rows)
      Dimension Eveca(Nao,Nmo)
      Double Precision Null

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
    
      Data Null,One /0.0D0,1.0D0/

      Call Getrec(20,"JOBARC","SCFEVECA",Nao*Nmo*Iintfp,Eveca)

      Ioff_ji = 1
      Ioff_xx = 1
 
      Do Irrep_i = 1, Nirrep 
         Irrep_j = Dirprd(Irrep_i,Irrep_l)

         Nsum = Pop(Irrep_j,1)
         Nrow = Iaopop(Irrep_j)
         Ncol = Pop(Irrep_i,1)
        
         Ioff_x = Ioffao(Irrep_j)
         Ioff_j = Pop(Irrep_j,1)

         If (Min(Nrow,Ncol,Nsum) .NE. 0) Then
            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Eveca(Ioff_x,Ioff_j),
     +                  Nao,JJ_respn_oo(Ioff_ji),Nsum,Null,
     +                  JJ_respn_xo,Nrow)
         Else
            Call Dzero(JJ_respn_xo,Nrow*Ncol)
         Endif 
         
         Nsum = Pop(Irrep_j,1)
         Ncol = Iaopop(irrep_j)
         Nao2 = Ncol

         Ioff_x = Ioffao(Irrep_i)
         Ioff_i = Pop(Irrep_i,1)

         If (Min(Nrow,Ncol,Nsum) .NE. 0) Then
            Call Xgemm("N","T",Nrow,Ncol,Nsum,One,JJ_respn_xo,Nrow,
     +                  Eveca(Ioff_x,Ioff_i),Nao,Null,
     +                  JJ_respn_xx(ioff_xx),Nrow)
         Else
            Call Dzero(JJ_respn_xx(ioff_xx),Nrow*Ncol)
         Endif 

         Ioff_ji = Ioff_ji + Nsum
         Ioff_xx = Ioff_xx + Nrow * Ncol

       Enddo

       Return
       End
       
