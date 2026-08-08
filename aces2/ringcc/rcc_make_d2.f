













































































































































































































      Subroutine Rcc_make_d2(E,D,Iwork,Nbasis,Isize,Laabb,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Dimension E(2*Nbasis),D(Laabb)
      Dimension Iwork(Isize)
      Integer A,B,A1

      Logical Rle
      COMMON/EXTRAPO/RLE



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

      NNP1O2(N)=(N*(N+1))/2
      NNM1O2(N)=(N*(N-1))/2 
      ILRG(IX) = INT(0.5*(1.0D0+SQRT(8.0D0*IX-7))+1.0D-08)+1 


      Call Getrec(20,"JOBARC","SCFEVALA",Nbasis*Iintfp,E)
      Call Getrec(20,"JOBARC","SCFEVALB",Nbasis*Iintfp,E(Nbasis+1))

      Irrepx = 1
      Call Newtyp2(Irrepx,48,5,7,.True.)
      Call Newtyp2(Irrepx,64,5,7,.True.)
      If (Iuhf .Ne. 0) Call Newtyp2(Irrepx,49,6,8,.True.)
      If (Iuhf .Ne. 0) Call Newtyp2(Irrepx,65,6,8,.True.)
      
      Do Ispin = 1, Iuhf+1
         Ioff = (Ispin-1) * Nbasis 
         Nrow = Nnp1o2(Nvrto(Ispin))
         Ncol = Nnp1o2(Nocco(Ispin))
         I000 = 1
         I010 = I000 + Nirrep
         I020 = I010 + Ncol 
         I030 = I020 + Nirrep
         I040 = I030 + Nrow 

         If (Ispin .EQ. 1) Then
            Call Getrec(20,"JOBARC",'SOAOA1X ',Nirrep,Iwork(I000))
            Call Getrec(20,"JOBARC",'SOAOA1  ',Ncol,  Iwork(I010))
            Call Getrec(20,"JOBARC",'SVAVA1X ',Nirrep,Iwork(I020))
            Call Getrec(20,"JOBARC",'SVAVA1  ',Nrow,  Iwork(I030))
         Else
            Call Getrec(20,"JOBARC",'SOBOB1X ',Nirrep,Iwork(I000))
            Call Getrec(20,"JOBARC",'SOBOB1  ',Ncol,  Iwork(I010))
            Call Getrec(20,"JOBARC",'SVBVB1X ',Nirrep,Iwork(I020))
            Call Getrec(20,"JOBARC",'SVBVB1  ',Nrow,  Iwork(I030))
         Endif 

          Ioffl = 0
          Ioffr = 0
          Do Irrep = 1, Nirrep
             Noo    = Iwork(I000+Irrep-1)
             Nvv    = Iwork(I020+Irrep-1)
             Icount = 0

             Do Ioo = 1, Noo
                Ioo_loc = Ioo + Ioffr
                Ioo_off = Iwork(I010+Ioo_loc-1) 
                I       = Ilrg(Ioo_off) - 1
                I1      = Ilrg(Ioo_off)
                J       = Ioo_off - Nnm1O2(I1-1) 
                Do Ivv = 1, Nvv
                   Icount = Icount + 1
                   Ivv_loc = Ivv + Ioffl
                   Ivv_off = Iwork(I030+Ivv_loc-1) 
                   A       = Nocco(Ispin) + Ilrg(Ivv_off) - 1
                   A1      = Ilrg(Ivv_off)
                   B       = Nocco(Ispin) + Ivv_off - Nnm1O2(A1-1) 
                   D(Icount) = 1.0D0/(E(Ioff+I)+E(Ioff+J)-
     +                                E(Ioff+A)-E(Ioff+B))
                Enddo
             Enddo 
             Call Putlst(D,1,Noo,1,Irrep,47+ispin)

             If (Rle) Then
                Do Iinv = 1, Icount
                    D(Iinv) = 1.0D0/D(Iinv)
                Enddo 
                Call Putlst(D,1,Noo,1,Irrep,63+ispin)
             Endif 

             Ioffr = Ioffr + Noo
             Ioffl = Ioffl + Nvv
          Enddo
      Enddo

      Return
      end
