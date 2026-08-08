













































































































































































































      Subroutine Rcc_make_t2(Work,Maxcor,Iuhf)

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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

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
      Write(6,"(a)") "--------Entering rcc_make_t2---------"
      Write(6,*)

      Irrepx = 1
      Imod   = 0

C Form lists to store <AB|IJ> integrals (ordered as AB;IJ)
C List 1 and 2.

      Call Inipck(Irrepx,19,21,1,Imod,0,1)
      If (Iuhf .Ne.0) Call Inipck(Irrepx,20,22,2,Imod,0,1)

       Eaa = 0.0D0
       Ebb = 0.0D0
       Do Ispin = 1, 1+Iuhf 
       
          Nsize = Idsymsz(Irrepx,Isytyp(1,4+Ispin),Isytyp(2,4+Ispin))
          Nrows = Nvrto(Ispin) * Nvrto(Ispin) +
     +            Nocco(Ispin) * Nocco(Ispin) +
     +            Nvrto(Ispin) * Nocco(Ispin) 

          I000 = 1
          I010 = I000 + Nsize
          I020 = I010 + Nsize 
          I030 = I020 + Nrows 
          Iend = I030
          If (Iend .Ge. Maxcor) Call Insmem("rcc_make_t2",Iend,
     +                                       Maxcor)
C Read in Coulomb integrals as Ispin=1: W(AI,BJ) Ispin=2 (ai,bj)

          Call Getall(Work(I010),Nsize,Irrepx,4+Ispin)

C Form the integrals (AB,IJ)by permuting I and B;Ispin=1: W(AB,IJ)
C Ispin=2 (ab,ij)

          Call Sstgen(Work(I010),Work(I000),Nsize,Vrt(1,Ispin),
     +                Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),
     +                Work(I020),Irrepx,"1324")

          Call Putall(Work(I000),Nsize,Irrepx,Ispin)

          Do Irrepr = 1, Nirrep
             Irrepl = Dirprd(Irrepr,Irrepx)

             Nrow_tring_ab = Irpdpd(Irrepl,4+Ispin)
             Ncol_tring_ij = Irpdpd(Irrepr,6+Ispin)
             Nrow_squar_ab = Irpdpd(Irrepl,18+Ispin)
             Ncol_squar_ij = Irpdpd(Irrepr,20+Ispin)

             I000 = 1
             I010 = I000 + Nrow_tring_ab * Ncol_tring_ij
             I020 = I010 + Nrow_squar_ab * Ncol_squar_ij
             I030 = I020 + Nrow_squar_ab * Ncol_squar_ij
             I040 = I030 + Nrow_squar_ab * Ncol_squar_ij
             Iend = I040

             If (Iend .Ge. Maxcor) Call Insmem("rcc_make_d2",Iend,
     +                                         Maxcor)
C Retrive D(A<=B,I<=J) and expand,

             Call Getlst(Work(I000),1,Ncol_tring_ij,1,Irrepr,47+Ispin)

C D(A<=B,I<=J) -->D(AB,I<=J)
             
             Call Symexp6(Irrepl,Vrt(1,ispin),Vrt(1,ispin),
     +                    Nrow_squar_ab,Nrow_tring_ab,
     +                    ncol_tring_ij,work(I010),work(I000),
     +                    work(I030))
             Call Transp(work(I010),work(I020),ncol_tring_ij,
     +                        nrow_squar_ab)
C D(IJ,A<=B) -->D(iJ,AB)

             Call Symexp6(Irrepr,Pop(1,ispin),Pop(1,ispin),
     +                    Ncol_squar_ij,Ncol_tring_ij,
     +                    nrow_squar_ab,work(I010),work(I020),
     +                    work(i030))
C D(IJ,AB) -->D(AB,IJ)

             Call Transp(work(I010),work(I020),Nrow_squar_ab,
     +                        Ncol_squar_ij)
C Retrive integrals I(AB;IJ) array

             Call Getlst(Work(I010),1,Ncol_squar_ij,1,Irrepr,Ispin)

             Do I = 1,Nrow_squar_ab * Ncol_squar_ij
                Work(I030+I-1)  = Work(I020+I-1) * Work(I010+I-1)
             Enddo

             Call Putlst(Work(I030),1,Ncol_squar_ij,1,Irrepr,
     +                   43+Ispin)

             If (Ispin .EQ. 1) Then
             Eaa = Eaa + Ddot(Nrow_squar_ab*Ncol_squar_ij,Work(I030),1,
     +                        Work(I020),1)
             Else
             Ebb = Ebb + Ddot(Nrow_squar_ab*Ncol_squar_ij,Work(I030),1,
     +                        Work(I020),1)
             Endif

          Enddo
      Enddo

      Return
      End
c-------------------------------------------------------------------

