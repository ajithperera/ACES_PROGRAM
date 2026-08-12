










      Subroutine Scrnc_sympack_aoints(W_in,W_out,Nbfns,Irrepx)

      Implicit Double Precision(A-H,O-Z)
      Dimension W_in(Nbfns,Nbfns), W_out(Nbfns*Nbfns)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
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

      Ndim2 = Nbfns * Nbfns 
      Call Dzero(W_out,Ndim2)

      Ithru=0
      Do Irrepr=1, Nirrep
         Irrepl=Dirprd(Irrepr,Irrepx)
          Do Nu=1, Iaopop(Irrepr)
             Do Mu=1,Iaopop(Irrepl) 
                Ithru = Ithru + 1
                Mu_off = Mu + IoffAO(Irrepl)-1 
                Nu_off = Nu + IoffAO(Irrepr)-1 
                W_out(Ithru) = W_in(Mu_off,Nu_off)
               Enddo
            Enddo
         Enddo 


      Return
      End
