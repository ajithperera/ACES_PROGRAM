










      Subroutine Pccd_form_htau_2d_pppp(Hvv_pq,Hvv_qp,Work,Maxcor,Nvrt,
     +                                  Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hvv_pq(Nvrt*Nvrt)
      Dimension Hvv_qp(Nvrt*Nvrt)

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      Data One,Two,Onem,Dnull,Ione /1.0D0,2.0D0,-1.0D0,0.0D0,1/

C V(ce,da)*Gamma(ce,db) (ignore the labels).

      Call Dzero(Hvv_qp,Nvrt*Nvrt)
      Call Dzero(Hvv_pq,Nvrt*Nvrt)

      Irrepx = Ione
      Do Irrep_ce = 1, Nirrep
         Irrep_db = Dirprd(Irrep_ce,Irrepx)
         Irrep_da = Dirprd(Irrep_ce,Irrepx)

         Nrow_ce = Irpdpd(Irrep_ce,19)
         Ncol_db = Irpdpd(Irrep_db,19)
         Nrow_da = Irpdpd(Irrep_da,19)
         Ncol_ce = Irpdpd(Irrep_ce,19)

         I000 = Ione
         I010 = I000 + Nrow_ce*Ncol_db
         I020 = I010 + Max(Ncol_db,Nrow_ce)
         I030 = I020 + Max(Ncol_db,Nrow_ce)
         Iend = I030 + Max(Ncol_db,Nrow_ce)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_pppp",
     +                                      Iend,Maxcor)
 
         Call Getlst(Work(I000),1,Ncol_db,1,Irrep_db,List_v)

         Call Spinad3(Irrep_ce,Vrt(1,1),Nrow_ce,Ncol_db,Work(I000),
     +                Work(I010),Work(I020),Work(I030))

         Iend = I010 + Nrow_da*Ncol_ce
         Call Getlst(Work(I010),1,Ncol_ce,2,Irrep_ce,List_g)

C V(ce,da)(t)*G(ce,db) -> S(b,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_a = 1, Nirrep
            Irrep_d  = Dirprd(Irrep_a,Irrep_da)
            Irrep_b  = Dirprd(Irrep_d,Irrep_db)

            Nb = Vrt(Irrep_b,1)
            Na = Vrt(Irrep_a,1)
            Nd = Vrt(Irrep_d,1)

            Nsum = Nrow_ce*Nd
            Nrow = Nb
            Ncol = Na
            Icheck = Min(Nsum,Nrow,Ncol)

            If (Icheck .Gt .0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                      Nsum,Work(Koff),Nsum,One,Hvv_qp(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Na*Nb
            Joff = Joff + Nsum*Nb
            Koff = Koff + Nsum*Na
 
         Enddo
      Enddo

      call pccd_check_htau("Htau_vv :",Hvv_qp,Nvrt,"VV","D")
CSSS      call checksum("Htau_vv :",Hvv_qp,Nvrt*Nvrt)
      Return
      End
