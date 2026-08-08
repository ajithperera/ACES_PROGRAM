










      Subroutine Psi4dbg_form_htau_2d_hhpp_pp(Hvv_pq,Hvv_qp,Work,
     +                                        Maxcor,Nvrt,Nbas,
     +                                        List_v,List_g) 
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

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D00,1/

C V(cb,km)*Gamma(ca,km) (Gamma stored as ca,km and V stored as cb,km)
   
      Irrepx = Ione
      Call Dzero(Hvv_qp,Nvrt*Nvrt)
      Call Dzero(Hvv_pq,Nvrt*Nvrt)

      Do Irrep_km = 1, Nirrep
         Irrep_ca = Dirprd(Irrep_km,Irrepx)
         Irrep_cb = Dirprd(Irrep_km,Irrepx)

         Nrow_cb = Irpdpd(Irrep_cb,13)
         Ncol_km = Irpdpd(Irrep_km,14)
         Nrow_ca = Irpdpd(Irrep_ca,13)
         Ncol_km = Irpdpd(Irrep_km,14)

         I000 = Ione
         I010 = I000 + Nrow_cb*Ncol_km
         I020 = I010 + Nrow_cb*Ncol_km
         I030 = I020 + Max(Ncol_km,Nrow_cb)
         Iend = I030 + Max(Ncol_km,Nrow_cb)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_pp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_km,1,Irrep_km,List_v)
         Call Spinad3(Irrep_cb,Vrt(1,1),Nrow_cb,Ncol_km,Work(I010),
     +                Work(I020),Work(I030))
C V(cb,km) -> V(km,cb)
         Call Transp(Work(I010),Work(I000),Ncol_km,Nrow_cb)

         I020 = I010 + Nrow_ca*Ncol_km
         Iend = I020 + Nrow_ca*Ncol_km
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ab",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_km,2,Irrep_km,List_g)
C G(cb,km) -> G(km,cb)
         Call Transp(Work(I020),Work(I010),Ncol_km,Nrow_ca)

C V(km,cb)(t)*G(km,ca) -> S(b,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_a = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_a,Irrep_ca)
            Irrep_b  = Dirprd(Irrep_c,Irrep_cb)

            Na = Vrt(irrep_a,1)
            Nb = Vrt(irrep_b,1)
            Nc = Vrt(irrep_c,1)

            Nsum = Ncol_km*Nc
            Nrow = Nb
            Ncol = Na

            Icheck =  Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hvv_qp(Ioff),Nrow)
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

