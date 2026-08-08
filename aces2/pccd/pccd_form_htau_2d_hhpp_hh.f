










      Subroutine Pccd_form_htau_2d_hhpp_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                                      Nocc,Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hoo_pq(Nocc*Nocc)
      Dimension Hoo_qp(Nocc*Nocc)

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

C V(ce,kj)*Gamma(ce,ki) (Gamma stored as ce,ki)
   
      Irrepx = Ione

      Do Irrep_ki = 1, Nirrep
         Irrep_ce = Dirprd(Irrep_ki,Irrepx)
         Irrep_kj = Dirprd(Irrep_ce,Irrepx)

         Nrow_ce = Irpdpd(Irrep_ce,19)
         Ncol_kj = Irpdpd(Irrep_kj,21)
         Nrow_ce = Irpdpd(Irrep_ce,19)
         Ncol_ki = Irpdpd(Irrep_ki,21)

         I000 = Ione
         I010 = I000 + Nrow_ce*Ncol_kj
         I020 = I010 + Max(Ncol_kj,Nrow_ce)
         Iend = I020 + Max(Ncol_kj,Nrow_ce)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_hh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_kj,1,Irrep_kj,List_v)
         Call Spinad3(Irrep_ce,Vrt(1,1),Nrow_ce,Ncol_kj,Work(I000),
     +                Work(I010),Work(I020))
         Iend = I010 + Nrow_ce*Ncol_ki
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ij",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_ki,2,Irrep_ki,List_g)

C V(ce,kj)(t)*G(ce,ki) (stored as (ce,ki)) -> S(j,i)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_i = 1, Nirrep
            Irrep_k  = Dirprd(Irrep_i,Irrep_ki)
            Irrep_j  = Dirprd(Irrep_k,Irrep_kj)

            Ni = Pop(irrep_i,1)
            Nj = Pop(irrep_j,1)
            Nk = Pop(irrep_k,1)

            Nsum = Nrow_ce*Nk
            Nrow = Nj
            Ncol = Ni

            Icheck =  Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hoo_qp(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Nsum*Nj
            Koff = Koff + Nsum*Ni
         Enddo 
      Enddo

      call pccd_check_htau("Htau_oo :",Hoo_qp,Nocc,"OO","D")

      Return
      End

