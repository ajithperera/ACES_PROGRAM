










      Subroutine Pccd_form_htau_2d_phph_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                                     Nocc,Nbas,List_v,List_g) 
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

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C V(ek,cj)*Gamma(ek,ci) 

      Irrepx = Ione
  
      Do Irrep_ek = 1, Nirrep
         Irrep_ci = Dirprd(Irrep_ek,Irrepx)
         Irrep_cj = Dirprd(Irrep_ek,Irrepx)

         Nrow_ek = Irpdpd(Irrep_ek,9)
         Ncol_cj = Irpdpd(Irrep_cj,9)
         Nrow_ci = Irpdpd(Irrep_ci,9)
         Ncol_ek = Irpdpd(Irrep_ek,9)

         I000 = Ione
         I010 = I000 + Nrow_ek*Ncol_cj 
         Iend = I010 + Nrow_ci*Ncol_ek
  
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_hh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_cj,1,Irrep_cj,List_v)
         Call Getlst(Work(I010),1,Ncol_ek,2,Irrep_ek,List_g)

C V(ek,cj)^(t)*G(ek,ci) -> S(j,i)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_i = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_i,Irrep_ci)
            Irrep_j  = Dirprd(Irrep_c,Irrep_cj)

            Nj = Pop(irrep_j,1)
            Ni = Pop(irrep_i,1)
            Nc = Vrt(irrep_c,1)
            
            Nsum = Nrow_ek*Nc
            Nrow = Nj
            Ncol = Ni
            Icheck = Min(Nsum,Nrow,Ncol)
            
            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hoo_qp(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Nsum*Nj
            Koff = Koff + Nsum*Ni
         Enddo
      Enddo

      call pccd_check_htau("Htau_oo :",Hoo_qp,Nocc,"OO","D")

      Return
      End

