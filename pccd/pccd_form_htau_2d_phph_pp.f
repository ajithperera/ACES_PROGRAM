










      Subroutine Pccd_form_htau_2d_phph_pp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                                     Nvrt,Nbas,List_v,List_g) 
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

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C V(mb,ne)*Gamma(ma,ne)
C V(bm,en)*G(am,en)=V(en,bm)G(en,am)

      Irrepx = Ione
      Ispin  = Ione
  
      Do Irrep_am = 1, Nirrep
         Irrep_en = Dirprd(Irrep_am,Irrepx)
         Irrep_bm = Dirprd(Irrep_en,Irrepx)

         Nrow_en = Irpdpd(Irrep_en,11)
         Ncol_bm = Irpdpd(Irrep_bm,11)
         Ncol_am = Irpdpd(Irrep_am,11)

         I000 = Ione
         I010 = I000 + Nrow_en*Ncol_bm
         I020 = I010 + Nrow_en*Ncol_am
         I030 = I020 + Nrow_en
         I040 = I030 + Nrow_en
         Iend = I040 + Nrow_en

         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_vv",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_bm,1,Irrep_bm,List_v)
         Call Getlst(Work(I010),1,Ncol_am,2,Irrep_am,List_g)

         Call Symtr1(Irrep_bm,Vrt(1,Ispin),Pop(1,Ispin),Nrow_en,
     +               Work(I000),Work(I020),Work(I030),Work(I040))
         Call Symtr1(Irrep_am,Vrt(1,Ispin),Pop(1,Ispin),Nrow_en,
     +               Work(I010),Work(I020),Work(I030),Work(I040))

C V(en,mb)^(t)*G(en,ma) -> S(b,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_a = 1, Nirrep
            Irrep_m  = Dirprd(Irrep_a,Irrep_am)
            Irrep_b  = Dirprd(Irrep_m,Irrep_bm)

            Na = Vrt(irrep_a,1)
            Nb = Vrt(irrep_b,1)
            Nm = Pop(irrep_m,1)
            
            Nsum = Nrow_en*Nm
            Nrow = Nb
            Ncol = Na
            Icheck = Min(Nsum,Nrow,Ncol)
            
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hvv_qp(Ioff),Nrow)
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

