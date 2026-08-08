










      Subroutine Pccd_form_htau_2d_hhhh(Hoo_pq,Hoo_qp,Work,Maxcor,Nocc,
     +                                  Nbas,List_v,List_g) 
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

      Data One,Two,Onem,Dnull,Ione /1.0D0,2.0D0,-1.0D0,0.0D0,1/

C V(mk,nj)*Gamma(mk,ni)
C Note that this is commented because of the OO gradient term added in 
C Pccd_form_htau_pq (Pccd_form_htau_oo). 

      Call Dzero(Hoo_qp,Nocc*Nocc)
      Call Dzero(Hoo_pq,Nocc*Nocc)

      Irrepx = Ione
      Do Irrep_ni = 1, Nirrep
         Irrep_mk = Dirprd(Irrep_ni,Irrepx)
         Irrep_nj = Dirprd(Irrep_mk,Irrepx)

         Nrow_mk = Irpdpd(Irrep_mk,21)
         Ncol_nj = Irpdpd(Irrep_nj,21)
         Nrow_ni = Irpdpd(Irrep_ni,21)
         Ncol_mk = Irpdpd(Irrep_mk,21)

         I000 = Ione
         I010 = I000 + Nrow_mk*Ncol_nj
         I020 = I010 + Max(Nrow_mk,Ncol_nj)
         Iend = I020 + Max(Nrow_mk,Ncol_nj)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_nj,1,Irrep_nj,List_v)
         Call Spinad1(Irrep_nj,Pop(1,1),Nrow_mk,Work(I000),
     +                Work(I010),Work(I020))

         Iend = I010 + Nrow_ni*Ncol_mk
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh",
     +                                      Iend,Maxcor)
         
         Call Getlst(Work(I010),1,Ncol_mk,2,Irrep_mk,List_g)

C V(mkn,j)^(t)*G(mkn,i) -> S(j,i)

         Ioff = Ione
         Joff = I000
         Koff = I010

         Do irrep_i = 1, Nirrep
            Irrep_n  = Dirprd(Irrep_i,Irrep_ni)
            Irrep_j  = Dirprd(Irrep_n,Irrep_nj)

            Nj = Pop(irrep_j,1)
            Ni = Pop(irrep_i,1)
            Nn = Pop(irrep_n,1)

            Nsum = Nrow_mk*Nn
            Ncol = Ni
            Nrow = Nj

            Icheck = Min(Nsum,Ni,Nj)
            If (Icheck .Ne. 0)  then
               Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                   Nsum,Work(Koff),Nsum,One,Hoo_qp(Ioff),
     +                   Nrow)
            Endif 

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Nrow_mk*Nn*Nj
            Koff = Koff + Ncol_mk*Nn*Ni
         Enddo
      Enddo

      Write(6,*)
      call pccd_check_htau("Htau_oo :",Hoo_qp,Nocc,"OO","D")
   
CSSS      call output(Hoo_qp,1,Nocc,1,Nocc,Nocc,Nocc,1)

      Return
      End
