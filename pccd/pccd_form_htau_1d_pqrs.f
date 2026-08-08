










      Subroutine Pccd_form_htau_1d_pqrs(Htau_pq,Htau_qp,Hpq,Dpq,Hoo_pq,
     +                                  Hoo_qp,Hvv_pq,Hvv_qp,Hvo,Hov,
     +                                  
     +                                  Htau_oo,Hoo,Doo,Htau_vv,Hvv,Dvv,
     +                                  Htau_ov,Htau_vo,H_ov,H_vo,D_ov,
     +                                  D_vo,Work,Nocc,Nvrt,Maxcor,
     +                                  Nbas)

      Implicit Double Precision(A-H,O-Z)

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Hpq(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Hoo_pq(Nocc,Nocc)
      Dimension Hoo_qp(Nocc,Nocc)
      Dimension Hvv_pq(Nvrt,Nvrt)
      Dimension Hvv_qp(Nvrt,Nvrt)
      Dimension Hvo(Nvrt*Nocc)
      Dimension Hov(Nocc*Nvrt)
      Dimension Work(Maxcor)

      Dimension Htau_oo(Nocc,Nocc)
      Dimension H_oo(Nocc,Nocc)
      Dimension D_oo(Nocc,Nocc)
      Dimension Htau_vv(Nvrt,Nvrt)
      Dimension H_vv(Nvrt,Nvrt)
      Dimension D_vv(Nvrt,Nvrt)
      Dimension Htau_ov(Nocc,Nvrt)
      Dimension Htau_vo(Nvrt,Nocc)
      Dimension H_ov(Nocc,Nvrt)
      Dimension H_vo(Nvrt,Nocc)
      Dimension D_ov(Nocc,Nvrt)
      Dimension D_vo(Nvrt,Nocc)
      Dimension Ioffo(8)
      Dimension Ioffv(8)

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

      Data Ione /1/

C Symmetry pack the Htau_pq,Hpq and Dpq matrices. 

      Call Pccd_symm_pack(Htau_qp,Htau_oo,Htau_vv,Htau_ov,Htau_vo,
     +                  Nocc,Nvrt,Nbas)
      Call Pccd_symm_pack(Hpq,H_oo,H_vv,H_ov,H_vo,Nocc,Nvrt,Nbas)
      Call Pccd_symm_pack(Dpq,D_oo,D_vv,D_ov,D_vo,Nocc,Nvrt,Nbas)


C Form Htau_pq contributions Htau_pqrs (This is four terns like 
C Htau_pq*Delta(rs) and all possible permutations. 

      I000 = Ione
      I010 = I000 + Nocc*Nocc
      Iend = I010 + Nvrt*Nvrt
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs",
     +                                   Iend,Maxcor)

      Maxcor = Maxcor - Iend 
      Call Pccd_form_htau_1d_pqrs_1(Htau_qp,Htau_oo,Htau_vv,Htau_ov,
     +                              Htau_vo,Work(I000),Work(I010),
     +                              Work(Iend),Ioffo,Ioffv,
     +                              Maxcor,Nocc,Nvrt,Nbas)
      Call Pccd_form_htau_1d_pqrs_2(Hpq,Dpq,H_oo,H_vv,H_ov,H_vo,D_oo,
     +                              Dvv,D_ov,D_vo,Work(Iend),
     +                              Ioffo,Ioffv,Maxcor,Nocc,Nvrt,Nbas)

      Return
      End 
   
   
