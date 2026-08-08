










      Subroutine Psi4dbg_form_htau_2d(Htau_pq,Htau_qp,Hoo_pq,Hoo_qp,
     +                                Hvv_pq,Hvv_qp,Hvo,Hov,Dpq,Dhf,
     +                                Dcc,Work,Nocc,Nvrt,Maxcor,Nbas,
     +                                E)
   
      Implicit Double Precision(A-H,O-Z)
      Character*4 Spin
      Logical pCCD,CCD,LCCD

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
      Dimension Dcc(Nbas,Nbas)
      Dimension Hoo_pq(Nocc,Nocc)
      Dimension Hoo_qp(Nocc,Nocc)
      Dimension Hvv_pq(Nvrt,Nvrt)
      Dimension Hvv_qp(Nvrt,Nvrt)
      Dimension Hvo(Nvrt*Nocc)
      Dimension Hov(Nocc*Nvrt)
      Dimension Work(Maxcor)
    
      Common /Meth/pCCD,CCD,LCCD

      Data One,Onem,Dnull,Half,Quart,Two /1.0D0,-1.0D0,0.0D0,0.50D0,
     +                                0.25D0,2.0D0/

C The OO block

      List_v = 16
      List_g = 116
      Call Psi4dbg_form_htau_2d_hhpp_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                                  Nocc,Nbas,List_v,List_g)

C The VV block

      List_v = 16
      List_g = 116
      Call Psi4dbg_form_htau_2d_hhpp_pp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                                  Nvrt,Nbas,List_v,List_g)


C The OV and VO blocks


      List_v1 = 30
      List_v2 = 10
      List_g  = 116
      Call Pccd_form_htau_2d_hhpp_ph(Hov,Hvo,Work,Maxcor,Nocc,
     +                               Nvrt,Nbas,List_v1,List_v2,
     +                               List_g)
      List_v = 10
      List_g = 113

      Call Pccd_form_htau_2d_hhhh_ph(Hvo,Work,Maxcor,Nocc,
     +                               Nvrt,Nbas,List_v,List_g)

      List_v  = 30
      List_g  = 133

      Call Pccd_form_htau_2d_pppp_ph(Hov,Work,Maxcor,Nocc,
     +                               Nvrt,Nbas,List_v,List_g)

      Call Pccd_sortgam(Work,Maxcor,0)

      List_v1 = 30
      List_v2 = 7
      List_g  = 123
      Fact    = One

      Call Pccd_form_htau_2d_phph_ph_0(Hov,Hvo,Work,Maxcor,Nocc,
     +                                 Nvrt,Nbas,List_v1,List_v2,
     +                                 List_g,Fact)

      List_v1 = 30
      List_v2 = 10
      List_g  = 125
      Fact    = One

      Call Pccd_form_htau_2d_phph_ph_1(Hov,Hvo,Work,Maxcor,Nocc,
     +                                 Nvrt,Nbas,List_v1,List_v2,
     +                                 List_g,Fact)

      Call Psi4dbg_form_htau(Htau_qp,Htau_pq,Hoo_qp,Hoo_pq,Hvv_qp,
     +                       Hvv_pq,Hov,Hvo,Work,Maxcor,Nocc,Nvrt,
     +                       Nbas)
      Return
      End 
