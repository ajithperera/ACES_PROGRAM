










      Subroutine Save_2alt_list(Work,Maxcor,T1ln,T2ln,T1ln_aa,
     +                          T1ln_bb,T2ln_aa,T2ln_bb,T2ln_ab,
     +                          Irrepx,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Integer T1ln,T2ln
      Integer T1ln_aa,T1ln_bb,T2ln_aa,T2ln_bb,T2ln_ab

      Dimension Work(Maxcor)
 
      Data Ione /1/

      List1_h = 490
      List1_t = 493

      I000 = Ione
      Call Getlst(Work(I000),1,1,1,1,List1_h)
      Call Putlst(Work(I000),1,1,1,1,List1_t)

      If (Iuhf .Ne. 0) Then
         Call Getlst(Work(I000),1,1,1,2,List1_h)
         Call Putlst(Work(I000),1,1,1,2,List1_t)
      Endif

      List2_h = 444
      List2_t = 454

      Call Getall(Work(I000),T2ln_aa,Irrepx,List2_h)
      Call Putall(Work(I000),T2ln_aa,Irrepx,List2_t)

      If (Iuhf .Ne. 0) Then
         Call Getall(Work(I000),T2ln_bb,Irrepx,List2_h+1)
         Call Putall(Work(I000),T2ln_bb,Irrepx,List2_t+1)
      Endif 

      Call Getall(Work(I000),T2ln_ab,Irrepx,List2_h+2)
      Call Putall(Work(I000),T2ln_ab,Irrepx,List2_t+2)



      Return
      End
      

