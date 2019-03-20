C----------------------------------------------------------------
C This include file combined with the SIP set_aces3_vars can be
C used to set variables for fractional occupation CC logic.
C 
      Integer O,V,OO(2),VV(2),OOO(3),VVV(3),OOOO(4),VVVV(4)
      Integer Contract,Maxbfns 
      Double Precision Occ_numsa,Occ_numsb
      Parameter (Maxbfns = 2000)
      common /Permutations/ O,V,OO,VV,OOO,VVV,OOOO,VVVV
      common /Frac_occ_dims/ Contract 
      Common /Frac_occs/ Occ_numsa(Maxbfns),Occ_numsb(Maxbfns)

C----------------------------------------------------------------


