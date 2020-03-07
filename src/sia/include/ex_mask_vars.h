C----------------------------------------------------------------
C This include file is used to built a excitation mask. This is
C used in search for core-states. The elements of the ibelong array
C if they are interest to represent core excitations.

      Integer Max_nocc,Ibelong_a,Ibelong_b,Max_pdim,Nsites 
      Double precision P_ovlp,Core_thres
      Parameter (Max_nocc =500)
      Parameter (Max_pdim =100)
      Common /Ee_mask/Ibelong_a(Max_nocc),Ibelong_b(Max_nocc)
      Common /P_ovlp/P_ovlp(Max_pdim),Nsites 

C----------------------------------------------------------------




