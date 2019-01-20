C----------------------------------------------------------------
C This include file contains data relevent to truncation of the 
C subspace. There several SIP that access this header file.
C
      Integer Ndim_old,Ndim_new,Maxdim,Iselect
      Double Precision Omega
      Parameter (Maxdim=100)
      Double Precision Eval_selct(Maxdim),Evec_selct(Maxdim,Maxdim)
      common /Truncate/Ndim_old,Ndim_new,Omega 
      common /eomdata/ Iselect

C----------------------------------------------------------------



