










      Subroutine get_ccoefs_block(Ccoefs,Ccoef_block,Ishell,
     +                            Npcoef,Ccoef_off,Length,Nshells)

      Implicit Double Precision (A-H, O-Z)

      Dimension Ccoefs(Npcoef)
      Dimension Ccoef_block(Length)
      Integer Ccoef_off(Nshells)

      ioff = Ccoef_off(Ishell) 

      Call Dcopy(Length,Ccoefs(ioff),1,Ccoef_block,1)
      

      
      Return
      End 
     
