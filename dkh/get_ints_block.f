










      Subroutine get_ints_block(Coreham,Shell_block,Ishell,Jshell,
     +                           Nprims,Prim_off,Cont_off,
     +                           Iprimfns,jprimfns,Nshells)

      Implicit Double Precision (A-H, O-Z)

      Dimension Coreham(Nprims, Nprims)
      Dimension Shell_block(iprimfns,jprimfns)
      Integer Prim_off(Nshells),Cont_off(Nshells)

      ioff = Prim_off(Ishell) - 1
      Joff = Prim_off(Jshell) - 1

      Do J = 1, Jprimfns
         Do I= 1, Iprimfns
            shell_block(I,J) = Coreham(Ioff+I,Joff+J)
         Enddo 
      Enddo 

      
      Return
      End 
     
