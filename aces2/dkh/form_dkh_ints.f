










      Subroutine form_dkh_ints(Dkh_block,Dkh_int,Ishell,Jshell,
     +                         Icntrfns,Jcntrfns,Ncnfns,
     +                         Prim_off,Cont_off,Nshells)

      Implicit Double Precision (A-H, O-Z)

      Dimension Dkh_block(Icntrfns,Jcntrfns)
      Dimension Dkh_int(Ncnfns,Ncnfns)
      Integer Prim_off(Nshells),Cont_off(Nshells)

      ioff = Cont_off(Ishell) - 1
      Joff = Cont_off(Jshell) - 1


      Jct = 0
      Do J = Joff+1, Jcntrfns + Joff
         Jct = Jct + 1
         Ict = 0
         Do I= Ioff+1,Icntrfns + Ioff
            Ict = Ict + 1
            Dkh_int(I,J) = Dkh_block(Ict,Jct)
         Enddo
      Enddo

      
      Return
      End 
     
