      Subroutine Sym_expand(out,In,Ndim)
 
      Implicit None

      Integer*8 I
      Integer*8 J
      Integer*8 Ndim
     
      Double Precision In(Ndim)
      Double Precision Out(Ndim)

      Call Dcopy(Ndim,In,1,out,1)

      Return
      End
