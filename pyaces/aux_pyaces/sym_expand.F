      Subroutine Sym_expand(outarr,In,Ndim)

      Implicit None

      Integer*8 I
      Integer*8 J
      Integer*8 Ndim

      Double Precision In(Ndim)
      Double Precision outarr(Ndim)

      Call Dcopy(Ndim,In,1,outarr,1)

      Return
      End
