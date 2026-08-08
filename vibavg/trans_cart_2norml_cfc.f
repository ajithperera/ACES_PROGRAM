










      subroutine Trans_Cart_2Norml_cfc(Cfc_cart,Cfc,Nrmlmodes,Work,
     +                                 Maxcor,Ndim,Nmodes,Iamlinear)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Nrmlmodes
      Dimension Cfc_cart(Ndim,Ndim)
      Dimension Cfc(Nmodes,Nmodes)
      Dimension Nrmlmodes(Ndim,Ndim)
      Dimension Work(Maxcor)

      Data Done,Dnull,Ione/1.0D0,0.0D0,1/

      Icol = 7
      If (Iamlinear .Gt. 0) Icol = 6
      I000 = Ione
      Iend = I000 + Ndim*Nmodes 
      If (Iend .Ge. Maxcor) Call Insmem("trans_Cart_2Norml_cfc",
     +                                   Iend,Maxcor)

      Call Dgemm("T","N",Nmodes,Ndim,Ndim,Done,Nrmlmodes(1,Icol),
     +            Ndim,Cfc_cart,Ndim,Dnull,Work(I000),Nmodes)
      Call Dgemm("N","N",Nmodes,Nmodes,Ndim,Done,Work(I000),Nmodes,
     +            Nrmlmodes(1,Icol),Ndim,Dnull,Cfc,
     +            Nmodes)

      Write(6,"(a)") " Transformed cubic force constants"
      call output(Cfc,1,Nmodes,1,Nmodes,Nmodes,Nmodes,1)

      Return
      End 

