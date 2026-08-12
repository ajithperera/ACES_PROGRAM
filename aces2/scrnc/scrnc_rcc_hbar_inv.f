










      Subroutine Scrnc_rcc_hbar_inv(work,Memleft,Iuhf,Imult,Iloc)
  
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft) 
      Dimension Iloc(8)

      Write(6,"(1x,a)") "-----Entering scnc_rcc_hbar_inv------"
      Call Scrnc_form_rcc_hbar(Work,Memleft,Iuhf,Imult)
      Call Scrnc_form_rcc_hbar_inv(Work,Memleft,Iuhf,Iloc) 
      
      Return
      End

