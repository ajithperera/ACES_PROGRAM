










      SUbroutine Tdee_return_mubar_dot(Work,Memleft,Iuhf,Irrepx,Iside)

      Implicit Double Precision(A-H,O-Z)
      Logical Source, Target

      Dimension Work(Memleft)

      Source = .False.
      Target = .True. 
      Call Tdee_init_lists(Work,Memleft,Iuhf,Irrepx,Source,
     +                     Target,Iside)

      Call Tdee_hbar_mult(Work,Memleft,Iuhf,Irrepx,Iside)


      Return
      End
