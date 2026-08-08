










      Subroutine Tdee_init_global_Vars(Irrepx,Iuhf,Iside,Nsize,
     +                                 Memleft)

      Implicit Double Precision(A-H,O-Z)


      Common /Tdee_vars/Nsize_dummy,Irrepx_dummy,Iside_dummy,
     +                  Iuhf_dummy,Memleft_Dummy



      Write(6,*)
      Write(6,"(a)") "  Initilizing global variables" 
     
      Irrepx_dummy  = Irrepx
      Iuhf_dummy    = Iuhf
      Nsize_dummy   = Nsize
      Iside_dummy   = Iside
      Memleft_dummy = Memleft

      Return
      End
