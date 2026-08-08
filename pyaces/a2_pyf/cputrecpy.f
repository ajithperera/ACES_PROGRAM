










      Subroutine Cputrecpy(Name,Length,Source)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Character*6 File
      Character*8 Name
      Character*8 Source

Cf2py intent (inout) Source

      Flag = 20
      File = "JOBARC"
      Call A2_putcrec(Name,Length,Source)
     
      Return
      End

      
