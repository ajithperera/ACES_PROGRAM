      Subroutine Dputrecpy(Name,Length,Source)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Character*6 File
      Character*8 Name
      Double Precision Source(Length)

      Flag = 20
      File = "JOBARC"
      Call A2_putrec(Flag,File,Name,Length,Source)
     
      Return
      End

      
