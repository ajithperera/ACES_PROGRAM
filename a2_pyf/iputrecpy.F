      Subroutine Iputrecpy(Name,Length,Isource)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Character*6 File
      Character*8 Name
      Integer*8 Isource(Length)

      Flag = 20
      File = "JOBARC"
      Call A2_putrec(Flag,File,Name,Length,Isource)
     
      Return
      End

      
