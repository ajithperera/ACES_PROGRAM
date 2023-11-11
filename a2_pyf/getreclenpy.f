










      Subroutine Getreclenpy(Name,Length)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Double Precision Dest
      Character*6 File
      Character*8 Name

      Flag = 0
      File = "JOBARC"
      Call A2_getrec(Flag,File,Name,Length,Dest)

      Return
      End

      
