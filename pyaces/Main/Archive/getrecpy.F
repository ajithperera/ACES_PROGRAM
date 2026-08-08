      Subroutine Getrecpy(Name,Length,Dest)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Character*6 File
      Character*8 Name
      Double Precision Dest(Length)

      include "machsp.com"
      include "jobarc.com"

      Flag = 20
      File = "JOBARC"
      Call A2_getrec(Flag,File,Name,Length,Dest)

      Return
      End

      
