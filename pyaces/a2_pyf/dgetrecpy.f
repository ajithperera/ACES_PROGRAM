










      Subroutine Dgetrecpy(Name,Length,Dest)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Character*6 File
      Character*8 Name
      Double Precision Dest(Length)

      Flag = 20
      File = "JOBARC"
CSSS      Write(6,"(a,i4)") " The requested record length", Length 
      Call A2_getrec(Flag,File,Name,Length,Dest)

      Return
      End

      
