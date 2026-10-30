      Subroutine Cputrecpy(Name,Length,Source)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Character*6 File
      Character*8 Name
C Source holds Length characters (assumed length: Character*8 here passed
C only the first 8 and stored whatever followed them in memory).
      Character*(*) Source

      Flag = 20
      File = "JOBARC"
      Call A2_putcrec(Flag,File,Name,Length,Source)
     
      Return
      End

      
