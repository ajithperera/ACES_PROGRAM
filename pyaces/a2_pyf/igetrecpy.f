










      Subroutine IGetrecpy(Name,Length,IDest)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Integer*8 IDest(Length)
      Character*6 File
      Character*8 Name

      Flag = 20
      File = "JOBARC"
      Call A2_getrec(Flag,File,Name,Length,IDest)
     
      Return
      End

      
