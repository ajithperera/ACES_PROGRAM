










      Subroutine Cgetrecpy(Name,Length,Dest)

      Implicit None

      Integer*8 Flag
      Integer*8 I
      Integer*8 Length
      Character*6 File
      Character*8 Name
      Character*8 Dest

      Write(6,"(a,a)")  "@-Entry to Getcrec,Rec name:", Name
      Write(6,"(a,I5)") "@-Entry to Getcrec,Rec len:", Length

      Flag = 20
      File = "JOBARC"
      Call A2_getcrec(Flag,File,Name,Length,Dest)

      Return
      End

      
