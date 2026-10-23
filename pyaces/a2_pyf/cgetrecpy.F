      Subroutine Cgetrecpy(Name,Length,Dest)

      Implicit None

      Integer*8 Flag
      Integer*8 Length
      Character*6 File
      Character*8 Name
      Character*1024 Dest
Cf2py intent(out) Dest

      Dest = " "
      If (Length .Gt. Len(Dest)) Then
         Write(6,"(a,a,a,I8,a,I5)") "@CGETRECPY: record ", Name,
     &      " length ", Length, " exceeds ", Len(Dest)
         Call aces_exit(1)
      Endif

      Flag = 20
      File = "JOBARC"
      Call A2_getcrec(Flag,File,Name,Length,Dest)

      Return
      End

      
