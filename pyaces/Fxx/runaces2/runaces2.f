










      Subroutine runaces2(Flag)

      Implicit None
      Integer*8 Flag

      If (Flag .Eq. 1) Then
         write(6,"(a)") "@-runaces2: Calling main"
         Call main()
      Else
         write(6,"(a)") "@-runaces2: Returing empty handed"
      Endif

      Return
      End


