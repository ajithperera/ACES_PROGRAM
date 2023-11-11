










      Subroutine Get2ehpy(Name,Dest,Ndim)

      Implicit None

      Integer*8 Flag 
      Integer*8 Length
      Integer*8 Ione
      Integer*8 Ndim
      Integer*8 Nao
      Integer*8 I
  
      Double Precision Dest(Ndim*Ndim*Ndim*Ndim)

      Character*80 Name
      Character*6 File
      Character*8 Label 

      Data Ione /1/

      Flag = 0
      File = "JOBARC"
      Call A2_getrec(Flag,File,"NBASTOT",Length,I)
      If (Length .Gt. 0) Then
         Flag = 1
         Call A2_getrec(Flag,File,"NBASTOT",Length,Nao)
      Else
         Write(6,"(a,a)")" The number of basis functions have not been",
     +                   " recorded. Most likely initilization"
         Write(6,"(a)")  " step has not been done."
         Return
      Endif 

      If (Ndim .Ne. Nao) Then
         Write(6,"(a,i3,2a,i3)") " The assigned dimension: ", Ndim, 
     +                           " does not match with the internally"
          Write(6,"(a,i3)")      " required:", Nao
         Return
      Endif

      If (Name(1:7) .Eq. "2elints") Then
         Label = 'TWOELSUP'
      Endif 

      Call Get2Eints(Label,Dest,Nao)

      Return
      End
 

    

