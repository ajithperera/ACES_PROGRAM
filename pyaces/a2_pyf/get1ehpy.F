      Subroutine Get1ehpy(Name,Dest,Ndim)

      Implicit None
      Integer*8 Flag 
      Integer*8 Length
      Integer*8 Memneed
      Integer*8 Nao
      Integer*8 J0
      Integer*8 J1
      Integer*8 J2
      Integer*8 I
      Integer*8 J
      Integer*8 K
      Integer*8 N
      Integer*8 Ndim
      Integer*8 Nirrep
      Integer*8 Ione
      Integer*8 Maxfns 

      Integer*8 Jcore(1)
      Integer*8 Nbfirr(8)
      Integer*8 Trioff(8)
      Integer*8 Bfoff(8)

      Double Precision Dest(*)

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
     +                           " does not match with the internally",
     +                           " required:", Nao
         Return
      Endif

      If (Name(1:7) .Eq. "overlap") Then
         Label = 'OVERLAP '
      Elseif (Name(1:4) .Eq. "oneh") Then
         Label = 'ONEHAMIL'
      Elseif (Name(1:7) .Eq. "kinetic") Then
         Label = 'KINETINT'
      Else
         Write(6,"(3a)") " @Get1ehpy: unknown integral name '",
     +                   Name(1:Len_trim(Name)),
     +                   "' (overlap, oneh or kinetic)"
         Return
      Endif 
  
      Call A2_getrec(Flag,File,"NIRREP  ",Ione,Nirrep)
      Call A2_getrec(Flag,File,"NUMBASIR",NIrrep,Nbfirr)
      Maxfns = 0 
      Do I = 1,Nirrep
         Maxfns = Max(Maxfns,Nbfirr(I)) 
      Enddo 

C Aces_malloc leaves Jcore(1) = 0 on failure. J0 can be negative (the
C heap lies below Jcore), so it says nothing about success.

      Jcore(1) = 0
      Memneed = Nao*(Nao+1)/2 + Maxfns*Maxfns
      Call Aces_malloc(Memneed,Jcore,J0)
  
      If (Jcore(1) .Eq. 0) Then
          Write(6,"(a,ES12.6E2,a)") " Unable to allocate", 
     +                Dble(Memneed)," integers of memory"
          Return
      Endif 

      Call Get1Eints(Label,Jcore(J0),Nao*(Nao+1)/2)

      If (Nirrep .Eq. Ione) Then

         Call Expnd2(Jcore(J0),Dest,Nao)
         Call C_free(Jcore)

      Else

C The IIII file holds the lower triangle of each irrep block, one block
C after the other. Dest is the full Nao x Nao matrix in the SO basis
C (functions in irrep order), block diagonal: block K occupies rows and
C columns Bfoff(K)+1 .. Bfoff(K)+Nbfirr(K).

         Trioff(1) = 0
         Bfoff(1)  = 0
         Do J = 1, Nirrep-1
            Trioff(J+1) = Trioff(J) + Nbfirr(J)*(Nbfirr(J)+1)/2
            Bfoff(J+1)  = Bfoff(J)  + Nbfirr(J)
         Enddo 
            
         J1 = J0
         J2 = J1 + Nao*(Nao+1)/2

         Call Zero(Dest,Nao*Nao)
         Do K = 1, Nirrep
            N = Nbfirr(K)
            If (N .Gt. 0) Then
               Call Expnd2(Jcore(J1+Trioff(K)),Jcore(J2),N)
               Do J = 1, N
                  Call Dcopy(N,Jcore(J2+(J-1)*N),1,
     +                       Dest(Bfoff(K)+1+(Bfoff(K)+J-1)*Nao),1)
               Enddo
            Endif
         Enddo
         Call C_free(Jcore)
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " @-Fortran-get1ehpy: The retrieved data"
      Write(6,"(a)") " Returned matrix (all irreps)"
      call output(Dest,1,nao,1,nao,nao,nao,1)
#endif 

      Return
      End
