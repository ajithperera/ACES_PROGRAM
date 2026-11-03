      Subroutine Putlstpy(Source,Column_Index,No_columns,Ndim,Irrep,
     +                    Left_index,Right_index)

C Write No_columns columns, starting at Column_Index, of MO list
C MOIO(Left_index,Right_index) from Source; the counterpart of Getlstpy,
C with the same meaning of the arguments (Irrep is not used).
C Nothing is written if the list does not exist, the column range is
C outside it or Source (Ndim long) is too short.

      Implicit None 

      Integer*8 Column_Index
      Integer*8 No_columns
      Integer*8 Ndim
      Integer*8 Irrep
      Integer*8 Left_index
      Integer*8 Right_index
      Integer*8 Obsolete
      Integer*8 Nrows
      Integer*8 Ncols
      Double Precision Source(Ndim) 

      Call Getlstdimpy(Left_index,Right_index,Nrows,Ncols)
      If (Nrows .Lt. 1 .Or. Ncols .Lt. 1) Then
         Write(6,"(a,i4,a,i4,a)") " @Putlstpy: list (",Left_index,",",
     +                            Right_index,") does not exist"
         Return
      Endif
      If (Column_Index .Lt. 1 .Or. No_columns .Lt. 1 .Or.
     +    Column_Index+No_columns-1 .Gt. Ncols) Then
         Write(6,"(a,i8,a,i8,a,i8)") " @Putlstpy: columns ",
     +         Column_Index," to ",Column_Index+No_columns-1,
     +         " requested, the list has ",Ncols
         Return
      Endif
      If (Ndim .Lt. Nrows*No_columns) Then
         Write(6,"(a,i12,a,i12)") " @Putlstpy: Ndim is ",Ndim,
     +         ", the columns need ",Nrows*No_columns
         Return
      Endif

      Obsolete = 1 
      Call A2_putlst(Source,Column_Index,No_columns,Obsolete,Left_index,
     +               Right_index)

      Return
      End 
