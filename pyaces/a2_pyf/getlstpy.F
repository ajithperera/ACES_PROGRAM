      Subroutine Getlstpy(Column_Index,No_columns,Irrep,Left_index,
     +                    Right_index,Destn,Ndim)

C Read No_columns columns, starting at Column_Index, of MO list
C MOIO(Left_index,Right_index) into Destn, as getlst(Z,first,n,1,irrep,
C list) does in ACES II.
C Left_index  : irrep of the columns (the MOIO row index)
C Right_index : list number (the MOIO column index, see pyaces.moio)
C Ndim        : length of Destn; must hold No_columns columns
C Irrep       : not used; the irrep is Left_index. Kept so that existing
C               callers keep working.
C Nothing is read, and Destn is left as it was, if the list does not
C exist, the column range is outside it or Destn is too short.

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
      Double Precision Destn(*)

      Call Getlstdimpy(Left_index,Right_index,Nrows,Ncols)
      If (Nrows .Lt. 1 .Or. Ncols .Lt. 1) Then
         Write(6,"(a,i4,a,i4,a)") " @Getlstpy: list (",Left_index,",",
     +                            Right_index,") does not exist"
         Return
      Endif
      If (Column_Index .Lt. 1 .Or. No_columns .Lt. 1 .Or.
     +    Column_Index+No_columns-1 .Gt. Ncols) Then
         Write(6,"(a,i8,a,i8,a,i8)") " @Getlstpy: columns ",
     +         Column_Index," to ",Column_Index+No_columns-1,
     +         " requested, the list has ",Ncols
         Return
      Endif
      If (Ndim .Lt. Nrows*No_columns) Then
         Write(6,"(a,i12,a,i12)") " @Getlstpy: Ndim is ",Ndim,
     +         ", the columns need ",Nrows*No_columns
         Return
      Endif

      Obsolete = 1 
      Call A2_getlst(Destn,Column_Index,No_columns,Obsolete,Left_index,
     +               Right_index)

      Return
      End
