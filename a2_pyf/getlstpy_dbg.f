      Subroutine Getlstpy_dbg(Column_Index,No_columns,Irrep,Left_index,
     +                        Right_index,Ndim)

C This assumes that the elements are stored in columns grouped together 
C according to irreducible representations 
C Destn: Destination array of Ndim size,
C Column_index : Index of the column that needs to retrived.
C No_columns : Number of columns retrived.
C Ndim : Length of the destination vector.
C Left_index  : Irreducible represention of the column
C Right_indes : Index unique to the data type.

      Implicit None 

      Integer*8 Column_Index
      Integer*8 No_columns
      Integer*8 Ndim
      Integer*8 I
      Integer*8 Irrep
      Integer*8 Left_index
      Integer*8 Right_index
      Integer*8 Obsolete
      Double Precision Destn(Ndim) 

      Obsolete=1
      Print*, Column_Index,No_columns,Ndim,Left_index,Right_index,
     +        Obsolete
C      Call A2_getlst(Destn,Column_Index,No_columns,Obsolete,Left_index,
C     +               Right_index)

      Call Getlst(Destn,1,1,1,1,91)
      Write(6,"(a)") "The destination @-Getlstpy_dbg"
      Write(6,"(6(1x,F10.6))") (Destn(I), I=1,Ndim)

      Return
      End 
