










      Subroutine Putlstpy(Source,Column_Index,No_columns,Ndim,Irrep,
     +                    Left_index,Right_index)

C This assumes that the elements are stored in columns grouped together 
C according to irreducibe
C Source: Source array of Ndim size,
C Column_index : Index of the column that needs to retrived.
C No_columns : Number of columns retrived.
C Ndim : Length of the destination vector.
C Left_index  : Irreducible represention of the column
C Right_indes : Index unique to the data type.

      Implicit None 

      Integer*8 Column_Index
      Integer*8 No_columns
      Integer*8 Ndim
      Integer*8 Irrep
      Integer*8 Left_index
      Integer*8 Right_index
      Integer*8 Obsolete
      Double Precision Source(Ndim) 

      Obsolete = 1 
      Call A2_putlst(Source,Column_Index,No_columns,Obsolete,Left_index,
     +               Right_index)

      Return
      End 
