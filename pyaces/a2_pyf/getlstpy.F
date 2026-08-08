      Subroutine Getlstpy(Column_Index,No_columns,Irrep,Left_index,
     +                    Right_index,Destn,Ndim)

C This assumes that the elements are stored in columns grouped together 
C according to irreducibe
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
      Double Precision Destn(*)

      Obsolete = 1 
      Call A2_getlst(Destn,Column_Index,No_columns,Obsolete,Left_index,
     +               Right_index)

C      write(6,*)
C      write(6,"(a)") " The destination array @-getlstpy"
C      write(6,"(6(1x,F10.5))") (Destn(i),i=1,Ndim) 
      Return
      End
