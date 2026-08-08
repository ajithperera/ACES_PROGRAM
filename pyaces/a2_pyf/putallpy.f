










      Subroutine Putallpy(Source,Ndim,Irrep,Right_index)

C This assumes that the elements are stored in columns grouped together 
C according to irreducibe
C Destn: Destination array of Ndim size,
C No_columns : Number of columns retrived.
C Ndim : Length of the destination vector.
C Irrep  : Irreducible represention of the data
C Right_indes : Index unique to the data type

      Implicit None 

      Integer*8 Ndim
      Integer*8 Irrep
      Integer*8 Right_index
      Double Precision Source(Ndim) 

      Call A2_putall(Source,Ndim,Irrep,Right_index)

      Return
      End 
