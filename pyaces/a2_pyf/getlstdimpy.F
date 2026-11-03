      Subroutine Getlstdimpy(Left_index,Right_index,Nrows,Ncols)

C The shape of MO list (Left_index,Right_index), i.e. MOIO(irrep,list):
C Nrows is the length of one column (moiosz) and Ncols the number of
C columns (moiods). Both are 0 for an index out of range or a list that
C has not been created.

      Implicit None

#include "lists.com"

      Integer*8 Left_index
      Integer*8 Right_index
      Integer*8 Nrows
      Integer*8 Ncols
Cf2py intent(out) Nrows
Cf2py intent(out) Ncols

      Nrows = 0
      Ncols = 0
      If (Left_index .Lt. 1 .Or. Left_index .Gt. _MAX_IO_GRPS .Or.
     +    Right_index .Lt. 1 .Or. Right_index .Gt. _MAX_IO_FAMS) Return
      If (moio(Left_index,Right_index) .Lt. 1) Return
      Nrows = moiosz(Left_index,Right_index)
      Ncols = moiods(Left_index,Right_index)

      Return
      End
