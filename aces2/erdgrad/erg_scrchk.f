










      Subroutine Erg_scrchk(Name,Imin,Liscr,Jmin,Lzscr)
C
C Stop when the integer or real scratch is smaller than the minimum
C an OED/ERD batch needs.
C
      Implicit Double Precision (A-H,O-Z)
      Character*(*) Name

      If (Imin .Gt. Liscr .Or. Jmin .Gt. Lzscr) Then
         Write(6,"(3a,4i12)") " @", Name, ": OED/ERD scratch too small",
     &                        Imin, Liscr, Jmin, Lzscr
         Call Errex
      Endif

      Return
      End
