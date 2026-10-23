      Subroutine Erd_vhord(L,Iord)
C
C Column order of VMOL's spherical table for an L shell (L >= 2):
C Iord(j), j = 1..2L+1, are the table columns that are harmonics,
C then Iord(2L+2..) the contaminants, each kind in table order.
C
C VMOL's table interleaves harmonics and contaminants (for f, x r^2
C is column 6); the harmonics are the columns with zero Laplacian.
C Column j of the table is the function VMOL labels with Cartesian
C component j (readin.F, JPRX).
C
      Implicit Double Precision (A-H,O-Z)
      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2)
      Dimension Cv(Mxyz,Mxyz), Q(Mxyz)
      Integer Iord(*)

      Ncrt = (L+1)*(L+2)/2
      Nsph = 2*L + 1
      Call Erd_vharm(L,Cv,Mxyz)

      Nh = 0
      Nc = 0
      Do J = 1, Ncrt
         Call Erd_lapl(L,Cv(1,J),Q)
         Qn = 0.0D0
         Do K = 1, (L-1)*L/2
            Qn = Qn + Dabs(Q(K))
         Enddo
         If (Qn .Eq. 0.0D0) Then
            Nh = Nh + 1
            Iord(Nh) = J
         Else
            Nc = Nc + 1
            Iord(Nsph+Nc) = J
         Endif
      Enddo
      If (Nh .Ne. Nsph) Then
         Write(6,"(a,2i3)") " @Erd_vhord: harmonic count in VMOL ",
     &                      "table ", L, Nh
         Call Errex
      Endif

      Return
      End
