










      Subroutine Erd_csyord(Csy,Ncrt,Nfao,Nirr,Key,Scr)
C
C Put the Cartesian SOs of Csy (columns, irrep blocks of Nfao each,
C entries +-1 on computational AO rows) in VMOL's order and phase
C (CSYMTRAN): within an irrep by the first AO row they contain, with
C that coefficient +1; the same rule Erd_soord applies to the SOs.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Csy(Ncrt,Ncrt), Scr(Ncrt,Ncrt)
      Integer Nfao(Nirr), Key(Ncrt)

      Do J = 1, Ncrt
         Key(J) = 0
         Do I = Ncrt, 1, -1
            If (Csy(I,J) .Ne. 0.0D0) Key(J) = I
         Enddo
         If (Csy(Key(J),J) .Lt. 0.0D0) Then
            Do I = 1, Ncrt
               Csy(I,J) = -Csy(I,J)
            Enddo
         Endif
      Enddo
      J0 = 0
      Do Ig = 1, Nirr
         Do J = J0+1, J0+Nfao(Ig)
            Irank = 0
            Do J2 = J0+1, J0+Nfao(Ig)
               If (Key(J2) .Lt. Key(J)) Irank = Irank + 1
            Enddo
            Call Dcopy(Ncrt,Csy(1,J),1,Scr(1,J0+Irank+1),1)
         Enddo
         J0 = J0 + Nfao(Ig)
      Enddo
      Call Dcopy(Ncrt*Ncrt,Scr,1,Csy,1)

      Return
      End
