      Subroutine Erd_consrt(U,Ncrt,Nso,Nconir,Nirr,Key,Scr)
C
C Sort the contaminant columns of U (Nso+1..Ncrt, in irrep blocks of
C Nconir each) within each irrep by their first nonzero row, VMOL's
C order for the deleted functions of FULLSOAO.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension U(Ncrt,Ncrt), Scr(Ncrt,*)
      Integer Nconir(Nirr), Key(Ncrt)

      Do J = Nso+1, Ncrt
         Key(J) = 0
         Do I = Ncrt, 1, -1
            If (U(I,J) .Ne. 0.0D0) Key(J) = I
         Enddo
      Enddo
      J0 = Nso
      Do Ig = 1, Nirr
         Do J = J0+1, J0+Nconir(Ig)
            Irank = 0
            Do J2 = J0+1, J0+Nconir(Ig)
               If (Key(J2) .Lt. Key(J)) Irank = Irank + 1
            Enddo
            Call Dcopy(Ncrt,U(1,J),1,Scr(1,Irank+1),1)
         Enddo
         Do J = 1, Nconir(Ig)
            Call Dcopy(Ncrt,Scr(1,J),1,U(1,J0+J),1)
         Enddo
         J0 = J0 + Nconir(Ig)
      Enddo

      Return
      End
