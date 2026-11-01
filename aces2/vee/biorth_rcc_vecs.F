      Subroutine Biorth_rcc_vecs(Evals,Evecr,Evecl,Matdim,Work,
     &                           Maxcor)
C
C DGEEV returns each left and right eigenvector scaled to unit length
C on its own. The doubles correction is formed as L(i)*W*R(i), which
C needs <L(i)|R(j)> = delta(ij). For a non-degenerate root this only
C means dividing L(i) by <L(i)|R(i)>. For a block of degenerate roots
C DGEEV picks the left and the right vectors of the block
C independently, so <L(i)|R(j)> is not diagonal there and L(i)*W*R(i)
C depends on that choice. Replace the left vectors of each block by
C the combinations that are biorthonormal to its right vectors,
C
C      L(block) <- L(block) * S**(-T),   S = L(block)**T * R(block).
C
C The eigenvalues must be sorted (Sort_rcc_eigs).
C
      Implicit Double Precision(A-H,O-Z)

      Dimension Evals(Matdim),Evecr(Matdim,Matdim)
      Dimension Evecl(Matdim,Matdim),Work(Maxcor)

      Data Degtol /1.0D-8/

      Ibgn = 1
      Do while (Ibgn .Le. Matdim)
         Iend = Ibgn
         Do while (Iend .Lt. Matdim)
            If (Dabs(Evals(Iend+1)-Evals(Ibgn)) .Gt. Degtol) Exit
            Iend = Iend + 1
         Enddo
         Nblk = Iend - Ibgn + 1

C S(k,l) = <L(k)|R(l)> in Work(I000), L(block)**T in Work(I010)

         I000 = 1
         I010 = I000 + Nblk*Nblk
         I020 = I010 + Nblk*Matdim
         Iend_w = I020 + Nblk
         If (Iend_w .Gt. Maxcor) Call Insmem("@-Biorth_rcc_vecs",
     &                                       Iend_w,Maxcor)

         Call Dgemm("T","N",Nblk,Nblk,Matdim,1.0D0,Evecl(1,Ibgn),
     &              Matdim,Evecr(1,Ibgn),Matdim,0.0D0,Work(I000),Nblk)
         Do K = 1, Nblk
            Call Dcopy(Matdim,Evecl(1,Ibgn+K-1),1,Work(I010+K-1),Nblk)
         Enddo

C S * X = L(block)**T gives X = L'(block)**T with L'**T * R = 1

         Call Dgesv(Nblk,Matdim,Work(I000),Nblk,Work(I020),Work(I010),
     &              Nblk,Info)
         If (Info .Ne. 0) Then
            Write(6,"(a,i4,a,i4)") " @-Biorth_rcc_vecs: <L|R> is",
     &                             " singular for the roots ",Ibgn
            Call Errex
         Endif

         Do K = 1, Nblk
            Call Dcopy(Matdim,Work(I010+K-1),Nblk,Evecl(1,Ibgn+K-1),1)
         Enddo

         Ibgn = Iend + 1
      Enddo

      Return
      End
