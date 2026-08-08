










      Subroutine Pccd_insist_no_oovv(Htau_pq,Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Dimension Htau_pq(Nbas,Nbas)
      Integer P,Q,A,B,I,J
 
      Data Dnull /0.0D0/

      Do I = 1, Nocc
         Do J = 1, Nocc
            Htau_pq(I,J) = Dnull
            Htau_pq(J,I) = Dnull
         Enddo
      Enddo 

      Do A = 1+Nocc,Nbas
         Do B = 1+Nocc,Nbas 
            Htau_pq(A,B) = Dnull
            Htau_pq(B,A) = Dnull
         Enddo
      Enddo 


      Return
      End 
