










      Subroutine Read_amps(Iunit,T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b,Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Work(Maxcor)

      Integer I,J,A,B

      Write(6,"(2a)") "The starting amplitudes are read from the",
     +                " TGUESS file."
 
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         Read(Iunit) T1aa(A,I)
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         Read(Iunit) T1bb(a,i)
      Enddo
      Enddo

      do i=1,Nocc_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,i,T1aa(a,i)
      enddo
      enddo

      Return
      End
