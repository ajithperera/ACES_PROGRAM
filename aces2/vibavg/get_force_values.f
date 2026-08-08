










      Subroutine Get_force_values(Data,D,Ndims,Iunit)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Data(Ndims,Ndims),D(Ndims)
      Character*80 Wrk

      Read(Iunit,*) Wrk
      Print*, Wrk

      Iblk_size = 6
      Nblocks = Ndims/Iblk_size 
      Nleft   = Ndims - Nblocks*Iblk_size
      Jc      = 2
      
      Do K = 1, Nblocks
         Ib = (K-1)*6 + 1
         Do I = Ib, Ndims+1
            Read(Iunit,"(6F14.10)") (D(J),J=Jc,Min(I,Jc+5))
            Do M = Jc, Min(I,Jc+5)
               Data(I-1,M-1) = D(M)
            Enddo
         Enddo
         Jc = Jc + 6
      Enddo

      If (Nleft .Ne. 0) Then
          Ib = Nblocks*Iblk_size + 1 
          Jc = Ib
          Do I = Ib, Ndims+1
             Read(Iunit,"(6F14.10)") (D(J),J=Jc,Min(I-1,Jc+5))
             Do M = Jc, Min(I-1,Jc+5)
               Data(I-1,M) = D(M)
            Enddo
          Enddo
      Endif

C Build the full force constant matrix (needs to have expanded version
C for reorientation transformations to work 

      Do I = 1, Ndims
         Do J = 1, Ndims
            Data(I,J) = Data(J,I)
         Enddo
      Enddo

      write(6,"(a,a)") "The force constant matrix:",Wrk
      call output(Data,1,Ndims,1,Ndims,Ndims,Ndims,1)

      Return
      End


