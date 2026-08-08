      Subroutine add(a,b,n)
      Implicit none
     
      Integer*8 n
      Integer*8 i
      Integer*8 m
      Integer*8 nset
      Logical Set_it

      parameter(nset=10)
      Double Precision a(n)
      Double Precision b(n)
      Double Precision a_1(nset)

      common /internal/A_1, set_it 
      save /internal/
      Do I = 1, n
         a(i) = a(i) + b(i)
      Enddo 
      write(6,"(a)") "@-add:The sum of the arrays"
      write(6,"(6(1x,F15.10))") (a(i),i=1,n)
      Write(6,"(i2)") n

      set_it = .True.
      If (set_it) Then
         Do i = 1, nset 
            a_1(i) = 1.0
         Enddo
      Endif 
      write(6,"(a)") "@-add:The common block vars"
      write(6,"(6(1x,F15.10))") (a_1(i),i=1,nset)
     
      write(6,"(a)") "@-add: calling substract"
      call substract(a,b,n)
      write(6,"(a)") "@-add: returned from substract"

      Write(6,"(a)") "@-add:The sum/substract of the arrays"
      write(6,"(6(1x,F15.10))") (a(i),i=1,n)

      Return
      End 
        
        


