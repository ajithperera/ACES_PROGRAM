      Subroutine Test_py(Work,Maxcor,Scr,Ndim,Iuhf)

      Implicit None
      Integer Iuhf
      Integer Maxcor
      Integer Noz
      Integer Ntw
      Integer Rec_len
      Integer Ioff
      Integer I
      Integer J
      Integer Ndim

      Double Precision Work(Maxcor)
      Double Precision Scr(Ndim)

      Data Noz,Ntw /0,20/

      Call Getrec(Noz,"JOBARC","COORD   ",Rec_len,Work)

      If (Rec_len .Gt. Noz) Then
          Write(6,"(a,1x,I3,1x,a)") "Record COORD with length",Rec_len,
     +                              "is found"
          
          Call Getrec(Ntw,"JOBARC","COORD   ",Rec_len,Work)
          Ioff = 0
          Do I = 1, Rec_len/3
             Write(6,"(3(1x,F15.9))") (Work(Ioff+J),J=1,3)
             Ioff = Ioff  + 3
          Enddo
          Call Dcopy(Ndim,Work,1,Scr,1)
      Endif 
      
      Return
      End 
