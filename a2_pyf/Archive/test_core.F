       Subroutine Test_core(Work,Maxcor)

       Implicit None
       Integer Maxcor
       Double precision Work
       Dimension Work(Maxcor)

       Write(6,"(a,i10)") "Double precision mem allocated", Maxcor
       Write(6,"(2(1x,F10.5))")  Work(1), Work(Maxcor)
       Write(6,"(2(1x,F10.5))") Work(Maxcor+5),Work(Maxcor+10)
      
       Return
       End 




