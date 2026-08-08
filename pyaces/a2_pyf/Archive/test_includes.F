       Subroutine Test_includes()

       Implicit None
       Integer  i

       include "machsp.com"
       include "jobarc.com"
       include "icore.com"
       include "istart.com"
       include "flags.com"
       include "flags2.com"

       write(6,"(a)") "machsp.com"
       Write(6,"(a,i2)") "iintln: ", iintln
       Write(6,"(a,i2)") "ifltln: ", ifltln
       Write(6,"(a,i2)") "iintfp: ", iintfp
       Write(6,"(a,i6)") "ialone: ", ialone 
       Write(6,"(a,i2)") "ibitwd: ", ibitwd 

       write(6,"(a)") "ja_flags"
       write(6,"(a,l)")  "bjaup: ", bjaup
       write(6,"(a,l)")  "bjamod:", bjamod
       write(6,"(a)") "istart"
       write(6,"(a,es12.5e2)")  "io: ", Dble(i0)
       write(6,"(a,es12.5e2)")  "icrsiz:", Dble(icrsiz)
       write(6,"(a)") "flags: iflags(100)"
       write(6,"(20(1x,i2))") (iflags(i),i=1,100)
       write(6,"(a)") "flags2: iflags2(500)"
       write(6,"(20(1x,i2))") (iflags2(i),i=1,500)
       

       Return
       End 




