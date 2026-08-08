










      function streq(str1,str2,insensitive)







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




      logical streq
      character *(*) str1,str2
      logical insensitive
      integer l1,l2,i,strlen
      character*1 c1,c2
c      callstack_curr='STREQ'
      l1=strlen(str1)
      l2=strlen(str2)
      if (l1.ne.l2) then
        streq=.false.
        return
      endif
      do i=1,l1
        c1=str1(i:i)
        c2=str2(i:i)
        if (insensitive) then
          if (c1.ge.'a' .and. c1.le.'z')
     &        c1=char(ichar(c1)+ichar('A')-ichar('a'))
          if (c2.ge.'a' .and. c2.le.'z')
     &        c2=char(ichar(c2)+ichar('A')-ichar('a'))
        endif
        if (c1.ne.c2) then
          streq=.false.
          return
        endif
      enddo
      streq=.true.
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:

