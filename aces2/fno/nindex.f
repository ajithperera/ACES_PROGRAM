










      function nindex(str,substr,ind)







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




      integer nindex
      character*(*) str,substr
      integer ind,strlen,sublen,strmax,i,ncindex,slen
c      callstack_curr='NINDEX'
      slen=strlen(str)
      sublen=strlen(substr)
      strmax=slen-sublen+1
      if (sublen.gt.slen .or. ind.ge.strmax .or. ind.lt.0) then
        nindex=0
        return
      end if
      i=ind
 10   nindex=ncindex(str,substr(1:1),i)
      if (nindex.eq.0) return
      if (nindex.gt.strmax) then
        nindex=0
        return
      end if
      if (str(nindex:nindex+sublen-1).eq.substr(1:sublen)) return
      i=nindex
      goto 10
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
