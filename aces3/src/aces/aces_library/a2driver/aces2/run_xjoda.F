      subroutine run_xjoda()
c----------------------------------------------------------------------------
c   Runs the ACES II xjoda in the work directory as "xjoda aces3", which
c   reads ZMAT (or continues an optimization or finite-difference series
c   from JOBARC) in C1 and in the input orientation. xjoda is taken from
c   $ACES_XJODA if it is set, otherwise from PATH. ACES III used to carry
c   its own copy of joda (geopt); 10/2026.
c----------------------------------------------------------------------------
      implicit none

      character*512 xjoda
      integer ilen, istat

      call get_environment_variable('ACES_XJODA', xjoda, ilen, istat)
      if (istat .ne. 0 .or. ilen .eq. 0) xjoda = 'xjoda'

      call c_flush_stdout()
      call execute_command_line(trim(xjoda)//' aces3', exitstat=istat)
      call c_flush_stdout()
      if (istat .ne. 0) then
         print *,'@RUN_XJODA: ',trim(xjoda),' aces3 failed, status ',
     *           istat
         call abort_job()
      endif

      return
      end
