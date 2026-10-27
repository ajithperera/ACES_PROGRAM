      subroutine run_symadapt()
c----------------------------------------------------------------------------
c   At the end of a job: runs "xa2proc symadapt", which gives the C1 SCF
c   orbitals (and the excited state natural orbitals, when their densities
c   are on JOBARC) the labels of the full point group and appends them to
c   summary.out. xa2proc is $ACES_XA2PROC, else next to $ACES_XJODA, else
c   on PATH. A failure only prints a warning. 10/2026.
c----------------------------------------------------------------------------
      implicit none

      character*512 xa2proc
      integer ilen, istat, i
      logical is_open

      call get_environment_variable('ACES_XA2PROC', xa2proc, ilen,
     *                              istat)
      if (istat .ne. 0 .or. ilen .eq. 0) then
         call get_environment_variable('ACES_XJODA', xa2proc, ilen,
     *                                 istat)
         if (istat .eq. 0 .and. ilen .gt. 0) then
            i = index(xa2proc(1:ilen), '/', back=.true.)
            xa2proc = xa2proc(1:i)//'xa2proc'
         else
            xa2proc = 'xa2proc'
         endif
      endif

c   summary.out is unit 66 in xaces3; close it so xa2proc can append.
      inquire(unit=66, opened=is_open)
      if (is_open) close(66)

      call c_flush_stdout()
      call execute_command_line(trim(xa2proc)//' symadapt',
     *                          exitstat=istat)
      call c_flush_stdout()
      if (istat .ne. 0) print *,'@RUN_SYMADAPT: warning: ',
     *   trim(xa2proc),' symadapt failed, status ',istat

      return
      end
