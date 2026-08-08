











      subroutine endscf
      implicit none
      integer iRc, iOption

c molcas.com : begin
      logical seward, petite_list
      character*8 fnord
      integer luord
      integer ipmat(8,2), lbbt, lbbs
      common /molcas_com/ seward, petite_list, fnord, luord,
     &                    ipmat, lbbt, lbbs
c molcas.com : end

c In pyaces's single-process model, do not finalize (close) the job
c archive/I-O subsystem here -- the caller (Python, via runscf.F) may
c still need to read records after Runscf returns, matching runscf.F's
c own existing "do not call aces_fin" decision for the same reason.

      if (seward) then
c      o $MOLCAS/src/scf/scf.f finalization routines... (fatal)
         iRc = -1
         iOption = 0
         call clsord(iRc,iOption)
         if (iRc.ne.0) then
            write(*,*) '@ENDSCF: There was an error closing the Molcas',
     &                 ' 2-electron integral file.'
            call errex
         end if
cSEWTRACE         call qexit("VSCF")
cSEWARD         call finish(itmp)
      end if

      return
      end

