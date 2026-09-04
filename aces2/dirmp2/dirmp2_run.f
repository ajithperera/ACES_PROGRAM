










      subroutine dirmp2_run (irnk)

c Everything xdirmp2 does between setting up the runtime and tearing it
c down again. Split out of the program body so it can also be driven in
c process: pyaces has already done aces_init_rte/MPW_Init and will do
c aces_ja_fin/MPW_Finalize itself, and calling the program is not an
c option from there. Same wrapper-split as joda/vscf/vtran/vprops.
c The classic xdirmp2 is unchanged -- its program unit now calls this.

      implicit double precision (a-h,o-z)
      common // icore(1)
c parallel_aces.com : begin

c This common block contains the MPI statistics for each MPI process. The values
c are initialized in the acescore library.

      external aces_bd_parallel_aces




      integer                nprocs, irank, icpuname

      character*(256) szcpuname

      common /parallel_aces/ nprocs, irank, icpuname,
     &                       szcpuname
      save   /parallel_aces/

c parallel_aces.com : end
      character*3 str

      call nr2asc (irnk,str)
      open (96,file='OUT.'//str)
      open (97,file='DUMP.'//str)
      write (96,'(a,2i3)') 'nprocs,irank:',nprocs,irnk

      call CRAPSI (icore(1),iuhf,-1)

      if (iuhf.eq.0) then
         write (96,'(a)') 'CLOSED SHELL DIRECT MP2 CALCULATION'
         call mp2_rhf
      else
         write (96,'(a)') 'OPEN SHELL DIRECT MP2 CALCULATION'
         call mp2_uhf
      end if

      return
      end
