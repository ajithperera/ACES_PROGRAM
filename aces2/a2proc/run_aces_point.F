      Subroutine Run_aces_point(Script)
C
C Runs one ACES calculation for the PES_scan module (IRC, grid and
C streaming scans). If Script (runaces2a for the reference point,
C runaces2b for the others) is on the PATH it is used, so that a site
C can stage the run elsewhere; otherwise xaces2 runs in the current
C directory with its output in out.N, N = 0, 1, 2, ... The only copies
C of runaces2a/b were the 2007 QTP-cluster scripts in save/. ISHELL
C takes at most 80 characters, so the logic is here. 10/2026.
C
      Character*(*) Script
      Character*79 Cmd
      Integer ishell
      External ishell
      Save Npoint
      Data Npoint /0/
C
      Cmd = 'command -v '//Script//' > /dev/null 2>&1'
      If (ishell(Cmd) .Eq. 0) Then
         Call Runit(Script)
      Else
         Write(Cmd, "(a,I0,a)") 'xaces2 > out.', Npoint, ' 2>&1'
         Call Runit(Cmd)
      Endif
      Npoint = Npoint + 1
C
      Return
      End
