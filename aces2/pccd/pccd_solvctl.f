










      Subroutine Pccd_solvctl(Alpha,Nodiis)
C
C Debug hook for the pCCD amplitude solver (Pccd_dodiis0): if the file
C pccd_solver.in exists its first line gives
C   Alpha Nodiis
C Alpha: damping of the Jacobi update, T = T(old) + Alpha*dT (default 1)
C Nodiis: 1 = no DIIS extrapolation (default 0).
C The file is read once; the values are reported on the first call.
C
      Implicit Double Precision(A-H,O-Z)
      Logical Exist,Done
      Save Done,Alpha0,Nodii0
      Data Done /.False./

      If (.Not. Done) Then
         Done   = .True.
         Alpha0 = 1.0D0
         Nodii0 = 0
         Inquire(File="pccd_solver.in",Exist=Exist)
         If (Exist) Then
            Open(Unit=41,File="pccd_solver.in",Status="Old")
            Read(41,*) Alpha0,Nodii0
            Close(41)
            Write(6,"(a,F8.4,a,i2)") " @-Pccd_solvctl: damping ",
     +            Alpha0,", no DIIS =",Nodii0
         Endif
      Endif
      Alpha  = Alpha0
      Nodiis = Nodii0
      Return
      End
