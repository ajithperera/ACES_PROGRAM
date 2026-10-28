      Subroutine Pccd_zt1(Work,Maxcor)
C
C Debug hook (file pccd_zt1.in, RHF): print max|T1| (list 90) and
C max|L1| (list 190) before the densities are built; with "ZERO" on the
C first line also set both to zero. Without _T1ON the T1 equations are
C not solved and T1 is the f(i,a)/D(i,a) guess of the current orbitals,
C which still enters Gamma(Ij,Ab) (list 116) through Pccd_h4x2all.
C
      Implicit Double Precision(A-H,O-Z)
      Character*8 Mode
      Logical Exist
      Dimension Work(Maxcor),Lists(2)
#include "sym.com"
      Data Lists /90,190/

      Inquire(File="pccd_zt1.in",Exist=Exist)
      If (.Not. Exist) Return
      Open(Unit=44,File="pccd_zt1.in",Status="Old")
      Read(44,"(a)") Mode
      Close(44)
      N = Nt(1)
      If (N .Gt. Maxcor) Call Insmem("Pccd_zt1",N,Maxcor)
      Do Il = 1, 2
         Call Getlst(Work,1,1,1,1,Lists(Il))
         Amax = 0.0D0
         Do I = 1, N
            Amax = Max(Amax,Dabs(Work(I)))
         Enddo
         Write(6,"(a,i4,a,1p,E12.4)") " @-Pccd_zt1: list",Lists(Il),
     +                                " max|amp| =",Amax
         If (Mode(1:4) .Eq. "ZERO") Then
            Call Dzero(Work,N)
            Call Putlst(Work,1,1,1,1,Lists(Il))
         Endif
      Enddo
      Return
      End
