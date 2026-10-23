      Subroutine Erd_irrfn(Iorb,Ipar,Irc,Isg)
C
C The SOs that one function of parity Ipar(3) (see Erd_cmppar) on the
C shells of orbit Iorb gives (same projection as Erd_sobld): the
C Norbm(Iorb) irreps Irc(j) (1 based) compatible with the stabilizer,
C and the sign Isg(m,j) of member m in the j-th of them.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Integer Ipar(3), Irc(8), Isg(8,8), Ipsh(0:7)
      Logical Compat

      Irep = Iorbsh(1,Iorb)
      Do Ir = 0, Nops-1
         Ipsh(Ir) = 1
         Do K = 1, 3
            If (Ipar(K) .Eq. 1) Ipsh(Ir) = Ipsh(Ir)*Iopsgn(K,Ir)
         Enddo
      Enddo
      Nc = 0
      Do Ig = 0, Nirr-1
         Compat = .True.
         Do Ir = 0, Nops-1
            Ichr = 1
            If (Mod(Popcnt(Iand(Ig,Ir)),2) .Eq. 1) Ichr = -1
            If (Ishmap(Irep,Ir) .Eq. Irep .And. Ichr .Ne. Ipsh(Ir))
     &         Compat = .False.
         Enddo
         If (Compat) Then
            Nc = Nc + 1
            Irc(Nc) = Ig + 1
            Do M = 1, Norbm(Iorb)
               Ir = Iorbop(M,Iorb)
               Ichr = 1
               If (Mod(Popcnt(Iand(Ig,Ir)),2) .Eq. 1) Ichr = -1
               Isg(M,Nc) = Ichr*Ipsh(Ir)
            Enddo
         Endif
      Enddo
      If (Nc .Ne. Norbm(Iorb)) Then
         Write(6,"(a,2i4)") " @Erd_irrfn: irrep count error ", Iorb, Nc
         Call Errex
      Endif

      Return
      End
