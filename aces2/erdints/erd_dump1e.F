      Subroutine Erd_dump1e(S,T,V,N)
C
C Debug: write the shell table and the ERD-order AO one-electron
C matrices to ERDAO1E.TXT / ERDAO1E.BIN (stream) for checking against
C VMOL.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension S(N,N), T(N,N), V(N,N)

      Open(Unit=31,File='ERDAO1E.TXT',Form='Formatted',
     &     Status='Unknown')
      Write(31,"(3i8,l3)") Natm, Nshl, N, Spher
      Do Is = 1, Nshl
         Write(31,"(7i8)") Is, Ishatm(Is), Lshl(Is), Nprsh(Is),
     &                     Nctsh(Is), Nfnsh(Is), Ioffn(Is)
      Enddo
      Write(31,"(2i8)") Norb, Nso
      Do Io = 1, Norb
         Write(31,"(10i8)") Io, Norbm(Io), (Iorbsh(M,Io),M=1,Norbm(Io)),
     &                      (Iorbop(M,Io),M=1,Norbm(Io))
      Enddo
      Do Ip = 1, Nso
         Io = Isoorb(Ip)
         Write(31,"(12i6)") Ip, Isoirr(Ip), Io, Isofn(Ip),
     &                      (Isosgn(M,Ip),M=1,Norbm(Io))
      Enddo
      Close(31)
      Open(Unit=31,File='ERDAO1E.BIN',Form='Unformatted',
     &     Access='Stream',Status='Replace')
      Write(31) S, T, V
      Close(31)

      Return
      End
