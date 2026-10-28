      Subroutine Erg_orbit(List,Nlist,Lrep,Nimg)
C
C Orbit of a shell pair (Nlist=2) or shell quartet (Nlist=4) under
C the point group. List is canonical: List(1) >= List(2) and, for a
C quartet, List(3) >= List(4) and pair(12) >= pair(34). Every image
C is brought to the same canonical form and coded as an integer.
C
C Lrep : the list is the representative of its orbit (smallest code)
C Nimg : number of distinct images (the orbit size)
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Integer List(Nlist), Icode(0:7)
      Logical Lrep

      Do Iop = 0, Nops-1
         I = Ishmap(List(1),Iop)
         J = Ishmap(List(2),Iop)
         Ij = Max(I,J)*(Max(I,J)-1)/2 + Min(I,J)
         If (Nlist .Eq. 2) Then
            Icode(Iop) = Ij
         Else
            K = Ishmap(List(3),Iop)
            L = Ishmap(List(4),Iop)
            Kl = Max(K,L)*(Max(K,L)-1)/2 + Min(K,L)
            Icode(Iop) = Max(Ij,Kl)*(Max(Ij,Kl)-1)/2 + Min(Ij,Kl)
         Endif
      Enddo

      Lrep = .True.
      Nimg = 0
      Do Iop = 0, Nops-1
         If (Icode(Iop) .Lt. Icode(0)) Lrep = .False.
         Lnew = 1
         Do Jop = 0, Iop-1
            If (Icode(Jop) .Eq. Icode(Iop)) Lnew = 0
         Enddo
         Nimg = Nimg + Lnew
      Enddo

      Return
      End
