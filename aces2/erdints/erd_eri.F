      Subroutine Erd_eri(List,Zscr,Lzscr,Iscr,Liscr,Alp,Cc,Icb,Ice,
     &                   Mxpk,Screen,Nbatch,Nfirst)
C
C One ERD batch (List(1) List(2)|List(3) List(4)) of contracted
C electron repulsion integrals, Zscr(Nfirst ...) with the functions of
C List(1) running fastest, over the basis functions the SOs are made
C of (Erd_fntr). The scratch needed is checked against
C ERD__MEMORY_ERI_BATCH every time.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Integer List(4), Iscr(Liscr), Icb(Mxpk), Ice(Mxpk)
      Dimension Zscr(Lzscr), Alp(Mxpk), Cc(Mxpk), X(4), Y(4), Z(4)
      Logical Screen

      Call Erd_pack(4,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,Mxpk)
      Do I = 1, 4
         Ia = Ishatm(List(I))
         X(I) = Xyzatm(1,Ia)
         Y(I) = Xyzatm(2,Ia)
         Z(I) = Xyzatm(3,Ia)
      Enddo
      I1 = List(1)
      I2 = List(2)
      I3 = List(3)
      I4 = List(4)

      Call Erd__memory_eri_batch(Nalp,Ncc,
     &     Nctsh(I1),Nctsh(I2),Nctsh(I3),Nctsh(I4),
     &     Nprsh(I1),Nprsh(I2),Nprsh(I3),Nprsh(I4),
     &     Lshl(I1),Lshl(I2),Lshl(I3),Lshl(I4),
     &     X(1),Y(1),Z(1),X(2),Y(2),Z(2),X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &     Alp,Cc,Spher,Imin,Iopt,Jmin,Jopt)
      If (Imin .Gt. Liscr .Or. Jmin .Gt. Lzscr) Then
         Write(6,"(a,4i12)") " @Erd_eri: ERD scratch too small ",
     &                       Imin, Liscr, Jmin, Lzscr
         Call Errex
      Endif
      Imax = Max(Imin,Min(Iopt,Liscr))
      Jmax = Max(Jmin,Min(Jopt,Lzscr))

      Call Erd__gener_eri_batch(Imax,Jmax,Nalp,Ncc,Ncsum,
     &     Nctsh(I1),Nctsh(I2),Nctsh(I3),Nctsh(I4),
     &     Nprsh(I1),Nprsh(I2),Nprsh(I3),Nprsh(I4),
     &     Lshl(I1),Lshl(I2),Lshl(I3),Lshl(I4),
     &     X(1),Y(1),Z(1),X(2),Y(2),Z(2),X(3),Y(3),Z(3),X(4),Y(4),Z(4),
     &     Alp,Cc,Icb,Ice,Spher,Screen,Iscr,Nbatch,Nfirst,Zscr)

C ERD functions -> the basis functions of the SOs (VMOL's).

      If (Nbatch .Gt. 0) Then
         N1 = Nfnsh(I1)
         N2 = Nfnsh(I2)
         N3 = Nfnsh(I3)
         N4 = Nfnsh(I4)
         Call Erd_fntr(Zscr(Nfirst),1,N1,N2*N3*N4,I1)
         Call Erd_fntr(Zscr(Nfirst),N1,N2,N3*N4,I2)
         Call Erd_fntr(Zscr(Nfirst),N1*N2,N3,N4,I3)
         Call Erd_fntr(Zscr(Nfirst),N1*N2*N3,N4,1,I4)
      Endif

      Return
      End
