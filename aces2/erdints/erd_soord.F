      Subroutine Erd_soord(Icmpao,Ncrt,Key,Inew,Iscr,Nsc)
C
C Put the SOs of each irrep in VMOL's order and phase: VMOL sorts them
C by the computational index of their function on the first member
C (orbit, then l, component, contraction, Erd_cmpord) and gives that
C member the + sign. Spherical functions (L >= 2) are ordered the same
C way with the harmonic (VMOL's order, Erd_vsph) in place of the
C Cartesian component. This reproduces VMOL's SYMTRANS (AO2SO).
C
C Rewrites Isoirr, Isoorb, Isofn, Isosgn and Lorbso (erdso common).
C Key, Inew are Nso long, Iscr is Nsc = 11*Nso long scratch.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Integer Icmpao(Ncrt), Key(Nso), Inew(Nso), Iscr(Nsc)
      Logical Lsph

C Sort key of every SO and the member that carries it.

      Do Ip = 1, Nso
         Io  = Isoorb(Ip)
         Is0 = Iorbsh(1,Io)
         L   = Lshl(Is0)
         Nk  = (L+1)*(L+2)/2
         Nc  = Nfnsh(Is0)/Nctsh(Is0)
         It  = (Isofn(Ip)-1)/Nc + 1
         Mc  = Mod(Isofn(Ip)-1,Nc) + 1
         Kmin = Ncrt + 1
         Mmin = 0
         Lsph = Spher .And. L .Ge. 2
         Do M = 1, Norbm(Io)
            If (Lsph) Then
               Iao = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk + 1
            Else
               Iao = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk + Mc
            Endif
            If (Icmpao(Iao) .Lt. Kmin) Then
               Kmin = Icmpao(Iao)
               Mmin = M
            Endif
         Enddo
         If (.Not. Lsph) Then
            Key(Ip) = Kmin*10000
         Else

C VMOL orders spherical functions like Cartesian ones: harmonic, then
C all contractions of that L on the atom. Key: start of the atom's L
C block, harmonic, contraction counted over the atom's L shells.

            Is  = Iorbsh(Mmin,Io)
            Ia  = Ishatm(Is)
            Nt  = 0
            Ito = 0
            Ifs = 0
            Do Js = 1, Nshl
               If (Ishatm(Js) .Eq. Ia .And. Lshl(Js) .Eq. L) Then
                  If (Ifs .Eq. 0) Ifs = Js
                  If (Js .Lt. Is) Ito = Ito + Nctsh(Js)
                  Nt = Nt + Nctsh(Js)
               Endif
            Enddo
            Kbl = Ncrt + 1
            Do M = 1, Norbm(Io)
               Ks = Ishmap(Ifs,Iorbop(M,Io))
               Kbl = Min(Kbl,Icmpao(Iofcrt(Ks)+1))
            Enddo
            Key(Ip) = Kbl*10000 + (Mc-1)*Nt + Ito + It
         Endif
         If (Isosgn(Mmin,Ip) .Lt. 0) Then
            Do M = 1, Norbm(Io)
               Isosgn(M,Ip) = -Isosgn(M,Ip)
            Enddo
         Endif
      Enddo

C New position of every SO: rank of its key within its irrep (the
C irreps keep their blocks).

      Iso0 = 0
      Do Ig = 1, Nirr
         Do Ip = Iso0+1, Iso0+Nbfirr(Ig)
            Irank = 0
            Do Iq = Iso0+1, Iso0+Nbfirr(Ig)
               If (Key(Iq) .Lt. Key(Ip)) Irank = Irank + 1
            Enddo
            Inew(Ip) = Iso0 + Irank + 1
         Enddo
         Iso0 = Iso0 + Nbfirr(Ig)
      Enddo

C Permute the SO tables (Iscr holds the old ones).

      Do Ip = 1, Nso
         Iscr(Ip)       = Isoirr(Ip)
         Iscr(Nso+Ip)   = Isoorb(Ip)
         Iscr(2*Nso+Ip) = Isofn(Ip)
         Do M = 1, 8
            Iscr(3*Nso+(Ip-1)*8+M) = Isosgn(M,Ip)
         Enddo
      Enddo
      Do Ip = 1, Nso
         Jp = Inew(Ip)
         Isoirr(Jp) = Iscr(Ip)
         Isoorb(Jp) = Iscr(Nso+Ip)
         Isofn(Jp)  = Iscr(2*Nso+Ip)
         Do M = 1, 8
            Isosgn(M,Jp) = Iscr(3*Nso+(Ip-1)*8+M)
         Enddo
      Enddo
      Do I = 1, Nso
         Lorbso(I) = Inew(Lorbso(I))
      Enddo

      Return
      End
