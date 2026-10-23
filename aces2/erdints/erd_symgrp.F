      Subroutine Erd_symgrp(Kgen)
C
C Set up the computational point group from the VMOL generators in
C MOL_ERDOED. ACES II works in D2h or one of its subgroups, so every
C operation is diag(+-1,+-1,+-1). A generator string lists the
C coordinates it flips ('X  ' is the yz mirror, 'XY ' is C2(z), 'XYZ'
C is the inversion). Operation r (0 <= r < 2**Ngen) is the product of
C the generators whose bit is set in r, the same numbering VMOL uses
C (readin.F, KA array), so irrep Gamma (0 based) has character
C (-1)**popcount(iand(Gamma,r)) and is ACES irrep Gamma+1.
C
C Also maps every atom and shell onto its image under each operation.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Character*3 Kgen(3)
      Character*1 Kxyz(3)
      Integer Igsgn(3,3), Ifirst(Mxatms+1)
      Data Kxyz /'X','Y','Z'/
      Data Tol /1.0D-6/

      Do Ig = 1, 3
         Do K = 1, 3
            Igsgn(K,Ig) = 1
         Enddo
      Enddo
      Do Ig = 1, Ngen
         Do Ic = 1, 3
            Do K = 1, 3
               If (Kgen(Ig)(Ic:Ic) .Eq. Kxyz(K)) Igsgn(K,Ig) = -1
            Enddo
         Enddo
      Enddo

      Nops = 2**Ngen
      Do Ir = 0, Nops-1
         Do K = 1, 3
            Iopsgn(K,Ir) = 1
            Do Ig = 1, Ngen
               If (Btest(Ir,Ig-1)) Iopsgn(K,Ir) = Iopsgn(K,Ir)*
     &                                            Igsgn(K,Ig)
            Enddo
         Enddo
      Enddo

C Atom images.

      Do Ir = 0, Nops-1
         Do Ia = 1, Natm
            Jmatch = 0
            Do Ja = 1, Natm
               Dist = 0.0D0
               Do K = 1, 3
                  Dist = Dist + Dabs(Xyzatm(K,Ja) -
     &                          Dble(Iopsgn(K,Ir))*Xyzatm(K,Ia))
               Enddo
               If (Dist .Lt. Tol .And. Zatm(Ja) .Eq. Zatm(Ia))
     &            Jmatch = Ja
            Enddo
            If (Jmatch .Eq. 0) Then
               Write(6,"(a,i4,a,i2)") " @Erd_symgrp: atom ", Ia,
     &               " has no image under operation ", Ir
               Call Errex
            Endif
            Iatmap(Ia,Ir) = Jmatch
         Enddo
      Enddo

C Shell images. The shells of an atom are contiguous and equivalent
C atoms carry the same basis, so the image of the k-th shell of atom
C a is the k-th shell of the image of a. Check that it really is.

      Do Ia = 1, Natm+1
         Ifirst(Ia) = 0
      Enddo
      Do Is = Nshl, 1, -1
         Ifirst(Ishatm(Is)) = Is
      Enddo
      Do Ir = 0, Nops-1
         Do Is = 1, Nshl
            Ia = Ishatm(Is)
            Ja = Iatmap(Ia,Ir)
            Js = Ifirst(Ja) + (Is - Ifirst(Ia))
            Ibad = 0
            If (Js .Gt. Nshl) Then
               Ibad = 1
            Else If (Ishatm(Js) .Ne. Ja .Or. Lshl(Js) .Ne. Lshl(Is)
     &          .Or. Nprsh(Js) .Ne. Nprsh(Is) .Or. Nctsh(Js) .Ne.
     &          Nctsh(Is)) Then
               Ibad = 1
            Else
               Do Ip = 1, Nprsh(Is)
                  If (Alpha(Iofalp(Js)+Ip) .Ne. Alpha(Iofalp(Is)+Ip))
     &               Ibad = 1
               Enddo
            Endif
            If (Ibad .Ne. 0) Then
               Write(6,"(a,i5,a)") " @Erd_symgrp: shell ", Is,
     &               " has no matching image; equivalent atoms"
               Write(6,"(a)") "  must carry the same basis set."
               Call Errex
            Endif
            Ishmap(Is,Ir) = Js
         Enddo
      Enddo

      Return
      End
