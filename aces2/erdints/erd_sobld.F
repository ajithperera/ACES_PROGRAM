      Subroutine Erd_sobld
C
C Build the SO basis from the shell orbits (fills the erdso common).
C
C For a function mu on the representative shell A of an orbit, the
C projector onto irrep G gives
C
C    P_G mu_A = sum_r chi_G(r) O_r mu_A = sum_r chi_G(r) p_mu(r) mu_r(A)
C
C where p_mu(r) = +-1 is the parity of mu under r (Erd_cmppar) and
C chi_G(r) = (-1)**popcount(iand(G,r)). It vanishes unless
C chi_G(h) = p_mu(h) for every h of the stabilizer H of A, which
C leaves exactly |orbit| irreps. For those, the coefficient of the
C member B = r_B(A) is |H| chi_G(r_B) p_mu(r_B); dropping |H| gives
C VMOL's unnormalized SOs, sum over members of +-mu_B.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Integer Ipar(3,64), Ipsh(0:7,64), Ichr(0:7,0:7), Nfill(8),
     &        Iso0(8)
      Dimension Fac(64), Cv(36,36)
      Integer Iparv(3,36)
      Logical Compat

      Nirr = Nops
      Nso  = Nfnerd
      If (Nso .Gt. Mxfn) Then
         Write(6,"(a,i6)") " @Erd_sobld: too many functions, max ",
     &                     Mxfn
         Call Errex
      Endif

C Characters.

      Do Ig = 0, Nirr-1
         Do Ir = 0, Nops-1
            Ichr(Ir,Ig) = 1
            If (Mod(Popcnt(Iand(Ig,Ir)),2) .Eq. 1) Ichr(Ir,Ig) = -1
         Enddo
      Enddo

C Shell orbits; members in the order of the operation generating
C them first.

      Do Is = 1, Nshl
         Korbsh(Is) = 0
      Enddo
      Norb = 0
      Do Is = 1, Nshl
         If (Korbsh(Is) .Eq. 0) Then
            Norb = Norb + 1
            Norbm(Norb) = 0
            Do Ir = 0, Nops-1
               Js = Ishmap(Is,Ir)
               If (Korbsh(Js) .Eq. 0) Then
                  Norbm(Norb) = Norbm(Norb) + 1
                  Iorbsh(Norbm(Norb),Norb) = Js
                  Iorbop(Norbm(Norb),Norb) = Ir
                  Korbsh(Js) = Norb
               Endif
            Enddo
         Endif
      Enddo

C First pass: count the SOs of each irrep. Second pass: assign them,
C irrep by irrep, in orbit order and then function order.

      Do Ipass = 1, 2
         Do Ig = 1, Nirr
            Nfill(Ig) = 0
         Enddo
         If (Ipass .Eq. 2) Then
            Iso0(1) = 0
            Do Ig = 2, Nirr
               Iso0(Ig) = Iso0(Ig-1) + Nbfirr(Ig-1)
            Enddo
         Endif
         Nlist = 0

         Do Iorb = 1, Norb
            Irep = Iorbsh(1,Iorb)
            Nm   = Norbm(Iorb)
            L    = Lshl(Irep)
            Nct  = Nctsh(Irep)
            Call Erd_cmppar(L,Ipar,Ncmp)
            Do Ic = 1, Ncmp
               Do Ir = 0, Nops-1
                  Ipsh(Ir,Ic) = 1
                  Do K = 1, 3
                     If (Ipar(K,Ic) .Eq. 1)
     &                  Ipsh(Ir,Ic) = Ipsh(Ir,Ic)*Iopsgn(K,Ir)
                  Enddo
               Enddo
            Enddo
            If (Ipass .Eq. 2) Iorbso(Iorb) = Nlist

C Within a shell ERD runs the components fastest (contraction major).

            Do If = 1, Nfnsh(Irep)
               Ic = Mod(If-1,Ncmp) + 1
               Ncomp = 0
               Do Ig = 0, Nirr-1
                  Compat = .True.
                  Do Ir = 0, Nops-1
                     If (Ishmap(Irep,Ir) .Eq. Irep .And.
     &                   Ichr(Ir,Ig) .Ne. Ipsh(Ir,Ic)) Compat = .False.
                  Enddo
                  If (Compat) Then
                     Ncomp = Ncomp + 1
                     Nfill(Ig+1) = Nfill(Ig+1) + 1
                     If (Ipass .Eq. 2) Then
                        Ip = Iso0(Ig+1) + Nfill(Ig+1)
                        Isoirr(Ip) = Ig + 1
                        Isoorb(Ip) = Iorb
                        Isofn(Ip)  = If
                        Do M = 1, Nm
                           Ir = Iorbop(M,Iorb)
                           Isosgn(M,Ip) = Ichr(Ir,Ig)*Ipsh(Ir,Ic)
                        Enddo
                        Nlist = Nlist + 1
                        Lorbso(Nlist) = Ip
                     Endif
                  Endif
               Enddo
               If (Ncomp .Ne. Nm) Then
                  Write(6,"(a,3i5)") " @Erd_sobld: irrep count error ",
     &                               Iorb, If, Ncomp
                  Call Errex
               Endif
            Enddo
         Enddo

         If (Ipass .Eq. 1) Then
            Do Ig = 1, Nirr
               Nbfirr(Ig) = Nfill(Ig)
            Enddo
         Endif
      Enddo

C Scale factors between ERD's Cartesian functions and VMOL's, and for
C a spherical basis the ERD -> VMOL harmonic transformations.

      If (Spher) Then
         Lmx = 0
         Do Is = 1, Nshl
            Lmx = Max(Lmx,Lshl(Is))
         Enddo
         Do L = 2, Lmx
            Call Erd_vsph(L,Cv,36,Wvsph(1,1,L),13,Iparv)
         Enddo
      Endif

      Do Is = 1, Nshl
         L = Lshl(Is)
         Call Erd_cmpfac(L,Spher,Fac,Ncmp)
         Do If = 1, Nfnsh(Is)
            Fnfac(Ioffn(Is)+If) = Fac(Mod(If-1,Ncmp)+1)
         Enddo
      Enddo

      If (Nlist .Ne. Nso) Then
         Write(6,"(a,2i6)") " @Erd_sobld: SO count error ", Nlist, Nso
         Call Errex
      Endif

      Return
      End
