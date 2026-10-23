      Subroutine Erd_cmppar(L,Ipar,Ncmp)
C
C Parity of each basis function of an angular momentum L shell under
C the sign flips x->-x, y->-y, z->-z: Ipar(k,i) = 1 when function i is
C odd in coordinate k, 0 when even. Every D2h operation therefore
C multiplies function i by prod_k Iopsgn(k,r)**Ipar(k,i): basis
C functions never mix, they only change sign.
C
C Cartesian functions are the monomials x^p y^q z^r in ERD's p>q>r
C order. Spherical ones (L >= 2) are VMOL's harmonics in VMOL's order
C (Erd_vsph), the functions the SO basis is built from.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2)
      Integer Ipar(3,*), Ixyz(3,Mxyz), Iparv(3,Mxyz)
      Dimension Cv(Mxyz,Mxyz), W(2*Lmax+1,2*Lmax+1)

      If (L .Gt. Lmax) Then
         Write(6,"(a,i3)") " @Erd_cmppar: angular momentum too high ",
     &                     L
         Call Errex
      Endif

      If (.Not. Spher .Or. L .Le. 1) Then
         Call Erd_monos(L,Ixyz,Ncmp)
         Do I = 1, Ncmp
            Do K = 1, 3
               Ipar(K,I) = Mod(Ixyz(K,I),2)
            Enddo
         Enddo
      Else
         Ncmp = 2*L + 1
         Call Erd_vsph(L,Cv,Mxyz,W,2*Lmax+1,Iparv)
         Do I = 1, Ncmp
            Do K = 1, 3
               Ipar(K,I) = Iparv(K,I)
            Enddo
         Enddo
      Endif

      Return
      End
