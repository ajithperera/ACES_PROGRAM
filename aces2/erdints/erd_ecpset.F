      Subroutine Erd_ecpset(Nsh,Iatm,Lsh,Npr,Nct,Alp,Ioa,Cc,Ioc,Nharm)
C
C Load the basis set into ecplib's shell tables (/nshel/ and friends),
C what ecplib's Ecp_basis_init does from VMOL's common blocks. Every
C atom is its own "unique" atom here (xerdints has all of them), and
C the contraction coefficients are made the way VMOL stores them
C (vmol/readin.F): normalized coefficient * (1/2pi)**(3/4) *
C (4 alpha)**((L+1)/2 + 1/4), i.e. including VMOL's primitive
C normalization.
C
      Implicit Double Precision (A-H,O-Z)
#include "../ecplib/ecp.com"

      Integer Iatm(Nsh), Lsh(Nsh), Npr(Nsh), Nct(Nsh), Ioa(Nsh),
     &        Ioc(Nsh)
      Dimension Alp(*), Cc(*)
      Logical Nharm
      Common /palcali/ nshdim
      Common /palcal/ toll
      Common /symshe/ mulsh(Mxtnsh)
      Common /modez / zetm(Mxtnsh)

      If (Nsh .Gt. Mxtnsh) Then
         Write(6,"(a,i6)") " @Erd_ecpset: too many shells for ecplib ",
     &                     Nsh
         Call Errex
      Endif
      Pippi = (0.5D0/(4.0D0*Datan(1.0D0)))**0.75D0
      Nshell = Nsh
      Nbf    = 0
      Iprim  = 0
      Do Is = 1, Nsh
         L  = Lsh(Is) + 1
         Katom(Is)  = Iatm(Is)
         Ktype(Is)  = L
         Kmini(Is)  = (L-1)*L*(L+1)/6 + 1
         Kmaxi(Is)  = L*(L+1)*(L+2)/6
         If (Nharm) Then
            Ndegen(Is) = 2*L - 1
         Else
            Ndegen(Is) = L*(L+1)/2
         Endif
         Kbfn(Is)   = Nbf + 1
         Nbf        = Nbf + Ndegen(Is)
         Mulsh(Is)  = 1
         Nprims(Is) = Npr(Is)
         Kprim(Is)  = Iprim + 1
         Numcon(Is) = Nct(Is)
         Zetm(Is)   = Alp(Ioa(Is)+1)
         Do Ip = 1, Npr(Is)
            A = Alp(Ioa(Is)+Ip)
            If (Iprim+Ip .Gt. Mxtnpr) Then
               Write(6,"(a)") " @Erd_ecpset: too many primitives"
               Call Errex
            Endif
            Expnt(Iprim+Ip) = A
            Zetm(Is) = Min(Zetm(Is),A)
            Fac = Pippi*(4.0D0*A)**(0.5D0*Dble(L)+0.25D0)
            Do Ic = 1, Nct(Is)
               Contr(Iprim+Ip,Ic) = Cc(Ioc(Is)+(Ic-1)*Npr(Is)+Ip)*Fac
            Enddo
         Enddo
         Iprim = Iprim + Npr(Is)
      Enddo
      Nshdim = Nsh
      Toll   = 25.0D0

      Return
      End
