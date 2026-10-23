      Subroutine Erd_pack(Nlist,List,Alp,Nalp,Cc,Ncc,Icb,Ice,Ncsum,
     &                    Mxpk)
C
C Collect exponents, contraction coefficients and the first/last
C nonzero primitive of each contraction for the Nlist (2 or 4) shells
C in List, concatenated in that order, as the OED/ERD batch routines
C want them.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Integer List(Nlist), Icb(Mxpk), Ice(Mxpk)
      Dimension Alp(Mxpk), Cc(Mxpk)

      Nalp  = 0
      Ncc   = 0
      Ncsum = 0
      Do Il = 1, Nlist
         Is = List(Il)
         Np = Nprsh(Is)
         Nc = Nctsh(Is)
         If (Nalp+Np .Gt. Mxpk .Or. Ncc+Np*Nc .Gt. Mxpk .Or.
     &       Ncsum+Nc .Gt. Mxpk) Then
            Write(6,"(a)") " @Erd_pack: packing buffer too small"
            Call Errex
         Endif
         Do Ip = 1, Np
            Alp(Nalp+Ip) = Alpha(Iofalp(Is)+Ip)
         Enddo
         Do Ic = 1, Np*Nc
            Cc(Ncc+Ic) = Ccoef(Iofcc(Is)+Ic)
         Enddo
         Do Ic = 1, Nc
            Icb(Ncsum+Ic) = Iccbeg(Iofcb(Is)+Ic)
            Ice(Ncsum+Ic) = Iccend(Iofcb(Is)+Ic)
         Enddo
         Nalp  = Nalp  + Np
         Ncc   = Ncc   + Np*Nc
         Ncsum = Ncsum + Nc
      Enddo

      Return
      End
