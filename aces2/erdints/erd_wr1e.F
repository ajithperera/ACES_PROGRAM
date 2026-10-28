      Subroutine Erd_wr1e(Lu,Label,B,Buf,Ibuf)
C
C Write a symmetry blocked one-electron SO matrix to unit Lu in VMOL's
C format (vmol/onel.F): a '********' record, the label, then records
C (Buf(600), Ibuf(600), Nut) and a final record with Nut = -1. The
C index is the lower triangle of each irrep block, irrep by irrep:
C offset(irrep) + i(i-1)/2 + j, i >= j counted within the irrep.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Parameter (Lnbuf = 600)
      Character*8 Label, Stars
      Dimension B(Nso,Nso), Buf(Lnbuf)
      Integer Ibuf(Lnbuf)
      Data Stars /'********'/

      Write(Lu) Stars
      Write(Lu) Label

      Nut  = 0
      Ipk0 = 0
      Iso0 = 0
      Do Ig = 1, Nirr
         Nb = Nbfirr(Ig)
         Do I = 1, Nb
            Do J = 1, I
               Nut = Nut + 1
               Buf(Nut)  = B(Iso0+I,Iso0+J)
               Ibuf(Nut) = Ipk0 + I*(I-1)/2 + J
               If (Nut .Eq. Lnbuf) Then
                  Write(Lu) Buf, Ibuf, Nut
                  Nut = 0
               Endif
            Enddo
         Enddo
         Ipk0 = Ipk0 + Nb*(Nb+1)/2
         Iso0 = Iso0 + Nb
      Enddo
      If (Nut .Ne. 0) Write(Lu) Buf, Ibuf, Nut
      Nut = -1
      Write(Lu) Buf, Ibuf, Nut

      Return
      End
