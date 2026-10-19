      Subroutine Pccd_rfo_disp(C0,X,Fac,F,U,Scr,G,Irp,Nbas,Npar)
C
C Write C = C0*exp(K) to SCFEVCA0 and SCFEVCB0 (RHF), where
C K(p,q) = Fac*X(pq) = -K(q,p) for the pairs p<q of the same irrep
C (Irp) in row order.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension C0(Nbas,Nbas),X(Npar),F(Nbas,Nbas),U(Nbas,Nbas)
      Dimension Scr(Nbas,Nbas),G(Nbas,Nbas),Irp(Nbas)
#include "machsp.com"
      Data Zero,One /0.0D0,1.0D0/

      Call Dzero(F,Nbas*Nbas)
      K = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               F(Ip,Iq) =  Fac*X(K)
               F(Iq,Ip) = -Fac*X(K)
            Endif
         Enddo
      Enddo
      Call Pccd_rfo_expm(F,U,Scr,G,Nbas)
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,C0,Nbas,U,Nbas,Zero,
     +           Scr,Nbas)
      Call Putrec(20,"JOBARC","SCFEVCA0",Nbas*Nbas*Iintfp,Scr)
      Call Putrec(20,"JOBARC","SCFEVCB0",Nbas*Nbas*Iintfp,Scr)
      Return
      End
