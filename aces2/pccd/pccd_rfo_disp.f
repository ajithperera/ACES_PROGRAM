










      Subroutine Pccd_rfo_disp(C0,X,Fac,F,U,Scr,G,Irp,Nbas,Npar)
C
C Write C = C0*exp(K) to SCFEVCA0 and SCFEVCB0 (RHF), where
C K(p,q) = Fac*X(pq) = -K(q,p) for the pairs p<q of the same irrep
C (Irp) in row order.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension C0(Nbas,Nbas),X(Npar),F(Nbas,Nbas),U(Nbas,Nbas)
      Dimension Scr(Nbas,Nbas),G(Nbas,Nbas),Irp(Nbas)


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



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
