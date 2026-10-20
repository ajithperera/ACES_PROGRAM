










      Subroutine Pccd_rfo_prb(Mode,Iact,C,E,Gv,X,Ev,Nroot,Step,F,U,Scr,
     +                        G,Irp,Nbas,Npar,Ifail)
C
C OO_METHOD=RFO soft-mode probe at a verified minimum. Very flat
C directions (exact Hessian eigenvalues between Xzero and Xsoft, the
C size of the negative-eigenvalue threshold of the stability check)
C can hide a lower minimum close by: N2 6-31G SYMMETRY=OFF from the HF
C orbitals, exact-Hessian steps stop at -108.970836 (physical
C eigenvalues 4.9e-7, 5.2e-5), 6.3 mEh above the minimum the SR1 path
C reaches. Along the 5.2e-5 mode E - E0 = -2.7e-6 at 0.2, -9.3e-5 at
C 0.5 (then -108.977154 in 19 gradients); flat along the other.
C The energy is evaluated at +-Dprb along the (up to Mxprb) softest
C modes, one macro iteration each, from the stationary point C0. The
C first probe lower than E0 - Eprb becomes the new starting point.
C   Mode "START": X(Npar,Nroot), Ev(Nroot) the lowest exact eigenpairs
C     (Pccd_rsp_stab) at the stationary point C, E, Gv. Iact = 1: no
C     soft mode, nothing done; 0: first probe set up (C displaced).
C   Mode "CONT ": E (Ifail = 0) at the probe in C. Iact = 0: next
C     probe set up; 1: done, C, E, Gv are those of the stationary
C     point again; 2: this probe is lower, continue from it.
C State in JOBARC (RFO_PHAS = 2 while probing, RFO_PB* records).
C Debug: file pccd_rfo_prb.in, Dprb (and Xsoft) instead of the defaults.
C
      Implicit Double Precision(A-H,O-Z)
      Character*5 Mode
      Logical Exist
      Dimension C(Nbas*Nbas),Gv(Npar),X(Npar,*),Ev(*),Step(Npar)
      Dimension F(*),U(*),Scr(*),G(*),Irp(Nbas),Pb(3)


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



      Parameter (Xzero=1.0D-8, Mxprb=2, Eprb=1.0D-5)
      Data Zero,One /0.0D0,1.0D0/

      Nsq   = Nbas*Nbas
      Dprb  = 0.5D0
      Xsoft = 1.0D-4
      Inquire(File="pccd_rfo_prb.in",Exist=Exist)
      If (Exist) Then
         Open(Unit=44,File="pccd_rfo_prb.in",Status="Old")
         Read(44,*,End=5,Err=5) Dprb,Xsoft
 5       Close(44)
      Endif

      If (Mode .Eq. "START") Then
         Iact = 1
         Nprb = 0
         Do Ir = 1, Nroot
            If (Dabs(Ev(Ir)) .Ge. Xzero .And. Ev(Ir) .Lt. Xsoft .And.
     +          Nprb .Lt. Mxprb) Then
               Nprb = Nprb + 1
               If (Nprb .Ne. Ir) Call Dcopy(Npar,X(1,Ir),1,X(1,Nprb),1)
               Ev(Nprb) = Ev(Ir)
            Endif
         Enddo
         If (Nprb .Eq. 0) Return
         Write(6,"(a,i2,a,F5.2,a,1p,2E11.3)") " @-RFO",Nprb,
     +      " soft mode(s): probing the energy at +-",Dprb,
     +      "; eigenvalue(s)",(Ev(Ir),Ir=1,Nprb)
         Pb(1) = E
         Pb(2) = Dble(Nprb)
         Pb(3) = Zero
         Call Putrec(20,"JOBARC","RFO_PBST",3*Iintfp,Pb)
         Call Putrec(20,"JOBARC","RFO_PBC0",Nsq*Iintfp,C)
         Call Putrec(20,"JOBARC","RFO_PBG0",Npar*Iintfp,Gv)
         Call Putrec(20,"JOBARC","RFO_PBVC",Npar*Nprb*Iintfp,X)
         Call Putrec(20,"JOBARC","RFO_PHAS",1,2)
         Kprb = 0
      Else
         Call Getrec(20,"JOBARC","RFO_PBST",3*Iintfp,Pb)
         E0   = Pb(1)
         Nprb = Nint(Pb(2))
         Kprb = Nint(Pb(3))
         If (Ifail .Eq. 0) Then
            Write(6,"(a,i2,a,F18.12,a,1p,E10.3)") " @-RFO probe",Kprb,
     +         ": E =",E,"  E - E0 =",E-E0
            If (E .Lt. E0-Eprb) Then
               Write(6,"(2a)") " @-RFO lower energy along a soft mode:",
     +              " the optimization continues from this point."
               Call Putrec(20,"JOBARC","RFO_PHAS",1,0)
               Iact = 2
               Return
            Endif
         Else
            Write(6,"(a,i2,a)") " @-RFO probe",Kprb,
     +         ": pCCD amplitudes unconverged or unphysical."
         Endif
         If (Kprb .Ge. 2*Nprb) Then
C All probes higher: back to the stationary point.
            Call Getrec(20,"JOBARC","RFO_PBC0",Nsq*Iintfp,C)
            Call Getrec(20,"JOBARC","RFO_PBG0",Npar*Iintfp,Gv)
            E = E0
            Call Putrec(20,"JOBARC","RFO_PHAS",1,0)
            Write(6,"(2a)") " @-RFO no lower energy along the soft",
     +           " modes."
            Iact = 1
            Return
         Endif
      Endif

C Next probe: mode (Kprb+2)/2, sign + for odd Kprb+1.
      Kprb = Kprb + 1
      Imod = (Kprb+1)/2
      Sgn  = One
      If (Mod(Kprb,2) .Eq. 0) Sgn = -One
      Call Getrec(20,"JOBARC","RFO_PBC0",Nsq*Iintfp,C)
      Call Getrec(20,"JOBARC","RFO_PBVC",Npar*Nprb*Iintfp,X)
      Vn = Dsqrt(Ddot(Npar,X(1,Imod),1,X(1,Imod),1))
      Do K = 1, Npar
         Step(K) = Sgn*Dprb*X(K,Imod)/Vn
      Enddo
C SCFEVCA0/B0 = C0*exp(K(Step)).
      Call Pccd_rfo_disp(C,Step,One,F,U,Scr,G,Irp,Nbas,Npar)
      Call Putrec(20,"JOBARC","GRADTEST",1,0)
      Pb(1) = E0
      If (Mode .Eq. "START") Pb(1) = E
      Pb(2) = Dble(Nprb)
      Pb(3) = Dble(Kprb)
      Call Putrec(20,"JOBARC","RFO_PBST",3*Iintfp,Pb)
      Iact = 0
      Return
      End
