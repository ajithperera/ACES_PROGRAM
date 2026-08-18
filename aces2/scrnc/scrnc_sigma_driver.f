










C Standalone external orchestration driver for the contour-deformation
C Sigma_c(p,q; i*omega) self-energy integral (see SCRNC_DEBUG_LOG.md,
C 2026-08-17 "self-energy architecture" / "continued 9" entries):
C
C   Sigma_c(p,q; i*omega) = -(1/2*pi) * Integral d(omega') Sum_m
C        [1/(i*omega+i*omega'-eps_m)] * W_c(p,m,q,m; i*omega')
C
C Design (Ajith, 2026-08-17): external orchestration (not an internal
C Fortran loop inside Scrnc_driver.F), in Fortran (not a shell/python
C wrapper), reusing a SINGLE shared working directory across all N
C frequency-point pipeline runs. This program is deliberately NOT part
C of the normal ACES2 aces_init/JOBARC-linked build -- it never calls
C Getrec/Putrec itself; instead it spawns the existing, unmodified
C xjoda...xscrnc pipeline as a subprocess once per Gauss-Legendre
C quadrature node (via the standard Fortran 2008 EXECUTE_COMMAND_LINE
C intrinsic) and reads back the plain text file "wcreduc.dat" that
C Scrnc_sigma_prep.F (compiled into xscrnc itself, normal ACES2 build)
C writes at the end of each such run, right before the NEXT iteration's
C fresh xjoda overwrites that directory's JOBARC (rerunning a pipeline
C step on top of stale JOBARC is NOT idempotent -- see debug log's
C "PATH trap" entry -- so this program always drives a full fresh
C xjoda...xscrnc chain per node, never a partial resume).
C
C Only N pipeline subprocess runs are needed in total, regardless of
C how many external evaluation points Nomega is set to: the expensive
C quantum-chemistry pipeline only depends on the INNER quadrature
C frequency omega' (the actual dynamical Freq the pipeline solves at),
C while the OUTER evaluation frequency omega only enters algebraically
C in the G0 denominator. So Wc_reduced(p,q,m) is read back once per
C inner node k and immediately folded into the running Sigma_c
C accumulation for every outer omega_j in one pass -- no need to store
C all N Wc_reduced_k simultaneously, and this same N-run cache
C architecture is ready to support a later real-axis Pade-continuation
C frequency scan (many omega_j) without any additional pipeline runs.
C
C Local copies of Scrnc_freq_grid/Scrnc_gaussl (originally
C scrnc_freq_grid.F) -- same multiply-duplicated-per-module convention
C already established throughout this codebase (see that file's own
C header comment) -- needed here because this program deliberately
C does NOT link against the rest of scrnc/ACES2's libraries.
C
C Build (standalone, no ACES2 library linkage -- see debug log for the
C exact command used):
C   ifx -i8 -o scrnc_sigma_driver scrnc_sigma_driver.F
C   cp scrnc_sigma_driver ../bin/xscrncsigma
C
C Inputs (plain text files in the shared working directory):
C   sigma_input:   line 1: N (quadrature order)
C                  line 2: Omega_c (quadrature scale, a.u.)
C                  line 3: Nomega (number of external omega points)
C                  next Nomega lines: Omega_ext(j) values (a.u.)
C   pipeline_cmds: single line, the full shell command to run ONE
C                  frequency point's pipeline pass in the current
C                  directory, e.g.
C                  "xjoda;xvmol;xvmol2ja;xvprops;xvscf;xvtran;xintprc;
C                   xringcc;xringcl;xscrnc" (ring-CC path) or
C                  "xjoda;xvmol;xvmol2ja;xvprops;xvscf;xvtran;xintprc;
C                   xvcc;xlambda;xscrnc" (CCSD path) -- deliberately not
C                  hardcoded here so this driver works for either path
C                  (or any future one) without editing Fortran source.
C
C Output: "sigma_c_output.dat", formatted, real+imaginary Sigma_c(p,q)
C for each requested Omega_ext(j).
C
      Program Scrnc_sigma_driver
C
      Implicit Double Precision (A-H,O-Z)
      Parameter (Maxfreq=200, Maxomega=50, Maxmo=64)
C
      Dimension Xnode(Maxfreq), Wnode(Maxfreq)
      Dimension Omega_ext(Maxomega)
      Dimension Wc_reduced(Maxmo,Maxmo,Maxmo)
      Dimension Eps(Maxmo)
      Dimension Sigma_re(Maxmo,Maxmo,Maxomega)
      Dimension Sigma_im(Maxmo,Maxmo,Maxomega)
C
      Character*500 Cmdline
      Logical Fexist
C
      Pi = 4.0D0 * Datan(1.0D0)
C
C --- Read control files.
C
      Inquire(File="sigma_input",Exist=Fexist)
      If (.NOT. Fexist) Then
         Write(6,"(a)") " @-Scrnc_sigma_driver: sigma_input not found"
         Stop 1
      Endif
      Open(Unit=21,File="sigma_input",Form="Formatted",Status="OLD")
      Read(21,*) N
      Read(21,*) Omega_c
      Read(21,*) Nomega
      If (N .GT. Maxfreq) Then
         Write(6,"(a,I5)") " @-Scrnc_sigma_driver: N exceeds Maxfreq, ",
     +                     Maxfreq
         Stop 1
      Endif
      If (Nomega .GT. Maxomega) Then
         Write(6,"(a,I5)")
     +      " @-Scrnc_sigma_driver: Nomega exceeds Maxomega, ",Maxomega
         Stop 1
      Endif
      Do J = 1, Nomega
         Read(21,*) Omega_ext(J)
      Enddo
      Close(Unit=21)
C
      Inquire(File="pipeline_cmds",Exist=Fexist)
      If (.NOT. Fexist) Then
         Write(6,"(a)")
     +      " @-Scrnc_sigma_driver: pipeline_cmds not found"
         Stop 1
      Endif
      Open(Unit=22,File="pipeline_cmds",Form="Formatted",Status="OLD")
      Read(22,"(a)") Cmdline
      Close(Unit=22)
C
      Write(6,"(a,I4,a,1PD12.4)")
     +   " @-Scrnc_sigma_driver: N=",N," Omega_c=",Omega_c
      Write(6,"(a,I4)") "   Nomega=",Nomega
      Write(6,"(a,a)")  "   pipeline: ",Trim(Cmdline)
C
C --- Build the imaginary-frequency quadrature grid.
C
      Call Scrnc_sigma_freq_grid(N,Omega_c,Xnode,Wnode)
C
C --- Zero the accumulators.
C
      Do J = 1, Nomega
       Do Iq = 1, Maxmo
        Do Ip = 1, Maxmo
         Sigma_re(Ip,Iq,J) = 0.0D0
         Sigma_im(Ip,Iq,J) = 0.0D0
        Enddo
       Enddo
      Enddo
      Nmo = 0
C
C --- Main loop: one pipeline subprocess run per quadrature node.
C
      Do K = 1, N
C
C        Clear the shared working directory back to ZMAT/GENBAS/this
C        driver's own control files before each node's pipeline run --
C        re-running xjoda alone (without deleting the previous node's
C        leftover generated files first) is NOT sufficient: hit
C        "@GETLST: Assertion failed. List (1,127) does not exist" on
C        the CCSD (non-ring) path's second node, xjoda apparently
C        doesn't fully reset every record-container file (e.g.
C        MOABCD) that the CCSD linear-response path writes dynamic
C        list records into between separate xscrnc invocations in the
C        same directory. Same root cause as the already-documented
C        "rerunning a pipeline step on top of stale JOBARC is NOT
C        idempotent" lesson (SCRNC_DEBUG_LOG.md) -- this just deletes
C        the leftover files explicitly rather than relying on xjoda to
C        do it. Stays within the single shared directory (never
C        creates a new one) -- only the generated-file CONTENTS get
C        cleared each iteration, not the directory itself.
C
         Call Execute_Command_Line(
     +      "find . -maxdepth 1 -type f "//
     +      "! -name ZMAT ! -name GENBAS ! -name sigma_input "//
     +      "! -name pipeline_cmds ! -name '*.log' -delete",
     +      Wait=.True.)
C
         Open(Unit=23,File="frequency",Form="Formatted",
     +        Status="UNKNOWN")
         Write(23,*) 1
         Write(23,*) 1
         Write(23,*) Xnode(K)
         Close(Unit=23)
C
         Write(6,"(a,I4,a,I4,a,1PD15.6)")
     +      " @-Scrnc_sigma_driver: node ",K," of ",N,
     +      ", omega'=",Xnode(K)
C
         Call Execute_Command_Line(Trim(Cmdline),Wait=.True.,
     +                              Exitstat=Iexit)
         If (Iexit .NE. 0) Then
            Write(6,"(a,I4,a,I6)")
     +         " @-Scrnc_sigma_driver: pipeline FAILED at node ",K,
     +         ", exit status ",Iexit
            Stop 1
         Endif
C
         Inquire(File="wcreduc.dat",Exist=Fexist)
         If (.NOT. Fexist) Then
            Write(6,"(a,I4)")
     +         " @-Scrnc_sigma_driver: wcreduc.dat missing after ",
     +         "node ",K
            Stop 1
         Endif
C
         Open(Unit=24,File="wcreduc.dat",Form="Formatted",
     +        Status="OLD")
         Read(24,*) Nmo_this
         If (K .EQ. 1) Then
            Nmo = Nmo_this
            If (Nmo .GT. Maxmo) Then
               Write(6,"(a,I5)")
     +            " @-Scrnc_sigma_driver: Nmo exceeds Maxmo, ",Maxmo
               Stop 1
            Endif
         Else If (Nmo_this .NE. Nmo) Then
            Write(6,"(a,I4)")
     +         " @-Scrnc_sigma_driver: Nmo changed between nodes at ",
     +         K
            Stop 1
         Endif
         Do Im = 1, Nmo
          Do Iq = 1, Nmo
           Do Ip = 1, Nmo
            Read(24,*) Wc_reduced(Ip,Iq,Im)
           Enddo
          Enddo
         Enddo
         Do Im = 1, Nmo
          Read(24,*) Eps(Im)
         Enddo
         Close(Unit=24)
C
C        Fold node K's contribution into Sigma_c(p,q;i*omega_j) for
C        every requested outer omega_j.
C
C        Full-axis convolution Int_{-inf}^{inf} dw' G0(iw+iw')Wc(iw'),
C        folded to [0,inf) using Wc(-iw')=Wc(iw') (even in the response
C        frequency, true for any real-Hamiltonian response function --
C        see SCRNC_DEBUG_LOG.md 2026-08-17 "continued 16"), gives TWO
C        propagator denominators per quadrature node, not one:
C
C        i*omega_j + i*omega'_K - eps_m = -eps_m + i*(omega_j+omega'_K)
C        i*omega_j - i*omega'_K - eps_m = -eps_m + i*(omega_j-omega'_K)
C
C        Both terms share the same Wc_reduced(.,.,Im) (that symmetry is
C        exactly what licenses the fold in the first place) and the
C        same Wnode(K) quadrature weight.
C
         Do J = 1, Nomega
          Dim  = Omega_ext(J) + Xnode(K)
          Dim2 = Omega_ext(J) - Xnode(K)
          Do Im = 1, Nmo
           Dre = -Eps(Im)
           Dmag2  = Dre*Dre + Dim*Dim
           Dmag22 = Dre*Dre + Dim2*Dim2
           Tre =  Dre / Dmag2   + Dre  / Dmag22
           Tim = -Dim / Dmag2   - Dim2 / Dmag22
           Factor = -Wnode(K) / (2.0D0*Pi)
           Do Iq = 1, Nmo
            Do Ip = 1, Nmo
             Sigma_re(Ip,Iq,J) = Sigma_re(Ip,Iq,J) +
     +                           Factor*Tre*Wc_reduced(Ip,Iq,Im)
             Sigma_im(Ip,Iq,J) = Sigma_im(Ip,Iq,J) +
     +                           Factor*Tim*Wc_reduced(Ip,Iq,Im)
            Enddo
           Enddo
          Enddo
         Enddo
C
      Enddo
C
C --- Write results.
C
      Open(Unit=25,File="sigma_c_output.dat",Form="Formatted",
     +     Status="UNKNOWN")
      Write(25,"(a,I4,a,I4)") "Nmo=",Nmo," Nomega=",Nomega
      Do J = 1, Nomega
         Write(25,"(a,I4,a,1PD15.6)")
     +      "Omega_ext(",J,") = ",Omega_ext(J)
         Write(25,"(a)") "  p   q         Re(Sigma_c)         "//
     +                    "Im(Sigma_c)"
         Do Iq = 1, Nmo
          Do Ip = 1, Nmo
           Write(25,"(2I4,2(1x,1PD20.10))") Ip,Iq,
     +           Sigma_re(Ip,Iq,J),Sigma_im(Ip,Iq,J)
          Enddo
         Enddo
      Enddo
      Close(Unit=25)
C
C --- Pade continuation (imaginary axis -> real axis) and on-shell
C     quasiparticle equation solve, diagonal elements only (standard
C     non-degenerate QP treatment -- off-diagonal Sigma_c not used
C     here). Needs Nomega>=2 imaginary-axis points to be meaningful;
C     skipped (matching this driver's original debug-mode behavior)
C     if Nomega=1.
C
C     HF-reference simplification (see this session's SCRNC_DEBUG_LOG.md
C     entry): for canonical HF orbitals, eps_HF(p) already equals
C     <p|h+J-Sigma_x|p> exactly (that IS the definition of "canonical" --
C     the Fock operator h+J-Sigma_x is diagonal in this basis). The GW
C     self-energy's static-exchange piece is the SAME Sigma_x operator,
C     so the usual DFT-style "Sigma_x - v_xc" double-counting correction
C     is EXACTLY ZERO for a HF reference (unlike DFT, where v_xc is not
C     exact exchange) -- standard result for HF-based G0W0, not new
C     physics, but flagged here since it means this driver deliberately
C     builds NO separate Sigma_x machinery. On-shell QP equation used:
C        eps_QP(p) = eps_HF(p) + Re[Sigma_c(p,p; eps_QP(p))]
C     solved by simple fixed-point iteration from eps_QP^(0)=eps_HF(p).
C
      If (Nomega .GE. 2) Call Scrnc_sigma_qp_solve(Nmo,Nomega,
     +                        Omega_ext,Sigma_re,Sigma_im,Eps,Maxmo,
     +                        Maxomega)
C
      Write(6,"(a)")
     +   " @-Scrnc_sigma_driver: done, see sigma_c_output.dat"
C
      Stop
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_qp_solve(Nmo,Nomega,Omega_ext,Sigma_re,
     +                                 Sigma_im,Eps,Maxmo,Maxomega)
C
C Diagonal Pade continuation + on-shell QP solve, one MO index at a
C time -- see the calling comment above for the physics.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Omega_ext(Maxomega)
      Dimension Sigma_re(Maxmo,Maxmo,Maxomega)
      Dimension Sigma_im(Maxmo,Maxmo,Maxomega)
      Dimension Eps(Maxmo)
C
      Double Complex Zpt(Maxomega), Fpt(Maxomega), Acoef(Maxomega)
      Double Complex Zeval, Fval
C
      Parameter (Maxiter=50, Tol=1.0D-8)
C
      Open(Unit=26,File="qp_energies.dat",Form="Formatted",
     +     Status="UNKNOWN")
      Write(26,"(a)") "  p      eps_HF          eps_QP     "//
     +                 "  Im[Sigma_c(eps_QP)]   iters  converged"
C
      Do Ip = 1, Nmo
C
         Do J = 1, Nomega
            Zpt(J) = Dcmplx(0.0D0,Omega_ext(J))
            Fpt(J) = Dcmplx(Sigma_re(Ip,Ip,J),Sigma_im(Ip,Ip,J))
         Enddo
C
         Call Scrnc_sigma_pade_coeff(Nomega,Zpt,Fpt,Acoef,Maxomega)
C
         Eps_qp   = Eps(Ip)
         Iconverged = 0
         Do Iter = 1, Maxiter
            Zeval = Dcmplx(Eps_qp,0.0D0)
            Call Scrnc_sigma_pade_eval(Nomega,Zpt,Acoef,Zeval,Fval,
     +                                  Maxomega)
            Eps_new = Eps(Ip) + Dble(Fval)
            If (Abs(Eps_new-Eps_qp) .LT. Tol) Then
               Eps_qp = Eps_new
               Iconverged = 1
               Goto 100
            Endif
            Eps_qp = Eps_new
         Enddo
  100    Continue
C
         Write(26,"(I3,2(1x,1PD15.6),1x,1PD15.6,3x,I3,5x,I3)")
     +      Ip,Eps(Ip),Eps_qp,Dimag(Fval),Iter,Iconverged
C
      Enddo
C
      Close(Unit=26)
C
      Write(6,"(a)")
     +   " @-Scrnc_sigma_driver: QP solve done, see qp_energies.dat"
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_pade_coeff(N,Zpt,Fpt,Acoef,Maxomega)
C
C Thiele's reciprocal-difference continued-fraction coefficients
C (standard Pade analytic-continuation algorithm, e.g. Vidberg-Serene
C 1977 -- textbook numerical-GW technique, not new to this project).
C g(1,i)=f(z_i); g(k,i)=(g(k-1,k-1)-g(k-1,i)) /
C ((z_i-z_(k-1))*g(k-1,i)) for k=2..N, i=k..N; a_k=g(k,k).
C
      Implicit Double Precision (A-H,O-Z)
      Double Complex Zpt(Maxomega), Fpt(Maxomega), Acoef(Maxomega)
      Double Complex G(Maxomega,Maxomega)
C
      Do I = 1, N
         G(1,I) = Fpt(I)
      Enddo
C
      Do K = 2, N
         Do I = K, N
            G(K,I) = (G(K-1,K-1) - G(K-1,I)) /
     +               ((Zpt(I)-Zpt(K-1)) * G(K-1,I))
         Enddo
      Enddo
C
      Do K = 1, N
         Acoef(K) = G(K,K)
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_pade_eval(N,Zpt,Acoef,Z,Fval,Maxomega)
C
C Evaluates the Thiele continued fraction at Z via the standard
C forward two-term recurrence (Vidberg-Serene): A(0)=0, A(1)=a(1),
C B(0)=1, B(1)=1; A(n)=A(n-1)+(z-z(n-1))*a(n)*A(n-2),
C B(n)=B(n-1)+(z-z(n-1))*a(n)*B(n-2); f(z)=A(N)/B(N).
C
      Implicit Double Precision (A-H,O-Z)
      Double Complex Zpt(Maxomega), Acoef(Maxomega), Z, Fval
      Double Complex A(0:Maxomega), B(0:Maxomega)
C
      A(0) = Dcmplx(0.0D0,0.0D0)
      A(1) = Acoef(1)
      B(0) = Dcmplx(1.0D0,0.0D0)
      B(1) = Dcmplx(1.0D0,0.0D0)
C
      Do N2 = 2, N
         A(N2) = A(N2-1) + (Z-Zpt(N2-1))*Acoef(N2)*A(N2-2)
         B(N2) = B(N2-1) + (Z-Zpt(N2-1))*Acoef(N2)*B(N2-2)
      Enddo
C
      Fval = A(N) / B(N)
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_freq_grid(N,Omega_c,Xnode,Wnode)
C
C Local copy of Scrnc_freq_grid (scrnc_freq_grid.F) -- see this file's
C own header comment for why this is duplicated rather than linked.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Xnode(N),Wnode(N)
      Dimension Xleg(N),Wleg(N)
C
      Call Scrnc_sigma_gaussl(-1.0D0,1.0D0,N,Xleg,Wleg)
C
      Do K = 1, N
         Xnode(K) = Omega_c * (1.0D0+Xleg(K)) / (1.0D0-Xleg(K))
         Wnode(K) = Wleg(K) * Omega_c * 2.0D0 /
     +              ((1.0D0-Xleg(K)) * (1.0D0-Xleg(K)))
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_gaussl(X1,X2,N,X,W)
C
C Local copy of aces2/vdint/gaussl.f, via scrnc_freq_grid.F's own
C Scrnc_gaussl.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension X(N),W(N)
      Parameter (Eps=3.D-14)
      M = (N+1)/2
      Xm = 0.5D0*(X2+X1)
      Xl = 0.5D0*(X2-X1)
      Do I = 1, M
         Z = Cos(3.141592654D0*(I-0.25D0)/(N+0.5D0))
    1    Continue
            P1 = 1.D0
            P2 = 0.D0
            Do J = 1, N
               P3 = P2
               P2 = P1
               P1 = ((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
            Enddo
            Pp = N*(Z*P1-P2)/(Z*Z-1.D0)
            Z1 = Z
            Z = Z1-P1/Pp
         If (Abs(Z-Z1) .GT. Eps) Goto 1
         X(I)     = Xm-Xl*Z
         X(N+1-I) = Xm+Xl*Z
         W(I)     = 2.D0*Xl/((1.D0-Z*Z)*Pp*Pp)
         W(N+1-I) = W(I)
      Enddo
      Return
      End
