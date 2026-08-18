C Standalone external orchestration driver for the contour-deformation
C Sigma_c(p,q; i*omega) self-energy integral (see SCRNC_DEBUG_LOG.md,
C 2026-08-17 "self-energy architecture" / "continued 9" entries):
C
C   Sigma_c(p,q; i*omega) = -(1/2*pi) * Integral d(omega') Sum_m
C        [1/(i*omega+i*omega'-eps_m)] * W_c(p,m,q,m; i*omega')
C
C Design (Ajith, 2026-08-17; restructured 2026-08-18 into a static/
C per-frequency split -- see SCRNC_DEBUG_LOG.md and
C project_scrnc_pipeline_restructure.md memory): external orchestration
C (not an internal Fortran loop inside Scrnc_driver.F), in Fortran (not
C a shell/python wrapper), reusing a SINGLE shared working directory
C across all N frequency-point pipeline runs. This program is
C deliberately NOT part of the normal ACES2 aces_init/JOBARC-linked
C build -- it never calls Getrec/Putrec itself; instead it spawns the
C existing, unmodified ACES2 binaries as subprocesses via the standard
C Fortran 2008 EXECUTE_COMMAND_LINE intrinsic, and reads back the plain
C text file "wcreduc.dat" that Scrnc_sigma_prep.F (compiled into xscrnc
C itself, normal ACES2 build) writes at the end of each xscrnc run.
C
C Two-phase pipeline, traced+confirmed against source 2026-08-18: only
C xscrnc reads the frequency file (scrnc_set_vars.F/vcceh.F) -- xjoda,
C xvmol, xvmol2ja, xvprops, xvscf, xvtran, xintprc, and the ground-state
C amplitude step (xvcc,xlambda for CCSD; xringcc,xringcl for ring-CC --
C confirmed zero frequency/omega references in ringcc/ringcl source,
C and rpa/'s superficially-similar -D_RCC-gated scrnc_set_vars.F copy
C builds into a wholly separate xrpa binary, never linked into
C xringcc/xringcl) are all frequency-independent. So:
C   - STATIC_CMDS runs ONCE, before the loop, to build the frequency-
C     independent state (integrals/SCF/ground-state amplitudes).
C   - PIPELINE_CMDS (now just the frequency-dependent step, e.g.
C     "xscrnc") runs once per quadrature node inside the loop.
C Naively looping PIPELINE_CMDS alone on top of the SAME on-disk static
C state repeatedly is NOT safe: xscrnc itself writes dynamic list
C records (e.g. into MOABCD, on the CCSD linear-response path) that
C collide on a second same-directory invocation -- the exact
C "@GETLST: Assertion failed" failure mode already hit and documented
C below for the old always-rerun-everything design, just one level
C removed (there it was xjoda not fully resetting xscrnc's leftover
C dynamic records; here it would be a second bare xscrnc call hitting
C its OWN previous invocation's leftover records with no xjoda in
C between at all). Fix: snapshot the directory's full generated-file
C state immediately after STATIC_CMDS completes (".static_snapshot/"),
C then each loop iteration clears the directory back to that pristine
C snapshot (never to bare ZMAT/GENBAS) before running PIPELINE_CMDS --
C so every xscrnc invocation always starts from the SAME clean static
C state, never from another xscrnc call's leftovers.
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
C                  Omega_ext (external evaluation grid) is NOT read
C                  from here -- it is generated automatically once
C                  Eps(1..Nmo) is known (end of node 1, see "Auto-
C                  build" below), so the grid always reaches this
C                  system's own orbital-energy range regardless of
C                  molecule/basis. See SCRNC_DEBUG_LOG.md 2026-08-18
C                  "continued 25".
C   static_cmds:   single line, the shell command that builds the
C                  frequency-INDEPENDENT state ONCE, before the loop,
C                  e.g.
C                  "xjoda;xvmol;xvmol2ja;xvprops;xvscf;xvtran;xintprc;
C                   xringcc;xringcl" (ring-CC path) or
C                  "xjoda;xvmol;xvmol2ja;xvprops;xvscf;xvtran;xintprc;
C                   xvcc;xlambda" (CCSD path).
C   pipeline_cmds: single line, the shell command for the frequency-
C                  DEPENDENT step only, run once per quadrature node,
C                  e.g. "xscrnc" -- deliberately not hardcoded here so
C                  this driver works for either path (or any future
C                  one) without editing Fortran source.
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
      Character*500 Cmdline, Staticcmd
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
      Close(Unit=21)
      If (N .GT. Maxfreq) Then
         Write(6,"(a,I5)") " @-Scrnc_sigma_driver: N exceeds Maxfreq, ",
     +                     Maxfreq
         Stop 1
      Endif
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
      Inquire(File="static_cmds",Exist=Fexist)
      If (.NOT. Fexist) Then
         Write(6,"(a)")
     +      " @-Scrnc_sigma_driver: static_cmds not found"
         Stop 1
      Endif
      Open(Unit=25,File="static_cmds",Form="Formatted",Status="OLD")
      Read(25,"(a)") Staticcmd
      Close(Unit=25)
C
      Write(6,"(a,I4,a,1PD12.4)")
     +   " @-Scrnc_sigma_driver: N=",N," Omega_c=",Omega_c
      Write(6,"(a,a)")  "   static (once):     ",Trim(Staticcmd)
      Write(6,"(a,a)")  "   per-frequency-node: ",Trim(Cmdline)
C
C --- Build the imaginary-frequency quadrature grid.
C
      Call Scrnc_sigma_freq_grid(N,Omega_c,Xnode,Wnode)
C
      Nmo = 0
      Nomega = 0
C
C --- Phase 1: build the frequency-independent state ONCE, then
C     snapshot it so every per-node xscrnc call below can restore a
C     pristine copy instead of running the whole chain again (see
C     header comment for why a bare restore-to-ZMAT/GENBAS-and-rerun-
C     everything is being replaced, not just the "static_cmds" part).
C
      Call Execute_Command_Line(
     +   "find . -maxdepth 1 -type f "//
     +   "! -name ZMAT ! -name GENBAS ! -name sigma_input "//
     +   "! -name pipeline_cmds ! -name static_cmds ! -name '*.log' "//
     +   "-delete",
     +   Wait=.True.)
C
      Write(6,"(a)") " @-Scrnc_sigma_driver: running static prefix"
      Call Execute_Command_Line(Trim(Staticcmd),Wait=.True.,
     +                           Exitstat=Iexit)
      If (Iexit .NE. 0) Then
         Write(6,"(a,I6)")
     +      " @-Scrnc_sigma_driver: static_cmds FAILED, exit status ",
     +      Iexit
         Stop 1
      Endif
C
      Call Execute_Command_Line(
     +   "rm -rf .static_snapshot && mkdir .static_snapshot && "//
     +   "cp -p $(find . -maxdepth 1 -type f "//
     +   "! -name ZMAT ! -name GENBAS ! -name sigma_input "//
     +   "! -name pipeline_cmds ! -name static_cmds "//
     +   "! -name '*.log') .static_snapshot/",
     +   Wait=.True.,Exitstat=Iexit)
      If (Iexit .NE. 0) Then
         Write(6,"(a,I6)")
     +      " @-Scrnc_sigma_driver: static-state snapshot FAILED, "//
     +      "exit status ",Iexit
         Stop 1
      Endif
C
C --- Phase 2: one xscrnc (frequency-dependent) subprocess run per
C     quadrature node, always starting from a fresh restore of the
C     Phase-1 snapshot -- never from another node's leftover state.
C
      Do K = 1, N
C
         Call Execute_Command_Line(
     +      "find . -maxdepth 1 -type f "//
     +      "! -name ZMAT ! -name GENBAS ! -name sigma_input "//
     +      "! -name pipeline_cmds ! -name static_cmds "//
     +      "! -name '*.log' -delete && "//
     +      "cp -p .static_snapshot/* .",
     +      Wait=.True.,Exitstat=Iexit)
         If (Iexit .NE. 0) Then
            Write(6,"(a,I4,a,I6)")
     +         " @-Scrnc_sigma_driver: snapshot restore FAILED at "//
     +         "node ",K,", exit status ",Iexit
            Stop 1
         Endif
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
         If (K .EQ. 1) Then
C
C           Auto-build the external evaluation grid, now that this
C           system's own Eps(1..Nmo) is known -- universal across any
C           molecule/basis, not tuned per test case. See
C           SCRNC_DEBUG_LOG.md 2026-08-18 "continued 25": a fixed
C           valence-scale sub-grid (proven across every case run so
C           far in this project -- CC-response energy scales in the
C           frontier-orbital region don't vary wildly with system)
C           plus a geometric extension (x1.6/step) computed from THIS
C           run's own Max|Eps|, reaching Safety*Emax regardless of how
C           deep or shallow that turns out to be. Cheap either way --
C           Sigma_c at any Omega_ext point reuses the same Wc_reduced,
C           no extra pipeline cost, so there is no reason to under-
C           reach.
C
            Nomega = 6
            Omega_ext(1) = 0.0D0
            Omega_ext(2) = 0.3D0
            Omega_ext(3) = 0.6D0
            Omega_ext(4) = 0.9D0
            Omega_ext(5) = 1.2D0
            Omega_ext(6) = 1.5D0
C
            Emax = 0.0D0
            Do Im = 1, Nmo
               Emax = Max(Emax,Abs(Eps(Im)))
            Enddo
            Etarget = 1.3D0 * Emax
            Omega_last = 1.5D0
C
            Do While (Omega_last.LT.Etarget .AND. Nomega.LT.Maxomega)
               Omega_last = Omega_last * 1.6D0
               Nomega = Nomega + 1
               Omega_ext(Nomega) = Omega_last
            Enddo
C
            Write(6,"(a,I4,a,1PD12.4,a,1PD12.4)")
     +         " @-Scrnc_sigma_driver: auto Omega_ext grid, Nomega=",
     +         Nomega," Emax=",Emax," reaching ",Omega_last
C
C           Zero the accumulators -- deferred to here since Nomega
C           wasn't known before this point.
C
            Do J = 1, Nomega
             Do Iq = 1, Maxmo
              Do Ip = 1, Maxmo
               Sigma_re(Ip,Iq,J) = 0.0D0
               Sigma_im(Ip,Iq,J) = 0.0D0
              Enddo
             Enddo
            Enddo
C
         Endif
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
C     solved by bracketed bisection on F(z)=z-eps_HF-Re[Sigma_c(z)]
C     (see Scrnc_sigma_qp_bracket below for why, not naive fixed-point
C     iteration).
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
C Each MO builds its Thiele fit from only the Omega_ext points closest
C to its OWN |eps_HF|, not the full shared Omega_ext array -- see
C SCRNC_DEBUG_LOG.md 2026-08-18 "continued 23": swept Nomega=8..36 on
C the same shared grid and found eps_QP scatters by up to ~0.25 Ha
C with N, non-monotonically; confirmed via independent arbitrary-
C precision (mpmath, 15-100 digit) arithmetic that this is NOT
C double-precision roundoff (bit-identical root at every precision)
C but the well-documented intrinsic ill-posedness of high-order
C EXACT-interpolation Pade continuation -- more points does not mean
C a better fit. The original hand-picked 6-point [0,1.5] grid gave
C clean, pilot-verified results for every valence/virtual MO; only a
C deep-core orbital far outside that range ever broke. Fix: the auto-
C built union grid (see main program, right after node 1's Eps is
C known) is cheap regardless of size -- Sigma_re/Sigma_im at any
C Omega_ext reuses the same already-computed Wc_reduced, no extra
C pipeline cost -- but each MO selects only its own nearby subset,
C recovering the original 6-point valence/virtual behavior
C automatically while giving deep orbitals their own locally-relevant
C points instead of either a wild extrapolation or an unstable
C all-points fit.
C
C Starting count is Nsel0=6 points (static + 5 nearest to |eps_HF|);
C if Scrnc_sigma_qp_bracket can't find a clean, non-pole-adjacent
C bracket at that count, ESCALATE to 8, then 10, ... up to Nselmax,
C adding more nearby (still LOCAL, not far-away) points before giving
C up -- see SCRNC_DEBUG_LOG.md 2026-08-18 "continued 25"/"continued
C 26": 6 points was not uniformly sufficient (H2O/STO-3G's own core
C orbital needed more, even though DZ's near-identical-energy core
C orbital was fine at 6) and there was no way to predict which MOs
C would need it. This escalation is NOT the same failure mode as
C "continued 23"'s huge-shared-grid instability -- these extra points
C stay local to the MO's own energy, not spread across the whole
C union grid's range.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Omega_ext(Maxomega)
      Dimension Sigma_re(Maxmo,Maxmo,Maxomega)
      Dimension Sigma_im(Maxmo,Maxmo,Maxomega)
      Dimension Eps(Maxmo)
C
      Double Complex Zpt(Maxomega), Fpt(Maxomega), Acoef(Maxomega)
      Double Complex Fval
C
      Parameter (Nsel0=6, Nselstep=2, Nselmax=20)
      Dimension Isel(Nselmax), Dsel(Nselmax)
C
      Open(Unit=26,File="qp_energies.dat",Form="Formatted",
     +     Status="UNKNOWN")
      Write(26,"(a)") "  p      eps_HF          eps_QP     "//
     +                 "  Im[Sigma_c(eps_QP)]   iters  converged"
C
      Istatic = 1
      Do J = 2, Nomega
         If (Omega_ext(J) .LT. Omega_ext(Istatic)) Istatic = J
      Enddo
C
      Do Ip = 1, Nmo
C
         Anchor = Abs(Eps(Ip))
         Iconverged = 0
         Nseltry = Nsel0
C
         Do While (Iconverged.EQ.0 .AND. Nseltry.LE.Nselmax .AND.
     +             Nseltry.LE.Nomega)
C
            Nuse = Min(Nseltry,Nomega)
C
C           The static (Omega_ext closest to 0) point is ALWAYS kept
C           -- dropping it can destabilize the fit even when the
C           remaining points are all individually close to Anchor
C           (confirmed on H2O/DZ MO 2, "continued 23"). Insertion-
C           select the remaining Nuse-1 Omega_ext indices closest to
C           Anchor, keeping Isel/Dsel sorted by distance ascending.
C
            Isel(1) = Istatic
            Dsel(1) = Abs(Omega_ext(Istatic)-Anchor)
            Nfound = 1
            Do J = 1, Nomega
               If (J .EQ. Istatic) Goto 50
               Dist = Abs(Omega_ext(J)-Anchor)
               If (Nfound .LT. Nuse) Then
                  Nfound = Nfound + 1
                  Isel(Nfound) = J
                  Dsel(Nfound) = Dist
                  Islot = Nfound
               Else If (Dist .LT. Dsel(Nuse)) Then
                  Isel(Nuse) = J
                  Dsel(Nuse) = Dist
                  Islot = Nuse
               Else
                  Islot = 0
               Endif
               Do While (Islot .GT. 2)
                  If (Dsel(Islot) .GE. Dsel(Islot-1)) Goto 40
                  Dtmp = Dsel(Islot)
                  Itmp = Isel(Islot)
                  Dsel(Islot) = Dsel(Islot-1)
                  Isel(Islot) = Isel(Islot-1)
                  Dsel(Islot-1) = Dtmp
                  Isel(Islot-1) = Itmp
                  Islot = Islot - 1
               Enddo
   40          Continue
   50          Continue
            Enddo
C
C           Re-sort the Nuse selected indices by Omega_ext value
C           ascending (standard Thiele-fit node ordering).
C
            Do Ii = 1, Nuse-1
               Do Jj = Ii+1, Nuse
                  If (Omega_ext(Isel(Jj)).LT.Omega_ext(Isel(Ii))) Then
                     Itmp = Isel(Ii)
                     Isel(Ii) = Isel(Jj)
                     Isel(Jj) = Itmp
                  Endif
               Enddo
            Enddo
C
            Fscale = 0.0D0
            Do Ii = 1, Nuse
               Jsrc = Isel(Ii)
               Zpt(Ii) = Dcmplx(0.0D0,Omega_ext(Jsrc))
               Fpt(Ii) = Dcmplx(Sigma_re(Ip,Ip,Jsrc),
     +                          Sigma_im(Ip,Ip,Jsrc))
               Fscale = Max(Fscale,Cdabs(Fpt(Ii)))
            Enddo
C
            Call Scrnc_sigma_pade_coeff(Nuse,Zpt,Fpt,Acoef,Maxomega)
C
            Call Scrnc_sigma_qp_bracket(Nuse,Zpt,Acoef,Eps(Ip),Fscale,
     +                                    Maxomega,Eps_qp,Fval,
     +                                    Iconverged,Iter)
C
            If (Iconverged .EQ. 0) Nseltry = Nseltry + Nselstep
C
         Enddo
C
         If (Iconverged .EQ. 0) Then
            Write(6,"(a,i3,a,i3,a)")
     +         " @-Scrnc_sigma_qp_solve: WARNING -- no bracketed root"
     +         //" found for MO",Ip," even after escalating to",
     +         Min(Nseltry,Nselmax)," points"
         Else If (Nseltry .GT. Nsel0) Then
            Write(6,"(a,i3,a,i3,a)")
     +         " @-Scrnc_sigma_qp_solve: MO",Ip," needed escalation"//
     +         " to",Nuse," points to find a clean bracket"
         Endif
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
      Subroutine Scrnc_sigma_qp_bracket(N,Zpt,Acoef,Eps_hf,Fscale,
     +                                    Maxomega,Eps_qp,Fval,
     +                                    Iconverged,Iter)
C
C Finds the root of F(z)=z-Eps_hf-Re[Sigma_c(z)] closest to Eps_hf,
C by scanning an expanding window centered on Eps_hf for a genuine
C sign change of F (rejecting any interval where the Pade-evaluated
C Sigma_c itself blows up relative to Fscale -- the magnitude actually
C spanned by the fit's own sample data -- the signature of a nearby
C pole rather than a root), then bisecting the bracket to Tol.
C
C Fscale-relative rejection, not a fixed absolute cutoff: a blanket
C |Sigma_c|<1000 threshold (the original implementation) missed a
C real case -- H2O/STO-3G MO 1's 6-point fit had a root sitting near
C a pole "only" mild enough to push |Sigma_c| to ~0.68, still under
C 1000 but ~190x every actual sampled value (max ~0.0036) -- i.e.
C obviously spurious relative to the data, just not in an absolute
C sense. See SCRNC_DEBUG_LOG.md 2026-08-18 "continued 25".
C

C Why not naive fixed-point iteration eps_new=Eps_hf+Re[Sigma_c(eps)]
C (the original implementation): F(z)=0 can have a perfectly good,
C isolated root even where the fixed-point map's local slope has
C magnitude greater than 1 at that root -- a property of the actual
C Sigma_c(z) function being continued, not a numerics artifact.
C Fixed-point iteration is only guaranteed to converge locally when
C that slope is below 1 in magnitude, and provably cannot converge
C otherwise (settles into a limit cycle instead), even though the
C root itself is fine. Confirmed on the H2O/STO-3G MO-5 case (the
C original 51-iteration non-convergence): standalone check found a
C single clean root at z=-0.4968 (no nearby pole) with local slope
C -2.14, exactly reproducing the observed oscillation. See
C SCRNC_DEBUG_LOG.md 2026-08-18 for the full pilot-code analysis.
C
      Implicit Double Precision (A-H,O-Z)
      Double Complex Zpt(Maxomega), Acoef(Maxomega), Fval, Fv1, Fv2
C
      Parameter (Nscan=200, Maxdoubling=8, Tol=1.0D-8,
     +           Whalf0=2.0D0, Rejfac=20.0D0, Fmaxfloor=1.0D-2)
C
      Iconverged = 0
      Eps_qp = Eps_hf
      Fval = Dcmplx(0.0D0,0.0D0)
      Iter = 0
      Whalf = Whalf0
      Fmax = Max(Rejfac*Fscale,Fmaxfloor)
C
      Do Idbl = 1, Maxdoubling
C
         Dz = 2.0D0*Whalf/Nscan
         Z1 = Eps_hf - Whalf
         Call Scrnc_sigma_pade_eval(N,Zpt,Acoef,Dcmplx(Z1,0.0D0),
     +                                Fv1,Maxomega)
         F1 = Z1 - Eps_hf - Dble(Fv1)
         Bestdist = 1.0D30
         Zbrlo = 0.0D0
         Zbrhi = 0.0D0
         Ifound = 0
C
         Do Iscan = 1, Nscan
            Z2 = Z1 + Dz
            Call Scrnc_sigma_pade_eval(N,Zpt,Acoef,Dcmplx(Z2,0.0D0),
     +                                   Fv2,Maxomega)
            F2 = Z2 - Eps_hf - Dble(Fv2)
C
            If ((F1.LT.0.0D0) .NEQV. (F2.LT.0.0D0)) Then
               If (Cdabs(Fv1).LT.Fmax .AND.
     +             Cdabs(Fv2).LT.Fmax) Then
                  Dist = Min(Abs(Z1-Eps_hf),Abs(Z2-Eps_hf))
                  If (Dist .LT. Bestdist) Then
                     Bestdist = Dist
                     Zbrlo = Z1
                     Zbrhi = Z2
                     Ifound = 1
                  Endif
               Endif
            Endif
C
            Z1 = Z2
            F1 = F2
            Fv1 = Fv2
         Enddo
C
         If (Ifound .EQ. 1) Goto 200
C
         Whalf = Whalf * 2.0D0
C
      Enddo
C
      Return
C
  200 Continue
C
      Call Scrnc_sigma_pade_eval(N,Zpt,Acoef,Dcmplx(Zbrlo,0.0D0),
     +                             Fv1,Maxomega)
      Fone = Zbrlo - Eps_hf - Dble(Fv1)
C
      Do Iter = 1, 200
         Zmid = 0.5D0*(Zbrlo+Zbrhi)
         Call Scrnc_sigma_pade_eval(N,Zpt,Acoef,Dcmplx(Zmid,0.0D0),
     +                                Fval,Maxomega)
         Fmid = Zmid - Eps_hf - Dble(Fval)
         If ((Fmid.LT.0.0D0) .EQV. (Fone.LT.0.0D0)) Then
            Zbrlo = Zmid
            Fone = Fmid
         Else
            Zbrhi = Zmid
         Endif
         If (Abs(Zbrhi-Zbrlo) .LT. Tol) Goto 300
      Enddo
C
  300 Continue
C
      Eps_qp = 0.5D0*(Zbrlo+Zbrhi)
      Call Scrnc_sigma_pade_eval(N,Zpt,Acoef,Dcmplx(Eps_qp,0.0D0),
     +                             Fval,Maxomega)
      Iconverged = 1
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
