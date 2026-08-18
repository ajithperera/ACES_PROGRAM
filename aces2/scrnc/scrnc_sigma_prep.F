      Subroutine Scrnc_sigma_prep(Work,Maxcor,Iuhf,Jj_lrespn,Nao_pairs)
C
C Builds the correlation-only screened interaction W_c = W - v in the
C MO particle-hole basis, reduced to the diagonal-in-second/fourth-
C index slice Wc_reduced(p,q,m) = Wc_mo(p,m,q,m) that the contour-
C deformation self-energy actually needs (see SCRNC_DEBUG_LOG.md,
C 2026-08-17 "continued 5" entry, standard equation):
C
C   Sigma_c(p,q; i*omega) = -(1/2*pi) * Integral d(omega') Sum_m
C        [1/(i*omega+i*omega'-eps_m)] * W_c(p,m,q,m; i*omega')
C
C Writes Wc_reduced plus the alpha-spin SCF orbital energies (needed
C for eps_m above) to a plain formatted file "wcreduc.dat" so a
C separate, non-ACES-linked standalone driver program
C (scrnc_sigma_driver.F) can read it back after each single-frequency-
C point pipeline subprocess run completes and JOBARC is about to be
C overwritten by the next point's fresh xjoda -- see that file's own
C header comment for the orchestration this feeds.
C
C W (Jj_lrespn, AO-pair basis) already contains the bare Coulomb piece
C (W = v + v*chi*v, scrnc_driver.F's own description) -- v is
C subtracted here at the AO-pair level (Jj_corr = Jj_lrespn - Jj_bare,
C valid since the AO->MO transform is linear) before the existing,
C already-verified Scrnc_wc_ao2mo transform is reused unchanged. Jj_bare
C is built from Scrnc_get2ints's raw AO two-electron integrals (same
C routine scrnc_prep_2ints.F uses to prep List_xxxx) -- this returns a
C matrix that is ALREADY fully symmetrized under (Mu<->Nu) and
C (Lam<->Sig) index-pair swaps (confirmed by inspection: it writes all
C 8 permutationally-equivalent index combinations for every unique
C integral read off disk), so re-indexing it into scrnc's own
C Iao_pair=(Mu-1)*Nao+Nu convention (scrnc_prep_2ints.F:68-77) is safe
C regardless of Scrnc_get2ints's own internal indx2(i,j)=i+(j-1)*Nao
C convention -- the two conventions give identical values for bare
C integrals specifically because of that pre-existing symmetrization,
C unlike Jj_lrespn itself (a genuine response function, not bare-
C integral-symmetric) where getting this backwards would silently
C transpose the result (same failure class as prior bugs in this
C project, see SCRNC_DEBUG_LOG.md).
C
C Alpha-spin only (SCFEVECA/SCFEVALA) -- inherited scope, not a new
C decision: Scrnc_wc_ao2mo_test.F already made this same simplification
C for Wc_mo itself.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Maxcor)
      Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)
#include "machsp.com"
#include "info.com"
C
      Character*8 Labelc, Labele
      Data Labelc /"SCFEVECA"/
      Data Labele /"SCFEVALA"/
C
      Nao = Nint(Sqrt(Dble(Nao_pairs)))
      Nmo = Nocco(1) + Nvrto(1)
      Nmo_pairs = Nmo * Nmo
C
      I_c       = 1
      I_jbare   = I_c       + Nao*Nmo
      I_g2      = I_jbare   + Nao_pairs*Nao_pairs
      I_wcmo    = I_g2      + Nao_pairs*Nao_pairs
      I_wcred   = I_wcmo    + Nmo_pairs*Nmo_pairs
      I_eps     = I_wcred   + Nmo*Nmo*Nmo
      I_scr     = I_eps     + Nmo
      Iend      = I_scr     + Nmo_pairs*Nao_pairs + Nao*Nao + Nao*Nmo +
     +                        Nmo*Nmo
C
      If (Iend .GT. Maxcor) Call Insmem("Scrnc_sigma_prep",Iend,Maxcor)
C
      Call Getrec(20,"JOBARC",Labelc,Iintfp*Nao*Nmo,Work(I_c))
      Call Getrec(20,"JOBARC",Labele,Iintfp*Nmo,Work(I_eps))
C
C G2 (Work(I_g2)) is transient scratch, fully consumed inside
C Scrnc_sigma_bare_ao -- not reused after it returns, but kept as its
C own non-overlapping block (rather than reusing its space for I_wcmo/
C I_scr below) for simplicity; trivial extra cost at this project's
C STO-3G toy-system scale (same "correctness first, optimize storage
C later" convention already used elsewhere here, e.g. Scrnc_wc_ao2mo.F's
C own header comment).
C
      Call Scrnc_sigma_bare_ao(Work(I_jbare),Work(I_g2),Nao,Nao_pairs)
C
C Jj_corr = Jj_lrespn - Jj_bare, done in place in the Jj_bare buffer.
C
      Call Scrnc_sigma_corr_ao(Jj_lrespn,Work(I_jbare),Nao_pairs)
C
      Call Scrnc_wc_ao2mo(Work(I_scr),Maxcor-I_scr+1,
     +                     Work(I_jbare),Nao,Nao_pairs,
     +                     Work(I_c),Nmo,Work(I_wcmo),Nmo_pairs)
C
      Call Scrnc_sigma_reduce(Work(I_wcmo),Nmo,Nmo_pairs,
     +                        Work(I_wcred))
C
      Call Scrnc_sigma_write(Nmo,Work(I_wcred),Work(I_eps))
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_bare_ao(Jj_bare,G2,Nao,Nao_pairs)
C
C Builds the bare AO two-electron Coulomb matrix in scrnc's own
C Iao_pair=(Mu-1)*Nao+Nu convention, from Scrnc_get2ints's raw output
C (own convention indx2(i,j)=i+(j-1)*Nao -- see this file's own header
C comment for why re-indexing is safe here despite the mismatch). G2
C is caller-supplied Work-array scratch (standard ACES2 convention --
C not a local automatic array).
C
      Implicit Double Precision (A-H,O-Z)
      Double Precision Jj_bare(Nao_pairs,Nao_pairs)
      Double Precision G2(Nao_pairs,Nao_pairs)
      Double Precision Buf(600)
      Dimension Ibuf(600)
C
      Call Scrnc_get2ints(G2,Buf,Ibuf,Nao)
C
      Do Mu = 1, Nao
       Do Nu = 1, Nao
        Iao_pair = (Mu-1)*Nao + Nu
        I_get2   = Mu + (Nu-1)*Nao
        Do Lam = 1, Nao
         Do Sig = 1, Nao
          Jao_pair = (Lam-1)*Nao + Sig
          J_get2   = Lam + (Sig-1)*Nao
          Jj_bare(Iao_pair,Jao_pair) = G2(I_get2,J_get2)
         Enddo
        Enddo
       Enddo
      Enddo
C
C Standing correctness self-check (cheap, kept permanently -- same
C convention as Scrnc_wc_ao2mo_bfcheck): bare 2-electron AO integrals
C must be exactly symmetric under (Iao_pair<->Jao_pair) swap.
C
      Xmaxasym = 0.0D0
      Do Iao_pair = 1, Nao_pairs
       Do Jao_pair = 1, Nao_pairs
        Diff = Abs(Jj_bare(Iao_pair,Jao_pair) -
     +             Jj_bare(Jao_pair,Iao_pair))
        If (Diff .GT. Xmaxasym) Xmaxasym = Diff
       Enddo
      Enddo
      Write(6,"(a,1PD15.6)")
     +   " @-Scrnc_sigma_bare_ao, max bare-v asymmetry = ", Xmaxasym
C
C Diagnostic (2026-08-18, sym=on/off Iao_pair permutation investigation,
C see SCRNC_DEBUG_LOG.md "continued 18"): dump Jj_bare in its own
C natural (Mu-1)*Nao+Nu order, plain one-value-per-line (Iao_pair
C fastest), so its values can be matched, by value, against the
C packed-order bare block recoverable from Scrnc_jj_rcc_assemble.F's
C own diagnostic dumps -- an empirical, code-reading-independent way to
C read off the actual packed<->natural permutation rather than
C hand-deriving it. Avoids the "output" pretty-printer's column-
C blocked/zero-row-skipping format, which is unreliable to parse back.
#ifdef _DEBUG_LVL0
      Open(Unit=88,File="jjbare_natural.dat",Form="Formatted",
     +     Status="UNKNOWN")
      Write(88,*) Nao_pairs
      Do Jao_pair = 1, Nao_pairs
       Do Iao_pair = 1, Nao_pairs
        Write(88,*) Jj_bare(Iao_pair,Jao_pair)
       Enddo
      Enddo
      Close(Unit=88)
#endif
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_corr_ao(Jj_lrespn,Jj_bare,Nao_pairs)
C
C In-place: Jj_bare := Jj_lrespn - Jj_bare (becomes W_c = W - v).
C
      Implicit Double Precision (A-H,O-Z)
      Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)
      Double Precision Jj_bare(Nao_pairs,Nao_pairs)
C
      Do Jao_pair = 1, Nao_pairs
       Do Iao_pair = 1, Nao_pairs
        Jj_bare(Iao_pair,Jao_pair) = Jj_lrespn(Iao_pair,Jao_pair) -
     +                                Jj_bare(Iao_pair,Jao_pair)
       Enddo
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_reduce(Wc_mo,Nmo,Nmo_pairs,Wc_reduced)
C
C Wc_reduced(p,q,m) = Wc_mo(p,m,q,m), i.e. Wc_mo(Ipm,Iqm) with
C Ipm=(p-1)*Nmo+m, Iqm=(q-1)*Nmo+m -- same compound-index convention
C Scrnc_wc_ao2mo itself uses (confirmed in that file's own header).
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Wc_mo(Nmo_pairs,Nmo_pairs)
      Dimension Wc_reduced(Nmo,Nmo,Nmo)
C
      Do Im = 1, Nmo
       Do Iq = 1, Nmo
        Iqm = (Iq-1)*Nmo + Im
        Do Ip = 1, Nmo
         Ipm = (Ip-1)*Nmo + Im
         Wc_reduced(Ip,Iq,Im) = Wc_mo(Ipm,Iqm)
        Enddo
       Enddo
      Enddo
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_sigma_write(Nmo,Wc_reduced,Eps)
C
C Plain formatted dump, one value per line -- read back by
C scrnc_sigma_driver.F, a separate standalone program with no ACES2
C library linkage (so no JOBARC/Getrec access), between subprocess
C pipeline runs. Loop order (Ip fastest, then Iq, then Im) must match
C exactly on the read side.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Wc_reduced(Nmo,Nmo,Nmo)
      Dimension Eps(Nmo)
C
      Open(Unit=87,File="wcreduc.dat",Form="Formatted",
     +     Status="UNKNOWN")
      Write(87,*) Nmo
      Do Im = 1, Nmo
       Do Iq = 1, Nmo
        Do Ip = 1, Nmo
         Write(87,*) Wc_reduced(Ip,Iq,Im)
        Enddo
       Enddo
      Enddo
      Do Im = 1, Nmo
       Write(87,*) Eps(Im)
      Enddo
      Close(Unit=87)
C
      Return
      End
