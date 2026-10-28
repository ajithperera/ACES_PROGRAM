










      Subroutine Scrnc_wc_ao2mo(Work,Maxcor,Jj_lrespn,Nao,Nao_pairs,
     +                          C,Nmo,Wc_mo,Nmo_pairs)
C
C Transforms scrnc's AO-pair screened-interaction matrix Jj_lrespn
C (Nao_pairs x Nao_pairs, Nao_pairs=Nao*Nao, FULL unpacked AO-pair
C basis) into the MO particle-hole basis needed for the GW self-energy
C contraction (see SCRNC_DEBUG_LOG.md, 2026-08-17 "self-energy
C architecture" entry):
C
C   Wc_mo(p,q,r,s) = Sum(Mu,Nu,Lam,Sig)
C        C(Mu,p)*C(Nu,q)*Jj_lrespn(Iao_pair,Jao_pair)*C(Lam,r)*C(Sig,s)
C
C using the AO-pair index convention Iao_pair=(Mu-1)*Nao+Nu (outer
C loop = Mu, inner loop = Nu), confirmed authoritative via
C scrnc_prep_2ints.F:68-77 (the Irrep_ao_pairs construction shared
C throughout scrnc). Output uses the identical convention on the MO
C side: Ipq=(p-1)*Nmo+q, so Wc_mo is Nmo_pairs x Nmo_pairs with
C Nmo_pairs=Nmo*Nmo.
C
C Done as two block transforms (avoids ever forming the full
C Nao_pairs x Nao_pairs Kronecker matrix C(x)C explicitly): first
C transform the row AO-pair (Mu,Nu)->(p,q) for every column AO-pair,
C giving an intermediate Wtmp(pq,Lam,Sig); then transform the column
C AO-pair (Lam,Sig)->(r,s) for every (already row-transformed) MO-pair.
C Deliberately written as plain small (Nao x Nao) / (Nmo x Nmo) block
C Xgemm calls, not a single large Nao_pairs-scale contraction --
C correctness-first for the STO-3G-scale systems this whole project has
C been verified against so far. Storage is O(Nao**4); NOT intended to
C scale to a realistic basis set as-is (flagged, not yet addressed,
C same debug-log entry). Standard ACES2 Work-array scratch convention
C (no ALLOCATABLE, consistent with the rest of this module). All
C loop/index variable names are kept in the I-N implicit-integer range
C (never P,Q,R,S,Sig, which fall in the O-Z implicit-DOUBLE range).
C
C ROOT CAUSE FOUND AND FIXED, 2026-08-17 (see SCRNC_DEBUG_LOG.md,
C "continued 8"): Jj_lrespn's name starts with 'J', which IS in the
C implicit I-N range -- and this file's own
C "Implicit Double Precision (A-H,O-Z)" does NOT cover I-N, so a bare
C "Dimension Jj_lrespn(...)" (no explicit type) silently made it
C Integer*8, not Double Precision. Confirmed directly via gdb: the
C address computed for Jj_lrespn was correct at every level of the call
C chain the whole time (ruling out the aliasing/overlap theory chased
C earlier), but every READ of an element converted its raw IEEE-754 bit
C pattern as if it were a 64-bit integer (matching the exact ~1e18-1e19
C garbage magnitude seen throughout this investigation -- e.g. a normal
C value like 0.7776 has an int64 bit-reinterpretation of
C ~4.6D+18). Fixed by declaring
C "Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)" explicitly in every
C routine in this file (and in Scrnc_wc_ao2mo_test.F) instead of relying
C on Dimension alone -- matching the rest of the codebase's own
C established convention for this exact array (see
C scrnc_jj_rcc_assemble.F:6, "Double Precision JJ_lrespn(...)").
C
C Given that, the "VERIFIED CORRECT IN ISOLATION" standalone tests this
C comment used to cite were never actually meaningful (their test
C harness happened to declare Jj_lrespn correctly, so they never
C exercised this bug) -- the fix above is the real verification.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Maxcor)
      Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)
      Dimension C(Nao,Nmo)
      Dimension Wc_mo(Nmo_pairs,Nmo_pairs)
C
      I_wtmp    = 1
      I_blk_ao  = I_wtmp    + Nmo_pairs*Nao_pairs
      I_blk_half= I_blk_ao  + Nao*Nao
      I_blk_mo  = I_blk_half+ Nao*Nmo
      Iend      = I_blk_mo  + Nmo*Nmo
C
      If (Iend .GT. Maxcor) Call Insmem("Scrnc_wc_ao2mo",Iend,Maxcor)
C
      Call Scrnc_wc_ao2mo_do(Jj_lrespn,Nao,Nao_pairs,C,Nmo,
     +                       Wc_mo,Nmo_pairs,Work(I_wtmp),
     +                       Work(I_blk_ao),Work(I_blk_half),
     +                       Work(I_blk_mo))
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_wc_ao2mo_do(Jj_lrespn,Nao,Nao_pairs,C,Nmo,
     +                             Wc_mo,Nmo_pairs,Wtmp,Blk_ao,
     +                             Blk_half,Blk_mo)
C
      Implicit Double Precision (A-H,O-Z)
      Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)
      Dimension C(Nao,Nmo)
      Dimension Wc_mo(Nmo_pairs,Nmo_pairs)
      Dimension Wtmp(Nmo_pairs,Nao_pairs)
      Dimension Blk_ao(Nao,Nao), Blk_half(Nao,Nmo), Blk_mo(Nmo,Nmo)
C
C --- Step 1: row transform (Mu,Nu) -> (Ip,Iq), one Nao x Nao block
C     per column AO-pair (Lam,Nsig=Jao_pair).
C
      Do Jao_pair = 1, Nao_pairs
C
         Do Mu = 1, Nao
            Do Nu = 1, Nao
               Iao_pair = (Mu-1)*Nao + Nu
               Blk_ao(Mu,Nu) = Jj_lrespn(Iao_pair,Jao_pair)
            Enddo
         Enddo
C
C        Blk_half(Nu,Ip) = Sum_Mu Blk_ao(Mu,Nu)*C(Mu,Ip) [=Blk_ao^T*C]
C
         Call Xgemm("T","N",Nao,Nmo,Nao,1.0D0,Blk_ao,Nao,
     +              C,Nao,0.0D0,Blk_half,Nao)
C
C        Blk_mo(Iq,Ip) = Sum_Nu C(Nu,Iq)*Blk_half(Nu,Ip) [=C^T*Blk_half]
C
         Call Xgemm("T","N",Nmo,Nmo,Nao,1.0D0,C,Nao,
     +              Blk_half,Nao,0.0D0,Blk_mo,Nmo)
C
         Do Ip = 1, Nmo
            Do Iq = 1, Nmo
               Ipq = (Ip-1)*Nmo + Iq
               Wtmp(Ipq,Jao_pair) = Blk_mo(Iq,Ip)
            Enddo
         Enddo
C
      Enddo
C
C --- Step 2: column transform (Lam,Nsig) -> (Ir,Is), one Nao x Nao
C     block per row MO-pair (Ip,Iq=Ipq), reading from Wtmp above.
C
      Do Ipq = 1, Nmo_pairs
C
         Do Lam = 1, Nao
            Do Nsig = 1, Nao
               Jao_pair = (Lam-1)*Nao + Nsig
               Blk_ao(Lam,Nsig) = Wtmp(Ipq,Jao_pair)
            Enddo
         Enddo
C
         Call Xgemm("T","N",Nao,Nmo,Nao,1.0D0,Blk_ao,Nao,
     +              C,Nao,0.0D0,Blk_half,Nao)
C
         Call Xgemm("T","N",Nmo,Nmo,Nao,1.0D0,C,Nao,
     +              Blk_half,Nao,0.0D0,Blk_mo,Nmo)
C
         Do Ir = 1, Nmo
            Do Is = 1, Nmo
               Irs = (Ir-1)*Nmo + Is
               Wc_mo(Ipq,Irs) = Blk_mo(Is,Ir)
            Enddo
         Enddo
C
      Enddo
C
      Return
      End
