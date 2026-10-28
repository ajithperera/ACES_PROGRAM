










      Subroutine Scrnc_jj_unpack_natural(Work,Maxcor,JJ_lrespn,
     +                                    Nao,Nao_pair)
C
C Un-packs JJ_lrespn from scrnc's internal irrep-grouped Iao_pair
C ordering (built by Scrnc_lrspn_rcc_main.F / Scrnc_lrspn_cc_main.F's
C own Npert_count running counter -- "visit all of irrep 1's AO-pairs,
C then irrep 2's, etc." -- and preserved by Scrnc_jj_rcc_assemble.F /
C Scrnc_jj_cc_lrspn.F's block-wise bare-v addition, which uses the same
C Isympert/Npert_count bookkeeping) back into the natural
C Iao_pair=(Mu-1)*Nao+Nu convention that Scrnc_wc_ao2mo.F and
C Scrnc_sigma_prep.F's Jj_bare (built independently via
C Scrnc_get2ints+Scrnc_sigma_bare_ao, genuinely natural-order) both
C assume. See SCRNC_DEBUG_LOG.md, 2026-08-18 "continued 18"/"continued
C 19": under sym=off (Nirrep=1) these two conventions coincide exactly
C (no permutation needed, invisible bug); under sym=on they diverge for
C most AO-pair slots, corrupting everything built from Jj_corr =
C Jj_lrespn - Jj_bare (Wc_mo, Wc_reduced, Sigma_c, eps_QP). This routine
C is the single fix point: called once, right after Jj_lrespn is built
C (either CC path), before Scrnc_wc_ao2mo_test/Scrnc_sigma_prep touch
C it, so neither of those (nor the already-verified sym=off path) needs
C any change.
C
C The permutation is rebuilt here from Iaopop/Dirprd alone (cheap,
C Nao_pair-sized) rather than reused from Scrnc_reord_cols.F's own
C Iordr_ao, which is local/transient to scrnc_prep_2ints.F and not
C saved anywhere reusable -- but mirrors that routine's exact
C construction, and was verified byte-for-byte against an independently
C captured bare-Coulomb matrix (2401/2401 AO-pair entries matched
C exactly for the H2O/STO-3G/C2v case, see debug log) before being
C written into this routine.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Maxcor)
      Double Precision JJ_lrespn(Nao_pair,Nao_pair)
C
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
C
      Dimension Irrep_aos(Nao), Irrep_ao_pairs(Nao_pair)
      Dimension Iordr_ao(Nao_pair)
C
C Nothing to do with only one irrep -- packed order already equals
C natural order (this is exactly why the bug was invisible under
C sym=off). Skip the scratch-buffer permute entirely in that case.
C
      If (Nirrep .Eq. 1) Return
C
      I_scratch = 1
      Iend      = I_scratch + Nao_pair*Nao_pair
      If (Iend .Gt. Maxcor) Call Insmem("Scrnc_jj_unpack_natural",
     +                                   Iend,Maxcor)
C
C Irrep_aos(Mu), Mu=1..Nao: ACES2's own AO numbering under symmetry is
C already irrep-blocked (first Iaopop(1) AOs = irrep 1, next Iaopop(2)
C = irrep 2, etc.) -- same assumption scrnc_prep_2ints.F's own
C Irrep_aos construction makes.
C
      Ithru = 0
      Do Irrep = 1, Nirrep
         Do Iaos = 1, Iaopop(Irrep)
            Ithru = Ithru + 1
            Irrep_aos(Ithru) = Irrep
         Enddo
      Enddo
C
C Irrep_ao_pairs(I), natural I=(Mu-1)*Nao+Nu (Mu outer, Nu inner) --
C matches scrnc_prep_2ints.F:68-77 exactly, including its
C Dirprd(Irrep_l,Irrep_k)=Dirprd(Irrep_of_Nu,Irrep_of_Mu) argument
C order (Dirprd is symmetric for this Abelian XOR-style group product,
C but kept identical for clarity/traceability).
C
      Ithru = 0
      Do Mu = 1, Nao
         Irrep_mu = Irrep_aos(Mu)
         Do Nu = 1, Nao
            Irrep_nu = Irrep_aos(Nu)
            Ithru = Ithru + 1
            Irrep_ao_pairs(Ithru) = Dirprd(Irrep_nu,Irrep_mu)
         Enddo
      Enddo
C
C Iordr_ao(Ipacked) = natural position -- matches Scrnc_reord_cols.F's
C own Iordr_ao construction exactly (scan natural order once per
C irrep, in increasing irrep order, appending matches in natural
C sub-order).
C
      Ithru = 0
      Do Irrepr = 1, Nirrep
         Do I = 1, Nao_pair
            If (Irrep_ao_pairs(I) .eq. Irrepr) Then
               Ithru = Ithru + 1
               Iordr_ao(Ithru) = I
            Endif
         Enddo
      Enddo
C
      If (Ithru .Ne. Nao_pair) Then
         Write(6,"(a,2I6)")
     +    " @-Scrnc_jj_unpack_natural: permutation count mismatch ",
     +    Ithru, Nao_pair
         Call Errex
      Endif
C
C Apply: JJ_natural(Iordr_ao(Ip),Iordr_ao(Iq)) = JJ_packed(Ip,Iq).
C Two-sided permute needs a scratch copy (can't safely do in place).
C
      Call Dcopy(Nao_pair*Nao_pair,JJ_lrespn,1,Work(I_scratch),1)
C
      Do Iq = 1, Nao_pair
       Iq_nat = Iordr_ao(Iq)
       Do Ip = 1, Nao_pair
        Ip_nat = Iordr_ao(Ip)
        JJ_lrespn(Ip_nat,Iq_nat) =
     +     Work(I_scratch-1+(Iq-1)*Nao_pair+Ip)
       Enddo
      Enddo
C
      Return
      End
