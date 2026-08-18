      Subroutine Scrnc_jj_rcc_assemble(Work,Maxcor,Iuhf,JJ_lrespn,
     +                                  Isympert,Ioffpert,Nao_pair)

      Implicit Double Precision (A-H,O-Z)

      Double Precision JJ_lrespn(Nao_pair,Nao_pair)
      Logical Realfreq
      Dimension Isympert(Nao_pair),Ioffpert(Nao_pair)
      Dimension Work(Maxcor)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

C Final JJ_lrespn assembly, split out of Scrnc_jj_rcc_lrspn.F
C (2026-08-15) so it can run ONCE after Scrnc_driver.F's dynamic
C +/-Freq Icount loop has populated lists 373/374 (the -Freq pass, or
C the single static solve) and, for a 2-pass dynamic run, 375/376 (the
C +Freq pass) -- mirroring vcceh/secprp.F's own ITOP/FACT combination
C (see scrnc_jj_cc_lrspn.F for the same pattern on the CCSD path).
C Freq is restored to its original sign by the time this runs (the
C driver's Icount loop negates it once per pass, an even number of
C times total for Itop=2).

      Itop = 1
      Fact = 2.0D0
      If (Freq .NE. 0.0D0 .AND. Realfreq) Then
         Itop = 2
         Fact = 1.0D0
      Endif

      Call Dzero(JJ_lrespn,Nao_pair*Nao_pair)

C Diagnostic confirming Isympert(1..Nao_pair)'s ordering convention --
C see SCRNC_DEBUG_LOG.md 2026-08-18 sym=on/off investigation for the
C printed evidence this produced (kept _DEBUG_LVL0-gated, matching this
C file's existing convention, rather than left permanently unguarded).
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-DIAG Iaopop(1..Nirrep):"
      Write(6,"(8(1x,I4))") (Iaopop(Ir),Ir=1,Nirrep)
      Write(6,"(a,I6)") " @-DIAG Nao_pair=",Nao_pair
      Write(6,"(a)") " @-DIAG Isympert(1..Nao_pair):"
      Write(6,"(10(1x,I3))") (Isympert(Ir),Ir=1,Nao_pair)
#endif

      I_pairs = Nao_pair
      J_pairs = Nao_pair

      Do Icount = 1, Itop

      Listin_r = 373
      Listin_l = 374
      If (Icount .EQ. 1 .AND. Itop .EQ. 2) Then
         Listin_r = 375
         Listin_l = 376
      Endif

      Do Iao_pair = 1, I_pairs

         Irrepx = Isympert(Iao_pair)

         Nsizec = Irpdpd(Irrepx,9)
         If (Iuhf .NE. 0) Then
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
         Endif

         I000 = 1
         I010 = I000 + Nsizec
         I020 = I010 + Nsizec
         Iend = I020

         Memleft = Maxcor - Iend
         If (Iend .Gt. Maxcor)  Call Insmem("Scrnc_jj_rcc_assemble",
     +                                      Iend,Memleft)

         Call Getlst(Work(I000),Iao_pair,1,1,Irrepx,Listin_r)
#ifdef _DEBUG_LVL0
         Write(6,"(a)") " <0|(1+L)M_bar and first-order e^T|0>"
         call checksum("373",work(i000),nsizec)
#endif
         If (Iuhf .eq. 0) Then
            Call Spntsing(Isympert(Iao_pair),Work(I000),Work(Iend),
     +                    Memleft)
         Endif

         Do Jao_pair = 1, J_pairs

            If (Isympert(Iao_pair) .eq. Isympert(Jao_pair)) Then

               Call Getlst(Work(I010),Jao_pair,1,1,Irrepx,
     +                     Listin_l)
#ifdef _DEBUG_LVL0
               call checksum("374",work(i010),nsizec)
#endif
               X = Fact*Ddot(Nsizec,Work(I010),1,Work(I000),1)

               JJ_lrespn(Iao_pair,Jao_pair) =
     +                  JJ_lrespn(Iao_pair,Jao_pair) + X
            Endif

         Enddo
      Enddo

      Enddo ! Enddo for Icount (dynamic +/-Freq loop)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Response-only JJ matrix (before bare v added)"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
C Diagnostic (2026-08-18, see scrnc_sigma_prep.F's own comment) -- plain
C dump of the packed-order response-only JJ_lrespn, before the bare-v
C add, for empirical permutation recovery.
      Open(Unit=89,File="jjlrespn_resp_packed.dat",Form="Formatted",
     +     Status="UNKNOWN")
      Write(89,*) Nao_pair
      Do Jao_pair = 1, Nao_pair
       Do Iao_pair = 1, Nao_pair
        Write(89,*) JJ_lrespn(Iao_pair,Jao_pair)
       Enddo
      Enddo
      Close(Unit=89)
#endif

C --------------------------------------------------------------------
C Add the bare (uncorrelated) Coulomb integral (mu nu|kappa lambda) so
C that JJ_lrespn holds the FULL screened Coulomb interaction
C    W = v + v*Chi*v
C rather than just the screening correction v*Chi*v computed above.
C List_xxxx (record 220, built in scrnc_prep_2ints.F) is the symmetry
C-packed bare AO two-electron integral list, indexed identically
C (irrep, within-irrep column) to JJ_lrespn's own Iao_pair bookkeeping,
C so it adds in block by block, irrep by irrep. 2026-08-12.
C --------------------------------------------------------------------

      List_xxxx   = 220
      I000        = 1
      Npert_count = 0
      Do Irrep = 1, Nirrep
         Nperts_4irrepr = Irpdpdao(Irrep)
         Irow_base      = Npert_count
         Do Iao_pair_4irrep = 1, Nperts_4irrepr
            Npert_count = Npert_count + 1
            Call Getlst(Work(I000),Iao_pair_4irrep,1,1,Irrep,
     +                  List_xxxx)
            Do Irow = 1, Nperts_4irrepr
               JJ_lrespn(Irow_base+Irow,Npert_count) =
     +              JJ_lrespn(Irow_base+Irow,Npert_count) +
     +              Work(I000-1+Irow)
            Enddo
         Enddo
      Enddo

C Diagnostic (2026-08-18) -- plain dump right after the bare-v add,
C BEFORE Symmet2, so (this dump) - (jjlrespn_resp_packed.dat) recovers
C the bare part alone, in packed order, unaffected by symmetrization.
#ifdef _DEBUG_LVL0
      Open(Unit=90,File="jjlrespn_plusbare_packed.dat",
     +     Form="Formatted",Status="UNKNOWN")
      Write(90,*) Nao_pair
      Do Jao_pair = 1, Nao_pair
       Do Iao_pair = 1, Nao_pair
        Write(90,*) JJ_lrespn(Iao_pair,Jao_pair)
       Enddo
      Enddo
      Close(Unit=90)
#endif

C Force-symmetrize the final W=v+v*chi*v matrix, matching the standard
C ACES2 EOM-effective-Hamiltonian property recipe (vcceh/secprp.F's
C Call Symmet2(XPROP,MAXPRT), J. Gauss 1990). 2026-08-15.

      Call Symmet2(JJ_lrespn,Nao_pair)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " RPA JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Write(6,*)
#endif

      Return
      End
