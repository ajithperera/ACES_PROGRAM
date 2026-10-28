










      Subroutine Scrnc_jj_rcc_assemble(Work,Maxcor,Iuhf,JJ_lrespn,
     +                                  Isympert,Ioffpert,Nao_pair,
     +                                  Nmult)

      Implicit Double Precision (A-H,O-Z)

      Double Precision JJ_lrespn(Nao_pair,Nao_pair)
      Logical Realfreq
      Logical Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Dimension Isympert(Nao_pair),Ioffpert(Nao_pair)
      Dimension Work(Maxcor)

      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +                Drccd
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

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
C
C SIGN (2026-10-15, SCRNC_DEBUG_LOG.md): secprp's pairing gives the
C polarizability-convention response, -chi (positive at w=0), but the
C bare-v add below forms W = v + v*chi*v and needs the physical
C chi(iw) <= 0. Fact is therefore negated. Before this fix every W_c,
C Sigma_c and QP energy had the correlation part with flipped sign
C (occupied QPs moved down while IP-EOM-CCSD moves them up); found by
C the Lehmann cross-check (Scrnc_lehmann.F), which matched the old
C W_c to ~0.2% in magnitude with opposite sign.
C
      Fact = -2.0D0
      If (Freq .NE. 0.0D0 .AND. Realfreq) Then
         Itop = 2
         Fact = -1.0D0
      Endif

      Call Dzero(JJ_lrespn,Nao_pair*Nao_pair)

C Diagnostic confirming Isympert(1..Nao_pair)'s ordering convention --
C see SCRNC_DEBUG_LOG.md 2026-08-18 sym=on/off investigation for the
C printed evidence this produced (kept 1-gated, matching this
C file's existing convention, rather than left permanently unguarded).
      Write(6,"(a)") " @-DIAG Iaopop(1..Nirrep):"
      Write(6,"(8(1x,I4))") (Iaopop(Ir),Ir=1,Nirrep)
      Write(6,"(a,I6)") " @-DIAG Nao_pair=",Nao_pair
      Write(6,"(a)") " @-DIAG Isympert(1..Nao_pair):"
      Write(6,"(10(1x,I3))") (Isympert(Ir),Ir=1,Nao_pair)

      I_pairs = Nao_pair
      J_pairs = Nao_pair

C 2026-08-19 (Ajith's call, see SCRNC_DEBUG_LOG.md): for the genuine
C 2-channel RHF RCCD case (Nmult=2, see Scrnc_driver.F), combine the
C singlet (Imult=1) and triplet (Imult=2) channels' response
C contributions with the standard 1-singlet/3-triplet spin-multiplicity
C weight, exactly matching Check_rhf_rcc.F's
C E_total=0.75*E_t+0.25*E_s reference formula. Loffset routes each
C channel to its own Scrnc_lrspn_rcc_main.F/Scrnc_jj_rcc_lrspn.F list
C numbers (373-376 singlet, 377-380 triplet). For Nmult=1 (UHF, or
C RHF DRCCD, unresolved separately), Wt=1 and Loffset=0 reproduce the
C original single-pass behavior exactly.

      Do Imult = 1, Nmult

      Wt = 1.0D0
      If (Nmult .EQ. 2) Then
         If (Imult .EQ. 1) Wt = 0.25D0
         If (Imult .EQ. 2) Wt = 0.75D0
      Endif
      Loffset = 4*(Imult-1)

      Do Icount = 1, Itop

      Listin_r = 373 + Loffset
      Listin_l = 374 + Loffset
      If (Icount .EQ. 1 .AND. Itop .EQ. 2) Then
         Listin_r = 375 + Loffset
         Listin_l = 376 + Loffset
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
         Write(6,"(a)") " <0|(1+L)M_bar and first-order e^T|0>"
         call checksum("373",work(i000),nsizec)
C 2026-08-19: RCCD vs DRCCD need DIFFERENT treatment here, confirmed
C empirically (see SCRNC_DEBUG_LOG.md) -- removing Spntsing helped
C RCCD (3.98%->2.45% vs UHF) but hurt DRCCD (1.44%->2.59%). Root
C cause: rcc_dwmbej.F's RHF list 56 (RCCD, now the only channel used,
C Scrnc_driver.F's Nmult=1) is the genuine, complete singlet RPA
C matrix (2J-K, confirmed by matching its own construction against the
C standard Casida singlet/triplet formulas) -- already fully spin-
C adapted, no further Spntsing-style transform needed or wanted.
C drcc_dwmbej.F's RHF list 56 (DRCCD) is different: direct-RPA drops
C the exchange term entirely (2J only, no K), and its Hbar-eigenvalue
C check (SCRNC_DEBUG_LOG.md) showed UHF's own direct-RPA Hbar still
C carries a SECOND eigenvalue RHF's single list-56 channel doesn't
C reproduce -- DRCCD's RHF response is missing real physics (an
C unbuilt second channel), a bigger, separate, not-yet-fixed problem.
C Spntsing was accidentally compensating for some of that gap, however
C incorrectly -- restore it for DRCCD only until DRCCD's real fix
C (a new Hbar piece) is built; RCCD never needs it.
         If (Iuhf .eq. 0 .AND. Drccd) Then
            Call Spntsing(Isympert(Iao_pair),Work(I000),Work(Iend),
     +                    Memleft)
         Endif

         Do Jao_pair = 1, J_pairs

            If (Isympert(Iao_pair) .eq. Isympert(Jao_pair)) Then

               Call Getlst(Work(I010),Jao_pair,1,1,Irrepx,
     +                     Listin_l)
               call checksum("374",work(i010),nsizec)
               X = Wt*Fact*Ddot(Nsizec,Work(I010),1,Work(I000),1)

               JJ_lrespn(Iao_pair,Jao_pair) =
     +                  JJ_lrespn(Iao_pair,Jao_pair) + X
            Endif

         Enddo
      Enddo

      Enddo ! Enddo for Icount (dynamic +/-Freq loop)

      Enddo ! Enddo for Imult (RHF RCCD singlet/triplet channels)

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
      Open(Unit=90,File="jjlrespn_plusbare_packed.dat",
     +     Form="Formatted",Status="UNKNOWN")
      Write(90,*) Nao_pair
      Do Jao_pair = 1, Nao_pair
       Do Iao_pair = 1, Nao_pair
        Write(90,*) JJ_lrespn(Iao_pair,Jao_pair)
       Enddo
      Enddo
      Close(Unit=90)

C Force-symmetrize the final W=v+v*chi*v matrix, matching the standard
C ACES2 EOM-effective-Hamiltonian property recipe (vcceh/secprp.F's
C Call Symmet2(XPROP,MAXPRT), J. Gauss 1990). 2026-08-15.

      Call Symmet2(JJ_lrespn,Nao_pair)

      Write(6,*)
      Write(6,"(a)") " RPA JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Write(6,*)

      Return
      End
