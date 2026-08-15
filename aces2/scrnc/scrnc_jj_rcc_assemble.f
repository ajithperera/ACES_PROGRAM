










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
      Fact = 2.0D0
      If (Freq .NE. 0.0D0 .AND. Realfreq) Then
         Itop = 2
         Fact = 1.0D0
      Endif

      Call Dzero(JJ_lrespn,Nao_pair*Nao_pair)

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
         Write(6,"(a)") " <0|(1+L)M_bar and first-order e^T|0>"
         call checksum("373",work(i000),nsizec)
         If (Iuhf .eq. 0) Then
            Call Spntsing(Isympert(Iao_pair),Work(I000),Work(Iend),
     +                    Memleft)
         Endif

         Do Jao_pair = 1, J_pairs

            If (Isympert(Iao_pair) .eq. Isympert(Jao_pair)) Then

               Call Getlst(Work(I010),Jao_pair,1,1,Irrepx,
     +                     Listin_l)
               call checksum("374",work(i010),nsizec)
               X = Fact*Ddot(Nsizec,Work(I010),1,Work(I000),1)

               JJ_lrespn(Iao_pair,Jao_pair) =
     +                  JJ_lrespn(Iao_pair,Jao_pair) + X
            Endif

         Enddo
      Enddo

      Enddo ! Enddo for Icount (dynamic +/-Freq loop)

      Write(6,*)
      Write(6,"(a)") " Response-only JJ matrix (before bare v added)"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)

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
