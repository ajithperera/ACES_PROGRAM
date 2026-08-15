










      Subroutine Scrnc_jj_rcc_lrspn(Work,Maxcor,Iuhf,JJ_lrespn,Isympert,
     +                             Ioffpert,Nao_pair,Hbarinv,Iloc)

      Implicit Double Precision (A-H,O-Z)

      Double Precision M_expect,Mbar_00
      Double Precision JJ_lrespn(Nao_pair,Nao_pair)
      Logical Incore
      Logical ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Logical Realfreq
      Integer Occ_pair,Occ_pair_4irrep
      Dimension Isympert(Nao_pair),Ioffpert(Nao_pair)
      Dimension Work(Maxcor)
      Dimension Hbarinv(*)
      Dimension Iloc(8)
      Character*8 Label(3),Work_label
      Data One, Zilch /1.0D0,0.0D0/

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)

      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +       Drccd
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs 
      Common /eomprop/Ieomprop

      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/



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

      Irreprl     = 1
      Npert_count = 0

      Write(6,"(1x,a)") "-----Entering scrnc_jj_rcc_lrspn------"

C The left-hand (L) response also needs multiplying by Hbar^(-1), same
C as the right-hand (R) response does in scrnc_lrspn_rcc_main.F -- this
C was previously entirely missing (L was just the raw, un-inverted
C Mbar operator), which broke JJ_lrespn's symmetry (large error for
C strongly-coupled AO pairs, negligible for weakly-coupled ones, since
C skipping the inversion is closer to correct when Hbar is close to
C diagonal). 2026-08-12.
C
C 2026-08-14: Hbar^(-1) is now built ONCE by the caller
C (Scrnc_driver.F) and passed in via Hbarinv/Iloc, instead of being
C independently rebuilt here. Scrnc_form_rcc_hbar.F hands its result to
C Scrnc_form_rcc_hbar_inv.F by writing it back into the same list it
C read the raw integrals from, so a second, independent build (as this
C routine used to do) silently re-added the one-particle Fock diagonal
C on top of the first build's already-dressed output -- breaking L/R
C Hbar consistency and hence the symmetry of the final W=v+v*chi*v
C matrix.

      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irreprl)
         Irrepx = Irrepl

         Nsizec = Irpdpd(Irrepx,9)
         If (Iuhf .NE. 0) Then
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
         Endif

         Ifnd = Iloc(irrepr)

         Call Scrnc_reset_lrspn_rcc_lists(Work,Maxcor,Iuhf,Irrepx)

C List_ljki, List_abki and List_bjki are the OO, VV and OV integral
C list
         Lenoo = Irpdpd(Irrepx,21) + Iuhf*Irpdpd(Irrepx,22)
         Lenvv = Irpdpd(Irrepx,19) + Iuhf*Irpdpd(Irrepx,20)
         Lenvo = Irpdpd(Irrepx,9)  + Iuhf*Irpdpd(Irrepx,10)

         Nperts_4irrepr = Irpdpdao(Irrepr)
 
         Iao_pair_4irrep = 0
         Do Iao_pair= 1, Nperts_4irrepr

         I000 = 1
         I010 = I000 + Lenoo
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Lenoo
         I050 = I040 + Lenvv
         I060 = I050 + Lenvo
         I070 = I060 + Lenvo
         Iend = I070
         Mxcor = Maxcor - Iend 

         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_jj_lrspn",Iend,
     +                                      Maxcor)
         Ioff_oo = I000
         Ioff_vv = I010 
         Ioff_vo = I020 

         Iao_pair_4irrep = Iao_pair_4irrep + 1 
         Npert_count     = Npert_count     + 1

         Call Scrnc_prep_culomb_ints(Work(Ioff_oo),Work(Ioff_vv),
     +                               Work(Ioff_vo),Iuhf,Lenoo,
     +                               Lenvv,Lenvo,Iao_pair_4irrep,
     +                               Irrepx,Irrepr)
         Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
         Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)

C We construct e^-TMue(pq)^T, similarity transformed Mu(p,q).
            
         Iside = 1
         Call Scrnc_form_mubar_00(Work(Ioff_oo),Work(Ioff_vv),
     +                            Work(Ioff_vo),Work(Iend),
     +                            Mxcor,Iuhf,Irrepx,M_expect,
     +                            Mbar_00)

         Call Scrnc_form_mubar_s(Work(Ioff_oo),Work(Ioff_vv),
     +                           Work(Ioff_vo),Work(I030),
     +                           Work(I040),Work(I050),
     +                           Work(I060),Work(Iend),Mxcor,
     +                           Iuhf,Irrepx,Iside,Mbar_00,
     +                           Ieomprop)

C R-side/L-routine cross-checkpoint (SCRNC_DEBUG_LOG.md step 1): this is
C this routine's own Iside=1 stage, on the same raw integrals for this
C Npert_count as Scrnc_lrspn_rcc_main.F's (proven-correct) R-side. Diff
C against that routine's matching checksum before trusting the Iside=2
C stage below.
      write(6,"(1x,a,i4)") "L-routine Npert_count=",Npert_count
      call checksum("L-routine Mbar_00:",Mbar_00,1)
      call checksum("L-routine Mubar_s(Iside=1):",Work(I030),
     +              Lenoo+Lenvv+2*Lenvo)
      Write(6,"(a,20F15.10)") " RAWTRACE Work(I050) right after "//
     +     "Iside=1 returns:",(Work(I050-1+I),I=1,MIN(Lenvo,20))

         Iside = 2
         Ibgn  = I070 
         I080  = I070 + Lenoo
         I090  = I080 + Lenvv
         I100  = I090 + Lenvo
         I110  = I100 + Lenvo
         Iend  = I110
         Mxcor = Maxcor - Iend 
         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
     +                                     Maxcor)

         Call Dcopy(Lenoo,Work(I000),1,Work(I070),1)
         Call Dcopy(Lenvv,Work(I010),1,Work(I080),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I090),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I100),1)
 
         Call Scrnc_form_mubar_s(Work(I030),Work(I040),
     +                           Work(I050),Work(I070),
     +                           Work(I080),Work(I090),
     +                           Work(I100),Work(Iend),Mxcor,   
     +                           Iuhf,Irrepx,Iside,Mbar_00-
     +                           M_expect,Ieomprop)

         If (Iside .Eq. 2) Then
             Mu_s    = 482
             Lenr1a = Irpdpd(Irrepx,9)
             Lenr1b = Irpdpd(Irrepx,10)

             I000  = 1
             Iloc2 = I000 + Lenr1a
             Iwrk  = Iloc2 + Lenr1b
             Call Getlst(Work(I000),1,1,1,1,Mu_s)
             If (Iuhf .NE. 0) Then
                 Call Getlst(Work(Iloc2),1,1,1,2,Mu_s)
             Endif

C Hbar is NOT symmetric (EOM-CC has distinct left/right eigenvectors),
C so the L-side needs (Hbar^-1)^T, not the same Hbar^-1 used for R.
C (A^-1)^T = (A^T)^-1 always, so transposing the already-computed
C Hbar^-1 here is exact, not an approximation. 2026-08-12.
             Call Xgemm("T","N",Nsizec,1,Nsizec,One,Hbarinv(Ifnd),
     +                  Nsizec,Work(I000),Nsizec,Zilch,Work(Iwrk),
     +                  Nsizec)

      write(6,*)
      call checksum("Jbar_L<(S,D)|pq>:",Work(Iwrk),Nsizec)

C Save the left amplitude straight to its real per-perturbation record.
C Bug fix 2026-08-14 (see SCRNC_DEBUG_LOG.md): this used to go through a
C "scratch" record 1 of list 374 (Putlst(...,1,...) here, then
C Getlst(record 1)+Putlst(record Npert_count) below) -- but record 1 of
C list 374 IS ALSO Npert_count=1's own permanent slot, so every later
C Npert_count's scratch write silently clobbered Npert_count=1's
C finished value with whatever the LAST perturbation's result happened
C to be. That's what broke JJ_lrespn's symmetry (record 1 ended up
C holding Npert_count=Nao_pair's value instead of its own). Writing
C directly to Npert_count (already known here) removes the collision.

          Call Putlst(Work(Iwrk),Npert_count,1,1,Irrepx,374)
          Endif

         Enddo ! Enddo for Iao_pair
      Enddo ! Enddo of irrepr

      Call Dzero(JJ_lrespn,Nao_pair*Nao_pair)

      I_pairs = Nao_pair
      J_pairs = Nao_pair

      Write(6,*)
      Write(6,"(a)") " Symmetry of perturbations"
      Write(6,"(5(1x,I5))") (Isympert(i),i=1,Nao_pair)
      Write(6,"(a)") " Offset for perturbations"
      Write(6,"(5(1x,I5))") (Ioffpert(i),i=1,Nao_pair)
      Write(6,*)
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
         If (Iend .Gt. Maxcor)  Call Insmem("Scrnc_jj_lrspn",Iend,
     +                                      Memleft)

         Call Getlst(Work(I000),Iao_pair,1,1,Irrepx,373)
         Write(6,"(a)") " <0|(1+L)M_bar and first-order e^T|0>"
         call checksum("373",work(i000),nsizec)
         Write(6,"(a,I3,a,10F15.10)") " 373 raw, Iao_pair=",Iao_pair,
     +                             ":",(Work(I000-1+I),I=1,Nsizec)
         If (Iuhf .eq. 0) Then
            Call Spntsing(Isympert(Iao_pair),Work(I000),Work(Iend),
     +                    Memleft) 
         Endif 

         Do Jao_pair = 1, J_pairs 

            If (Isympert(Iao_pair) .eq. Isympert(Jao_pair)) Then
 
               Call Getlst(Work(I010),Jao_pair,1,1,Irrepx,
     +                     374)
               call checksum("374",work(i010),nsizec)
               Write(6,"(a,I3,a,10F15.10)") " 374 raw, Jao_pair=",
     +                             Jao_pair,":",
     +                             (Work(I010-1+I),I=1,Nsizec)
               X = Ddot(Nsizec,Work(I010),1,Work(I000),1)

               JJ_lrespn(Iao_pair,Jao_pair) = 
     +                  JJ_lrespn(Iao_pair,Jao_pair) + X
            Endif 
              
         Enddo
      Enddo

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
C Call Symmet2(XPROP,MAXPRT), J. Gauss 1990): the L/R response
C construction above (T2 vs ringcl's ground-state Lambda2, which are
C never identical -- distinct left/right EOM/effective-Hamiltonian
C quantities) leaves JJ_lrespn with a small numerical asymmetry by
C construction, same as vcceh's XPROP before its own Symmet2 call.
C vcceh's established practice is to symmetrize the finished property
C matrix explicitly rather than require exact symmetry to fall out of
C the L/R construction itself. 2026-08-15.

      Call Symmet2(JJ_lrespn,Nao_pair)

      Write(6,*)
      Write(6,"(a)") " RPA JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Write(6,*)

      Return
      End
