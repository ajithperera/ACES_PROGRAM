      Subroutine Scrnc_jj_cc_lrspn(Work,Maxcor,Iuhf,JJ_lrespn,Isympert,
     +                             Ioffpert,Nao_pair)

      Implicit Double Precision (A-H,O-Z)

#include "maxbasfn.par"
      
      Double Precision M_expect,Mbar_00
      Double Precision JJ_lrespn(Nao_pair,Nao_pair)
      Logical Incore
      Logical ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Logical Realfreq
      Integer Occ_pair,Occ_pair_4irrep
      Dimension Isympert(Nao_pair),Ioffpert(Nao_pair)
      Dimension Work(Maxcor)
#ifdef _DEBUG_LVL0
      Dimension D(Maxbasfn,Maxbasfn)
#endif 
      Character*8 Label(3),Work_label

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      Common/Sym_pqxx_ints/Irrep_aos(Maxbasfn),
     &                     Irrep_ao_pairs(Maxbasfn*Maxbasfn)
      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +       Drccd
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs 
      Common /eomprop/Ieomprop

      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

      Irreprl     = 1
      Npert_count = 0
      Do Irrep = 1, Nirrep
         Nperts_4irrepr = Irpdpdao(Irrep)

         Iao_pair_4irrep = 0
         Do Iao_pair= 1, Nperts_4irrepr

         Iao_pair_4irrep = Iao_pair_4irrep + 1 
         Npert_count     = Npert_count     + 1
         Irrep_ao_pair = Irrep_ao_pairs(Npert_count)
         Irrepr = Irrep_ao_pair
         Irrepl = Dirprd(Irrepr,Irreprl)
         Irrepx = Irrepl

#ifdef _DEBUG
      Nbfns   = Nocco(1) + Nvrto(1)
      Call scrnc_get_pert_type(Work,Maxcor,Nbfns,Label,Isympert,
     +                        Iuhf)
      Do Occ_pair= 1, 3
      Irrepx     = Isympert(occ_pair)
      Ipert      = occ_pair
      Work_label = Label(occ_pair)
#endif
         Nsizec = Irpdpd(Irrepx,9)
         If (Iuhf .EQ. 0) Then
            Nsizec = Nsizec  + Idsymsz(Irrepx,13,14)
         Else
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
            Nsizec = Nsizec + Idsymsz(Irrepx,1,3)
            Nsizec = Nsizec + Idsymsz(Irrepx,2,4)
            Nsizec = Nsizec + Idsymsz(Irrepx,13,14)
         Endif

         Call Scrnc_reset_lrspn_cc_lists(Work,Maxcor,Iuhf,Irrepx)

C List_ljki, List_abki and List_bjki are the OO, VV and OV integral
C list
         Lenoo = Irpdpd(Irrepx,21) + Iuhf*Irpdpd(Irrepx,22)
         Lenvv = Irpdpd(Irrepx,19) + Iuhf*Irpdpd(Irrepx,20)
         Lenvo = Irpdpd(Irrepx,9)  + Iuhf*Irpdpd(Irrepx,10)

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


         Call Scrnc_prep_culomb_ints(Work(Ioff_oo),Work(Ioff_vv),
     +                                Work(Ioff_vo),Iuhf,Lenoo,
     +                                Lenvv,Lenvo,Iao_pair_4irrep,
     +                                Irrepx,Irrepr)
#ifdef _DEBUG
      Nbfns   = Nocco(1) + Nvrto(1)
      Call Getrec(20, "JOBARC","NBASTOT ", 1, Naobfns)
      Length   = Naobfns * (Naobfns+1)/2
      Length21 = Nbfns * Nbfns
      Length22 = Nbfns * Naobfns
      I080 = I070 + Length
      I090 = I080 + Length21
      I100 = I090 + Length22
      I110 = I100 + Length22
      Iend = I110
      Mxcor = Maxcor - Iend
      If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
     +                                  Maxcor)
      Call Scrnc_prep_dipole_ints(Work(I070),Work(I080),Work(I090),
     +                           Work(I100),Work(I000),Work(I010),
     +                           Work(I020),Work(Iend),
     +                           Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
     +                           Naobfns,Iuhf,Ipert,Work_label,
     +                           Irrepx)
      Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
      Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
      Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
      Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)
#endif 
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

         If (Iuhf .eq. 0) Then 
            Call Scrnc_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                  Work(I050),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00,Ieomprop)
                                          
         Else
            Call Scrnc_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                  Work(I050),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00,Ieomprop)
         Endif 

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

         If (Iuhf .eq. 0) Then   
            Call Scrnc_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                  Work(I060),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00-Mexpect,
     +                                  Ieomprop)
         Else
            Call Scrnc_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                  Work(I060),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00-Mexpect,
     +                                  Ieomprop)
         Endif   

          If (Iside .Eq. 2) Then
             Mu_s    = 482
             Mu_d    = 453
             Ioff_s  = 0
             Ioff_d  = 0
             Call Scrnc_load_vec(Irrepx,Work(I000),Maxcor,Mu_s,
     +                           Ioff_s,Mu_d,Ioff_d,Iuhf,.False.)

#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("Jbar_L<(S,D)|pq>:",Work(I000),Nsizec)
#endif

C Save the left amplitude straight to its real per-perturbation record.
C Bug fix 2026-08-14 (see SCRNC_DEBUG_LOG.md): this used to go through a
C "scratch" record 1 of list 374 (Putlst(...,1,...), then
C Getlst(record 1)/Putlst(record Npert_count) below) -- but record 1 of
C list 374 is ALSO Npert_count=1's own permanent slot, so every later
C Npert_count's scratch write would silently clobber Npert_count=1's
C finished value. Same collision as the ring-CC bug fixed in
C scrnc_jj_rcc_lrspn.F; fixed identically here as a defensive measure
C -- this CCSD path's own control case (H2/STO-3G) never exposed it, but
C hadn't been tested thoroughly enough to trust on less-symmetric
C systems.

             Call Putlst(Work(I000),Npert_count,1,1,Irrepx,374)
          Endif

         Enddo ! Enddo for Iao_pair 
      Enddo ! Enddo of irrepr

      Call Dzero(JJ_lrespn,Nao_pair*Nao_pair)

      I_pairs = Nao_pair
      J_pairs = Nao_pair

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Symmetry of perturbations"
      Write(6,"(5(1x,I5))") (Isympert(i),i=1,Nao_pair)
      Write(6,"(a)") " Offset for perturbations"
      Write(6,"(5(1x,I5))") (Ioffpert(i),i=1,Nao_pair)
      Write(6,*)
#endif 
C Dynamic (frequency-dependent) response assembly, mirroring
C vcceh/secprp.F's ITOP/FACT combination exactly (2026-08-15): for a
C nonzero real frequency, Scrnc_lrspn_cc_main solved the R-side
C amplitude at +Freq (saved to list 375) and -Freq (left in list 373),
C separately; sum the two contributions here (FACT=1 each). For the
C static case (Freq=0) or an imaginary frequency, this collapses to a
C single pass reading list 373 with FACT=2 -- reproduces the original
C static code path exactly at Freq=0 (Itop=1, single Ddot, times 2).

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

      Do Icount = 1, Itop

      Listin = 373
      If (Icount .EQ. 1 .AND. Itop .EQ. 2) Listin = 375

      Do Iao_pair = 1, I_pairs

         Irrepx = Isympert(Iao_pair)
         Call Scrnc_newsize(Work,Maxcor,Iuhf,Irrepx,Nsizec)

         I000 = 1
         I010 = I000 + Nsizec
         I020 = I010 + Nsizec
         Iend = I020

         Memleft = Maxcor - Iend
         If (Iend .Gt. Maxcor)  Call Insmem("Scrnc_jj_lrspn",Iend,
     +                                      Memleft)

         Call Getlst(Work(I000),Iao_pair,1,1,Irrepx,Listin)
#ifdef _DEBUG_LVLM
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
     +                     374)
#ifdef _DEBUG_LVLM
               call checksum("374",work(i010),nsizec)
#endif
               X = Fact*Ddot(Nsizec,Work(I010),1,Work(I000),1)

#ifdef _DEBUG_LVL0
      If ((Iao_pair.EQ.9 .AND. Jao_pair.EQ.41) .OR.
     +    (Iao_pair.EQ.41 .AND. Jao_pair.EQ.9)) Then
         Nvopart = Min(20,Nsizec)
         Ndblpart = Nsizec - Nvopart
         Xvo  = Ddot(Nvopart,Work(I010),1,Work(I000),1)
         Xdbl = 0.0D0
         If (Ndblpart .GT. 0) Xdbl = Ddot(Ndblpart,Work(I010+Nvopart),
     +                                     1,Work(I000+Nvopart),1)
         Write(6,"(a,i3,a,i3,a,i5,a,F20.14,a,F20.14,a,F20.14)")
     +      " CROSSCHECK Iao_pair=",Iao_pair," Jao_pair=",Jao_pair,
     +      " Nsizec=",Nsizec," X=",X," Xvo(1:20)=",Fact*Xvo,
     +      " Xdbl(21:end)=",Fact*Xdbl
         call checksum("CROSSCHECK 373(Iao_pair)",Work(I000),Nsizec)
         call checksum("CROSSCHECK 374(Jao_pair)",Work(I010),Nsizec)
      Endif
#endif

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
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " CI-like JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Write(6,*)
#endif 
      
      If (Ieomprop .eq. 2) Then
          Itop = 1
          Call Calcquad1(Work(I000),JJ_lrespn,Maxcor,Isympert,
     +                   Ioffpert,Nao_pair,Iuhf,1.0D0,
     +                   Nao_pair,1,Itop,.False.)
      Endif

C Force-symmetrize the final W=v+v*chi*v matrix, matching the standard
C ACES2 EOM-effective-Hamiltonian property recipe (vcceh/secprp.F's
C Call Symmet2(XPROP,MAXPRT), J. Gauss 1990). 2026-08-15.

      Call Symmet2(JJ_lrespn,Nao_pair)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " LR JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Nbfns   = Nocco(1) + Nvrto(1)
      i000=1
      i010=i000+Nbfns*Nbfns 
      call getrec(20,"JOBARC","HFDENSTY",Nbfns*Nbfns,work(i000))
      ind=0
      do i=1,nbfns 
      do j=1,nbfns
         ind=ind+1
         d(j,i) = Work(i000-1+ind)
      enddo
      enddo
      ind=0
      do l=1,Nbfns
      do k=1,Nbfns
      do j=1,Nbfns 
      do i=1,Nbfns
         ind= ind + 1
         Work(i010-1+ind)= D(k,l)* D(i,j)
      enddo
      enddo
      enddo
      enddo
      nbfns2=Nbfns*Nbfns 
      write(6,*) "The SCF density matrix product"
      call output(work(i010),1,nbfns2,1,nbfns2,nbfns2,nbfsn2,1)
      e2 = Ddot(Nbfns2*Nbfns2,JJ_lrespn,1,Work(i010),1)
      Write(6,"(a,F15.9)") "The pseudo energy: ", e2
#endif 
    
      Return
      End 
