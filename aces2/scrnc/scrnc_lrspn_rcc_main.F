      Subroutine Scrnc_lrspn_rcc_main(Work,Maxcor,Iuhf,Isympert,
     +                                Ioffpert,Nperts,Hbarinv,Iloc,
     +                                Icount,Itop)

      Implicit Double Precision (A-H,O-Z)

#include "maxbasfn.par"

      Double Precision M_expect,Mbar_00
      Logical Incore
      Logical ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Logical Realfreq
      Logical Rcc
      Integer Occ_pair,Occ_pair_4irrep
      Dimension Isympert(Nperts),Ioffpert(Nperts)
      Dimension Work(Maxcor)
      Dimension Hbarinv(*)
      Dimension Iloc(8)
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
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs 
      Common /eomprop/Ieomprop

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data One, Zilch /1.0D0,0.0D0/
     
#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering scrnc_lrspn_rcc_main------"
#endif
C In the case of RPA like situations, we do not necessarily need
C to solve linear equations (smaller problem). We have Hbar C = Y.
C We can simply write C = Hbar^(-1) Y. Hbar^(-1) is perturbation
C independent, so it is built ONCE by the caller (Scrnc_driver.F) and
C shared with Scrnc_jj_rcc_lrspn.F -- see Hbarinv/Iloc arguments.
C 2026-08-14: this used to be built independently here AND in
C Scrnc_jj_rcc_lrspn.F; since Scrnc_form_rcc_hbar.F hands its result to
C Scrnc_form_rcc_hbar_inv.F by writing it back into the same list it
C read the raw integrals from, the second independent build silently
C re-added the one-particle Fock diagonal on top of the first build's
C already-dressed output, breaking L/R Hbar consistency.

      Irreprl     = 1
      Npert_count = 0

C Next step is to prepare the integrals. We have the repulsion
C integrals prepared in the form X(p(1)q(1),i(2)j(2)) where p and q can
C be occ-occ,vrt-vrt or occ-vrt and i(2) and j(2) always kept as AO
C We process the p(1)q(1) integrals for a given pair of i(2)j(2)
C and lets designate them as Mu(pq).

C 2026-08-15: list creation (Scrnc_form_lrspn_rcc_lists) moved to
C Scrnc_driver.F, called ONCE before the dynamic +/-Freq Icount loop,
C instead of here -- this routine is now called once PER Icount pass,
C and re-touching/re-opening the underlying list files on the second
C pass was found to corrupt the (otherwise Freq-independent) raw Hbar
C data that Scrnc_form_rcc_hbar_inv.F reads (same class of bug as the
C original list-374 collision, one level removed: a shared underlying
C file getting reset by a second Updmoi/list-creation call). See
C SCRNC_DEBUG_LOG.md, 2026-08-15 dynamic-ring-CC entry.

      Do Irrepr = 1, Nirrep
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
         If (Iuhf .NE. 0) Then
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
         Endif

         Call Scrnc_reset_lrspn_rcc_lists(Work,Maxcor,Iuhf,Irrepx)

C The R(DR)CCD HBAR^(-1) is saved in first Nrow * Ncol position of Work
C array. Form the 1^st order perturbed amplitudes by HBAR(^-1)

         If (Iuhf .EQ. 0) Then
            Nrow  = Irpdpd(Irrepx,9)
            Ncol  = Nrow
            Nmult = 2
         Else
            Nrow  = Irpdpd(Irrepx,9) + Irpdpd(Irrepx,10)
            Ncol  = Nrow
            Nmult = 1
         Endif 

         Ifnd = Iloc(irrepr)
         Ibgn = 1
         I000 = Ibgn + Nrow * Ncol

         Call Dcopy(Nrow*Ncol,Hbarinv(Ifnd),1,Work(Ibgn),1)
#ifdef _DEBUG_LVL0
         write(6,*) " @-Scrnc_lrspn_rcc_main Hbar^(-1)"
         call output(work(Ibgn),1,nrow,1,ncol,nrow,ncol,1)
#endif 

C List_ljki, List_abki and List_bjki are the OO, VV and OV integral
C list
         Lenoo = Irpdpd(Irrepx,21) + Iuhf*Irpdpd(Irrepx,22)
         Lenvv = Irpdpd(Irrepx,19) + Iuhf*Irpdpd(Irrepx,20)
         Lenvo = Irpdpd(Irrepx,9)  + Iuhf*Irpdpd(Irrepx,10)

         Nperts_4irrepr = Irpdpdao(Irrepr)

         Iao_pair_4irrep = 0
         Do Iao_pair= 1, Nperts_4irrepr

         Npert_count           = Npert_count     + 1
         Iao_pair_4irrep       = Iao_pair_4irrep + 1
         Isympert(Npert_count) = Irrepx
         Ioffpert(Npert_count) = Npert_count 
          
         I010 = I000 + Lenoo
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Lenoo
         I050 = I040 + Lenvv
         I060 = I050 + Lenvo
         I070 = I060 + Lenvo
         Iend = I070

         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
     +                                      Maxcor)
         Ioff_oo = I000
         Ioff_vv = I010
         Ioff_vo = I020

          Call Scrnc_prep_culomb_ints(Work(Ioff_oo),Work(Ioff_vv),
     +                                Work(Ioff_vo),Iuhf,Lenoo,
     +                                Lenvv,Lenvo,Iao_pair_4irrep,
     +                                Irrepx,Irrepr)

          Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
          Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
          Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
          Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)

C We construct e^-TMue(pq)^T, similarity transformed Mu(p,q).

          Do ISide = 1, 1
               
             Call Scrnc_form_mubar_00(Work(Ioff_oo),Work(Ioff_vv),
     +                                Work(Ioff_vo),Work(Iend),
     +                                Maxcor,Iuhf,Irrepx,M_expect,
     +                                Mbar_00)

             Call Scrnc_form_mubar_s(Work(Ioff_oo),Work(Ioff_vv),
     +                               Work(Ioff_vo),Work(I030),
     +                               Work(I040),Work(I050),
     +                               Work(I060),Work(Iend),Maxcor,
     +                               Iuhf,Irrepx,Iside,Mbar_00,
     +                               Ieomprop)

C R-side/L-routine cross-checkpoint (SCRNC_DEBUG_LOG.md step 1): the
C same raw integrals for this Npert_count feed an Iside=1 mubar_s stage
C in Scrnc_jj_rcc_lrspn.F too, before that routine's own Iside=2 stage
C diverges. R-side is proven correct (RPA-sourced-T2 elimination test),
C so these two checksums should match theirs exactly for the same
C Npert_count -- diff to confirm before narrowing into Iside=2.
#ifdef _DEBUG_LVL0
      write(6,"(1x,a,i4)") "R-side Npert_count=",Npert_count
      call checksum("R-side Mbar_00:",Mbar_00,1)
      call checksum("R-side Mubar_s(Iside=1):",Work(I030),
     +              Lenoo+Lenvv+2*Lenvo)
#endif
             If (Iside .eq. 1) Then

                 Mu_s   = 480
                 Lenr1a = Irpdpd(Irrepx,9)
                 Lenr1b = Irpdpd(Irrepx,10)

                 I010 = I000
                 I020 = I010 + Lenr1a 
                 Call Getlst(Work(I010),1,1,1,1,Mu_s)
                 If (Iuhf .NE. 0) Then
                    I030 = I020 + Lenr1b
                    Call Getlst(Work(I020),1,1,1,2,Mu_s)
                 Endif 

                 I040 = I030 + Lenr1a + (Iuhf-1) * Lenr1b
                 Iend = I040
                 If (Iend .Ge. Maxcor) 
     +              Call Insmem("Scrnc_lrspn_rcc_main",Iend,Maxcor)

                 Call Xgemm("N","N",Nrow,1,Nrow,One,Work(Ibgn),
     +                       Nrow,Work(I010),Nrow,Zilch,Work(I030),
     +                       Nrow)
                 
C Save the 1^st order perturbed amplitude for all perturbations

                 Call Putlst(Work(I030),Npert_count,1,1,Irrepx,
     +                       372+Iside)

C Dynamic (frequency-dependent) response, 2026-08-15: on the +Freq pass
C of a 2-pass dynamic run (see Scrnc_driver.F's Icount loop), also save
C to list 375 (secprp's own convention, mirrors the CCSD path) before
C the -Freq pass overwrites list 372+Iside=373 in place.

                 If (Icount .EQ. 1 .AND. Itop .EQ. 2) Then
                    Call Putlst(Work(I030),Npert_count,1,1,Irrepx,375)
                 Endif

#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("1^st order amps,Iside=1:",Work(I030),Nsizec)
#endif
             Else

                 Mu_s   = 482
                 Lenr1a = Irpdpd(Irrepx,9)
                 Lenr1b = Irpdpd(Irrepx,10)

                 I010 = I000
                 I020 = I000 + Lenr1a 
                 Call Getlst(Work(I010),1,1,1,1,Mu_s)
                 If (Iuhf .NE. 0) Then
                    I030 = I020 + Lenr1b
                    Call Getlst(Work(I020),1,1,1,2,Mu_s)
                 Endif 

                 I040 = I030 + Lenr1a + (Iuhf-1) * Lenr1b
                 Iend = I040
                 If (Iend .Ge. Maxcor) 
     +              Call Insmem("Scrnc_lrspn_rcc_main",Iend,Maxcor)

                 Call Xgemm("N","N",Nrow,1,Nrow,One,Work(Ibgn),
     +                       Nrow,Work(I010),Nrow,Zilch,Work(I030),
     +                       Nrow)

C Save the 1^st order perturbed amplitude for all perturbations

                 Call Putlst(Work(I030),Npert_count,1,1,Irrepx,
     +                       372+Iside)

#ifdef _DEBUG_LVL0
      write(6,*)
      call checksum("1^st order amps,Iside=2",Work(I030),Nsizec)
#endif
             Endif 

          Enddo ! Enddo for Iside
        Enddo ! Enddo for Iao_pair 

      Enddo ! Enddo of irrepr

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") " Symmetry of perturbations"
      Write(6,"(5(1x,I5))") (Isympert(i),i=1,Nperts)
      Write(6,"(a)") " Offset for perturbations"
      Write(6,"(5(1x,I5))") (Ioffpert(i),i=1,Nperts)
      stop
#endif
    
      Return
      End 
