













































































































































































































      SUBROUTINE Scrnc_driver(Work,Maxmem,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Logical Ccsd,Mbpt,Parteom,Nodavid,Rccd,Drccd
      Logical Realfreq
      Dimension Work(Maxmem)
      Dimension Iloc(8)
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)


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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
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
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +                Drccd

      Irrepx = 1
      Call scrnc_set_vars(Iuhf)
      Call Scrnc_getaoinf(Iuhf,Irrepx)

      Call Scrnc_bcktran_test(Work,Maxmem,Iuhf)

      Call Scrnc_make_lists(Iuhf,Irrepx,List_aibj,List_xibj,
     +                      List_xxbj,List_xxxx,List_ijxx,
     +                      List_abxx,List_aixx,List_iaxx)

C Read and process atomic orbital integrals and write them to
C list_xxx.

      Call Scrnc_prep_2ints(Work,Maxmem,Irrepx,List_xxxx) 

C (XX,XX)->(IX,XX)->(IJ,XX) (ISPIN=1)
C (XX,XX)->(iX,XX)->(ij,XX) (ISPIN=2)

      Listao_src = List_xxxx
      Listmo_tar = List_ijxx
      Icase      = 3

      Call Scrnc_xxxx_pqxx(Work,Maxmem,Iuhf,Listao_src,Listmo_tar,
     +                     Irrepx,Icase)

C (XX,XX)->(AX,XX)->(AB,XX) (ISPIN=1)
C (XX,XX)->(BX,XX)->(ab,XX) (ISPIN=2)

      Listmo_tar = List_abxx
      Icase      = 4

      Call Scrnc_xxxx_pqxx(Work,Maxmem,Iuhf,Listao_src,Listmo_tar,
     +                     Irrepx,Icase)

C (XX,XX)->(AX,XX)->(AI,XX) (ISPIN=1)
C (XX,XX)->(ax,XX)->(AI,XX) (ISPIN=2)

      Listmo_tar = List_aixx
      Icase      = 5

      Call Scrnc_xxxx_pqxx(Work,Maxmem,Iuhf,Listao_src,Listmo_tar,
     +                     Irrepx,Icase)

C (AI,XX) = (IA,XX) (ISPIN=1)
C (ai,XX) = (ai,XX) (ISPIN=2)

      Do Ispin = 1, Iuhf+1
         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Nrow   = Irpdpd(Irrepl,8+Ispin)
            Ncol   = Irpdpdao(Irrepr)
            Nsize  = Nsize + Nrow * Ncol
            If (Nsize .Gt. Maxmem) Call Insmem("scrnc_driver",
     +                                          Nsize,Maxmem)
            Call Getlst(Work,1,Ncol,1,Irrepr,List_aixx+Ispin)
            Call Putlst(Work,1,Ncol,1,Irrepr,List_iaxx+Ispin)
         Enddo 
      Enddo 

C If Rccd or DrCCD is used form T1 and L1 lists and set them to zero,

      If (Rccd .OR. Drccd) Then
         Call Scrnc_rccd_dummy_inits(Work,Maxmem,Iuhf)
         Call Scrnc_form_rccd_g2(Work,Maxmem,Iuhf)
      Endif 

      If (Iuhf.eq.0) Call Spinad56(Work,Maxmem,Iuhf)

      Call Modf(Work,Maxmem,Iuhf,1)

      Write(6,"(a)") " Hbar elements"
      Call checkhbar(Work,Memleft,Iuhf)
      Call Scrnc_built_respn_den(Work,Maxmem,Iuhf)

C The following codes compute the linear response (lsrnc) of two Coulomb 
c operators. 

      Nao = 0
      Do Irrepr = 1, Nirrep
         Nao = Nao + Iaopop(Irrepr)
      Enddo

      Nao_pairs = Nao * Nao
   
      Write(6,*)
      Write(6,"(a,I6)") " The number of ao pairs:", Nao_pairs
      
      I000 = 1
      Iend = Maxmem

C Build Hbar^(-1) ONCE here and share it between the R-side and L-side
C response routines below, instead of letting each of them build it
C independently. Scrnc_form_rcc_hbar.F hands its result to
C Scrnc_form_rcc_hbar_inv.F by writing it back into the SAME list it
C read the raw (ai,bj) integrals from -- so a second, independent
C rebuild would silently re-add the one-particle F(ab)d(ij)-F(ij)d(ab)
C diagonal on top of the first build's already-dressed output, breaking
C L/R Hbar consistency and hence the symmetry of the final
C W=v+v*chi*v matrix (confirmed via instrumentation: the second build's
C diagonal came out exactly one SCF orbital-energy gap too large).
C Building it once and sharing it is both cheaper and correct by
C construction. 2026-08-14.

      Ihbarinv_bgn = 1
      If (Rccd .OR. Drccd) Then
         Ndim  = 0
         Iloc(1) = 1
         Do Irrepr = 2, Nirrep
            Length = Irpdpd(Irrepr-1,9)
            If (Iuhf .NE. 0) Length = Length + Irpdpd(Irrepr-1,10)
            Length = Length * Length
            Iloc(Irrepr) = Iloc(Irrepr-1) + Length
         Enddo
         Do Irrepr = 1, Nirrep
            Ndim = Ndim + Irpdpd(Irrepr,9)
            If (Iuhf .NE. 0) Ndim = Ndim + Irpdpd(Irrepr,10)
         Enddo
C Scrnc_form_rcc_hbar/_inv need scratch space beyond the persistent
C Ndim*Ndim+Ndim Hbar^(-1) storage itself (diagonal buffers, SSTGEN
C scratch, MINV pivot arrays) -- give them the same generous headroom
C the old per-caller Maxcor budget provided, not just the bare
C Ndim*Ndim+Ndim. Only the first Ineed words (read back via
C Hbarinv(Iloc(irrepr)) in the callees below) need to survive after
C this call returns; the padding above that is free scratch that gets
C overwritten before it would ever be read as data.
         Ineed        = Ndim * Ndim + Ndim
         Ipad         = 10000 + 10 * Ineed
         Ihbarinv_bgn = Iend - (Ineed + Ipad) + 1
         Iend         = Ihbarinv_bgn - 1

C Separate, pristine copy of the raw (dressed) Hbar (2026-08-15), built
C ONCE by Scrnc_form_rcc_hbar+Scrnc_rcc_hbar_read before the dynamic
C Icount loop and never rewritten again -- Scrnc_form_rcc_hbar_inv now
C only COPIES from this into Work(Ihbarinv_bgn) each pass. The actual
C bug (found via isolation testing, see SCRNC_DEBUG_LOG.md): the
C original Scrnc_rcc_hbar_inv is a WRAPPER (Scrnc_rcc_hbar_inv.F, easy
C to miss -- named almost identically to Scrnc_form_rcc_hbar_inv.F)
C that calls Scrnc_form_rcc_hbar (destructively adds the one-particle
C Fock diagonal IN PLACE via GETALL-then-PUTALL to the SAME list) THEN
C Scrnc_form_rcc_hbar_inv. Calling that wrapper more than once (which
C the dynamic Icount loop needs to do, once per +/-Freq pass) added the
C Fock diagonal again on every call. Fix: call Scrnc_form_rcc_hbar
C exactly once (below, sharing this same scratch buffer, same as the
C old wrapper did), then Scrnc_rcc_hbar_read once to snapshot the
C now-correctly-dressed Hbar -- same Ineed+Ipad sizing as Ihbarinv_bgn
C since Scrnc_form_rcc_hbar needs it as scratch too.
         Ihbar_bgn    = Iend - (Ineed + Ipad) + 1
         Iend         = Ihbar_bgn - 1
      Endif

      Isym_pert_bgn  = Iend - Nao_pairs + 1
      Ioff_pert_bgn  = Isym_pert_bgn - Nao_pairs
      Ijj_lrspn_bgn  = Ioff_pert_bgn - Nao_pairs * Nao_pairs
      Memleft        = Ijj_lrspn_bgn - 1

      If (Memleft .Le. 0) Then
         Length = 2*Nao_pairs + Nao_pairs * Nao_pairs
         Call Insmem("scrnc_driver",Length,Maxmem)
      Endif

      If (Rccd .OR. Drccd)  Then

C List creation (2026-08-15): moved out of Scrnc_lrspn_rcc_main.F to
C run exactly once here, before the dynamic Icount loop below re-calls
C that routine per +/-Freq pass -- see the comment in
C Scrnc_lrspn_rcc_main.F for why calling it twice corrupted things.

         Call Scrnc_form_lrspn_rcc_lists(Work(I000),Memleft,Iuhf,
     +                                   Nao_pairs)

C Dynamic (frequency-dependent) response, 2026-08-15: mirror the CCSD
C path's Icount/Itop +/-Freq loop (see scrnc_lrspn_cc_main.F). Unlike
C CCSD, ring-CC's Hbar^(-1) itself is Freq-dependent (built via
C Scrnc_form_rcc_hbar_inv as (Hbar+Freq)^(-1)), so it must be rebuilt fresh
C for each pass -- both the R-side and L-side per-pass routines below
C consume whatever Hbarinv this pass just built. The final JJ assembly
C (reading lists 373/374 and, for a dynamic run, 375/376) only happens
C ONCE, after both passes, in the new Scrnc_jj_rcc_assemble.

         Itop = 1
         If (Freq .NE. 0.0D0 .AND. Realfreq) Itop = 2

C Imaginary-frequency axis, 2026-08-16: square Freq exactly once here
C (never inside Scrnc_form_rcc_hbar_inv.F's per-irrep loop, which would
C re-square it once per irrep) -- see that routine's own header comment
C and scrnc_lrspn_cc_main.F's identical fix for the CCSD path (mirrors
C vcceh/secprp.F's own "IF(.NOT. REALFREQ) FREQ = FREQ*FREQ" convention).

         If (.NOT. Realfreq) Freq = Freq * Freq

         Call Scrnc_form_rcc_hbar(Work(Ihbar_bgn),Ineed+Ipad,Iuhf,1)
         Call Scrnc_rcc_hbar_read(Work(Ihbar_bgn),Ineed+Ipad,Iuhf,
     +                             Iloc)

         Do Icount = 1, Itop

         Call Scrnc_form_rcc_hbar_inv(Work(Ihbar_bgn),
     +                           Work(Ihbarinv_bgn),Ineed+Ipad,
     +                           Iuhf,Iloc)

          Call Scrnc_lrspn_rcc_main(Work(I000),Memleft,Iuhf,
     +                              Work(Isym_pert_bgn),
     +                              Work(Ioff_pert_bgn),Nao_pairs,
     +                              Work(Ihbarinv_bgn),Iloc,
     +                              Icount,Itop)

          Call Scrnc_jj_rcc_lrspn(Work(I000),Memleft,Iuhf,
     +                            Work(Isym_pert_bgn),
     +                            Work(Ioff_pert_bgn),Nao_pairs,
     +                            Work(Ihbarinv_bgn),Iloc,
     +                            Icount,Itop)

         Freq = -Freq

         Enddo ! Enddo for Icount (dynamic +/-Freq loop)

         Call Scrnc_jj_rcc_assemble(Work(I000),Memleft,Iuhf,
     +                              Work(Ijj_lrspn_bgn),
     +                              Work(Isym_pert_bgn),
     +                              Work(Ioff_pert_bgn),Nao_pairs)

      Else

          Call Scrnc_lrspn_cc_main(Work(I000),Memleft,Iuhf,
     +                             Work(Isym_pert_bgn),
     +                             Work(Ioff_pert_bgn),Nao_pairs)

          Call Scrnc_jj_cc_lrspn(Work(I000),Memleft,Iuhf,
     +                           Work(Ijj_lrspn_bgn),
     +                           Work(Isym_pert_bgn),
     +                           Work(Ioff_pert_bgn),Nao_pairs)
      Endif

C This block is obsolete. 
C At this point Work(Ijj_lrspn_bgn) contains the what I called the "response
C matrix". It is of the type (IJ|KL) (all alpha and occupied). We need to 
C transform this to (munu|lamdasigma) basis. First we do (IJKJL) --> (XX|KL) 
C followed (XX|KL) --> (KL|XX) --> (XX|XX). The intermediate (XX|KL)
C is kept in memory. 
C
CSSS      Call Scrnc_bcktran(Work(I000),Memleft,Iuhf,Work(Ijj_lrspn_bgn),
CSSS     +                   Nocc_pairs,Irrepx)

C AO->MO transform of Wc into the MO particle-hole basis, first step of
C the self-energy construction (route decided 2026-08-17: contour-
C deformation, see SCRNC_DEBUG_LOG.md). Always required -- not a
C diagnostic, so not 1-gated (Ajith, 2026-08-17): every
C downstream self-energy step needs Wc_mo, so this transform must run
C unconditionally, unlike Scrnc_bcktran_test.F's genuinely-optional
C validation-harness convention.
C
      Call Scrnc_wc_ao2mo_test(Work(I000),Memleft,Iuhf,
     +                         Work(Ijj_lrspn_bgn),Nao_pairs)

C Correlation-only W_c = W - v, reduced to Wc_reduced(p,q,m) =
C Wc_mo(p,m,q,m), written to "wcreduc.dat" for the external
C scrnc_sigma_driver orchestration to pick up after this pipeline pass
C completes (Ajith, 2026-08-17: do the bare-Coulomb subtraction now,
C not deferred -- see SCRNC_DEBUG_LOG.md). Always required for the
C same reason the AO2MO transform above is unconditional -- not
C 1-gated.

      Call Scrnc_sigma_prep(Work(I000),Memleft,Iuhf,
     +                       Work(Ijj_lrspn_bgn),Nao_pairs)

C Imaginary-frequency quadrature grid for the Sigma_c(i*omega) contour-
C deformation integral (Ajith, 2026-08-17: Gauss-Legendre + rational map
C to [0,Infinity), over a minimax grid -- see Scrnc_freq_grid.F's header
C comment). Pure-numerics self-check against closed-form integrals,
C genuinely optional (unlike the AO2MO transform above), so
C 1-gated like Scrnc_bcktran_test.F. NOT yet used to drive the
C actual per-frequency-point pipeline runs / Sigma_c accumulation --
C that orchestration is the next step, see SCRNC_DEBUG_LOG.md.
C
      Call Scrnc_freq_grid_test

      RETURN
      END
