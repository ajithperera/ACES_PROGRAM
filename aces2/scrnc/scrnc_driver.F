#include "flags.h"
      SUBROUTINE Scrnc_driver(Work,Maxmem,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Logical Ccsd,Mbpt,Parteom,Nodavid,Rccd,Drccd
      Dimension Work(Maxmem)
      Dimension Iloc(8)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +                Drccd

      Irrepx = 1
      Call scrnc_set_vars(Iuhf)
      Call Scrnc_getaoinf(Iuhf,Irrepx)
  
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

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Hbar elements"
      Call checkhbar(Work,Memleft,Iuhf)
#endif
      Call Scrnc_built_respn_den(Work,Maxmem,Iuhf)

C The following codes compute the linear response (lsrnc) of two Coulomb 
c operators. 

      Nao = 0
      Do Irrepr = 1, Nirrep
         Nao = Nao + Iaopop(Irrepr)
      Enddo

      Nao_pairs = Nao * Nao
   
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,I6)") " The number of ao pairs:", Nao_pairs
#endif 
      
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
            If (Iuhf .NE. 0) Ndim = Ndim + Irpdpd(Irrepr-1,10)
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
         Call Scrnc_rcc_hbar_inv(Work(Ihbarinv_bgn),Ineed+Ipad,
     +                           Iuhf,1,Iloc)
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

          Call Scrnc_lrspn_rcc_main(Work(I000),Memleft,Iuhf,
     +                              Work(Isym_pert_bgn),
     +                              Work(Ioff_pert_bgn),Nao_pairs,
     +                              Work(Ihbarinv_bgn),Iloc)

          Call Scrnc_jj_rcc_lrspn(Work(I000),Memleft,Iuhf,
     +                            Work(Ijj_lrspn_bgn),
     +                            Work(Isym_pert_bgn),
     +                            Work(Ioff_pert_bgn),Nao_pairs,
     +                            Work(Ihbarinv_bgn),Iloc)

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

      RETURN
      END
