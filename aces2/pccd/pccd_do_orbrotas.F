#include "flags.h"
      Subroutine Pccd_do_orbrots(Work,Maxcor,Iuhf,IGrad_calc,IHess_calc,
     +                           Scale,OOmicroItr,OOmacroItr,OOtotalItr)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Grdcnv
      Logical pCCD,CCD,LCCD
      Logical OneP_ONLY,Hdump,Onep_save,Ohess

      Common /CALC/pCCD,CCD,LCCD
      Common /ORBR_HESS/ONEP_ONLY

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"

      Data Ione,Itwo,Izero,One,Half /1,2,0,1.0D0,0.50D0/

      Grdcnv = .False.
      Ncycle = 0
      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)
C RFO_FAIL: unconverged amplitudes, no densities; Pccd_rfo_step only
C rejects the step.
      Call Getrec(0,"JOBARC","RFO_FAIL",Length,Junk)
      Ifail = 0
      If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_FAIL",1,Ifail)

      If (Ifail .Eq. 1) Then
         Continue
      Else If (Iuhf .EQ. 0 .AND. Idrop .EQ. 0) Then
C Analytic frozen-amplitude orbital Hessian (pCCD): reads Gamma lists
C 123/125 before Pccd_form_htau_pq reorders them. Debug switch
C pccd_ohess.in (gradient check and pccd_ohess.dat).
         Inquire(File="pccd_ohess.in",Exist=Ohess)
C OO_METHOD=RFO unless pccd_rfo_h0.in says CONST (not at the displaced
C points of the stability check).
         Ih0 = 0
         If (Iflags2(h_IFLAGS2_oo_method) .Eq. 6) Then
            Call Pccd_rfo_h0mode(Ih0)
            Call Getrec(0,"JOBARC","RFO_PHAS",Length,Junk)
            Iphase = 0
            If (Length .Gt. 0) Call Getrec(20,"JOBARC","RFO_PHAS",1,
     +                                     Iphase)
            If (Iphase .Eq. 1) Ih0 = 0
         Endif
         Iok = 0
         If (pCCD .And. (Ohess .Or. Ih0 .Ge. 1))
     +      Call Pccd_ohess(Work,Maxcor,Iok)
         Call Putrec(20,"JOBARC","OHESOK  ",1,Iok)
C Debug (file pccd_hess.in): dump the two-particle density and all
C other lists 1-250 as they are (GA, AL), then form the gradient and
C the analytic orbital Hessian lists once, in the production order
C (Pccd_form_htau_pq reorders Gamma lists 123/125 in place, so it must
C not run twice), and dump the Hessian lists (FU).
         Inquire(File="pccd_hess.in",Exist=Hdump)
         If (Hdump) Then
            Call Pccd_hess_dump(Work,Maxcor,"GA  ")
            Call Pccd_hess_dump(Work,Maxcor,"AL  ")
            Call Pccd_hess_dump(Work,Maxcor,"ZERO")
            Onep_save = ONEP_only
            ONEP_only = .False.
            Mxc = Maxcor
            Call Pccd_form_htau_pqrs(Work,Mxc,Iuhf)
            Mxc = Maxcor
            Call Pccd_form_htau_pq(Work,Mxc,Iuhf,Igrad_calc,1,Scale)
            Call Pccd_hess_dump(Work,Maxcor,"FU  ")
            ONEP_only = Onep_save
         Else
         If(.Not. ONEP_only) Call Pccd_form_htau_pqrs(Work,Maxcor,Iuhf)
C OO_METHOD=RFO does not use the (unfinished, undebugged) analytic
C Hessian lists that Pccd_form_htau_1D_pqrs builds after the gradient.
         Ihess = IHess_calc
         If (Iflags2(h_IFLAGS2_oo_method) .Eq. 6) Ihess = 0
         Call Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          Ihess,Scale)
         Endif
         If (pCCD .And. Ohess) Call Pccd_ohess_chk(Work,Maxcor)
      Else

C Uhf only or UHF and RHF for forzen core orbitals.

         Call Pccd_form_uhtau(Work,Maxcor,Iuhf)
      Endif 

      Call GetSDInfo(nocca,noccb,nvrta,nvrtb,nbas,nvecDim,Naobfns)

      If (Iflags2(h_IFLAGS2_oo_method) .EQ. 1) then ! BB Steep. Desc.
         call Drive_BB_SD(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                    nocca,nvrta,nvecDim)
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 2) then !BFGS
        print*,'Have not include L-BFGS yet; to do'
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 3) then !FullNR
        print*,'Have not included FullNR yet; to do'
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 4) then !AMSgrad
          call Drive_AMSgrad(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                        nocca,nvrta,nvecDim)
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 6) then !RFO
          Tol  = 10.0D0**(-Iflags(h_IFLAGS_bruck_conv))
          Etol = 10.0D0**(1-Iflags(h_IFLAGS_cc_conv))
          Call Pccd_rfo_step(Work,Maxcor,Nbas,Nocca,Iuhf,Grdcnv,Tol,
     +                       Etol,Grd_max,Grd_rms,Ncycle)
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 5) then !AMSgrad
          Call Getrec(-20,'JOBARC','ORBOPITR',Ione,Ncycle)
          Call Pccd_rotgrd(Work,Maxcor,Nbas,Nocca,Noccb,Nvrta,Nvrtb,
     +                     Iuhf,Grdcnv,Tol,Grd_max,Grd_rms,Ncycle)
      endif
C 
      If (Grdcnv .And. .Not. Ncycle .Eq. 0) Call Pccd_finalize(Work,
     +                                      Maxcor,Iuhf,Tol,Grd_max,
     +                                      Grd_rms)
      If (IFLAGS(h_IFLAGS_reference).GT.0) 
     +IFLAGS(h_IFLAGS_reference) = Ione
      IFLAGS(h_IFLAGS_scf_maxcyc)= Izero
      IFLAGS(h_IFLAGS_nonhf)     = Ione
C Do not turn on non-HF key-word sinply because this indirectly add singles 
      IFLAGS(h_IFLAGS_qrhf_orb)  =Izero
      IFLAGS(h_IFLAGS_guess)     =Izero

      Call Putrec(20,'JOBARC','IFLAGS  ',100,IFLAGS)

      Return
      End 
 
