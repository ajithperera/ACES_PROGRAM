










      Subroutine rcl_hbar(Work,Maxcor,Iuhf,Solve_4lambda)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Maxcor)

      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,QCISD,
     +        UCC,CC2,RCCD,DRCCD
      LOGICAL Get_rpa 
      LOGICAL Solve_4lambda 
      LOGICAL EOM,CIS,RPA,EOM_SDRCCD,EOM_SRCCD,EOM_SFDRCCD,EOM_SFRCCD,
     &        EOM_SDXRCCD,EOM_SFDXRCCD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/EXCITE/EOM,CIS,RPA,EOM_SRCCD,EOM_SDRCCD,EOM_SFRCCD,
     &              EOM_SFDRCCD,EOM_SDXRCCD,EOM_SFDXRCCD

C This routine construct the Hbar for ring-CC or direct CC methods. 
C At this stage, we assume that both T2 and L2 vectors corresponding
C to each method is available on the disc (list 44-6,144-46).

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

      
      Call Rcl_setmet(.True.)

C If we all care about getting the RPA using EOM-RCCD then set get_rpa
C to true. This will limit Hbar(a,b) and Hbar(i,j) to their lead term.
C (I am not sure what is the point since the original formulation of
C RPA was  EOM like; May be it is usefull to see RPA as excitaions
C from aring-CCD ground state). 
C 
      Call Rcl_setlst(Work,Maxcor,Iuhf)

      Get_rpa = .False. 
      IF (Rpa .OR. Cis .OR. Eom_sfrccd .OR. Eom_sfdrccd) 
     +   Get_rpa = .True. 

      Write(6,*)
      If (Rpa .or. Eom_sfrccd .or. Eom_sfdrccd) Then
         Write(6,"(2a)") " Hbar for rCCD or drCCD is constructed",
     +                    " when Hbar(a,b) and H(i,j) is limited"
         Write(6,"(2a)")   " to Fock diagonals: this gives EOM(SF)",
     +                    "-RCCD (RPA) or EOM(SF)-DRCCD (DRPA)"
      Else if (Cis) Then
         Write(6,"(a)") " Hbar for CIS is constructed" 
      Elseif (Eom_srccd .or. Eom_sdrccd) Then
         Write(6,"(2a)") " Hbar for EOM(S)-RCCD or EOM(S)-DRCCD is",
     +                    " constructed"
      Elseif (Eom_sfdxrccd .or. Eom_sdxrccd) Then
         Write(6,"(2a)") " Hbar for EOM(S)-DXRCCD or EOM(SF)-DXRCCD is",
     +                   " constructed"
      Endif 
C These Hbar terms originaly used only for EOM(SF)-CCD (non RPA) methods
C were previously always built regardless of Get_rpa (see git history) --
C restored the original guard: skip this amplitude correction for pure
C RPA/DRPA, where Hbar(a,b)/Hbar(i,j) must stay at the bare Fock diagonal.

      If (.Not. Get_rpa) Then
         Scale = 0.50D0
         Write(6,*)
         Write(6,"(a,a)") " Hbar(i,j) and Hbar(a,b) including",
     +                    " -1/2T2(mn,ae)<mn||be> -> Hbar(a,b) and"
         Write(6,"(a,a)") " 1/2T2(im,ef)<jm||ef> -> Hbar(i,j) is",
     +                    " built.`"
         Write(6,*)

C COnstruct W(ab,ij) from W(a<b,i<j) and save

         Call Rcl_expanded_abij(WORK(INEXT),MAXCOR,IUHF)

C Construct Hbar(a,b) = f(a,b) - 1/2T2(mn,ae)<mn||be>

         Call Rcl_form_hbar_ab(Work,Maxcor,Iuhf,Scale)

C Construct Hbar(j,i) = f(j,i) + 1/2T2(im,ef)<jm||ef>

         Call Rcl_form_hbar_ji(Work,Maxcor,Iuhf,Scale)
      Endif

C Construc Hbar(aj,ib) = <aj||ib> - T(im,ea)<mj||eb>
C for various ring methods.

      If (.NOT. Solve_4lambda) Then
         If (Drccd) Then
            If (Eom_sdxrccd .OR. Eom_sfdxrccd) Then 
               Do Ilist = 54,59,2-Iuhf
                  Call Zersym(Work,Ilist)
               Enddo 
               C4 = -1.0D0
               If (Iuhf .Ne. 0) Then
                  Call Dxrcl_dwmbej(Work,Maxcor,"AAAA",Iuhf,C4)
                  Call Dxrcl_dwmbej(Work,Maxcor,"BBBB",Iuhf,C4)
                  Call Dxrcl_dwmbej(Work,Maxcor,"BABA",Iuhf,C4)
                Endif 
                Call Dxrcl_dwmbej(Work,Maxcor,"ABAB",Iuhf,C4)
            Else 
               CALL Drcl_dwmbej(Work,Maxcor,Iuhf)
            Endif 
         Elseif (Rccd) Then
            If (Iuhf .Eq.0) Then
               CALL Rcl_dwmbej_r(Work,Maxcor,Iuhf)
            Else
               CALL Rcl_dwmbej_u(Work,Maxcor,Iuhf)
            Endif
         Endif 
      Endif 

C Hbar(ij,ab) = <ij||ab> and and Hbar(mb,ej) are already on lists (54-59)
C There are no other terms for ring-CC hbar.To be on the safe
C side, lets zero out the rest of the Hbar list. This is to 
C aoid the existing EOM code getting confused.
      
      Call rcl_set_hbar2zero(Work,Maxor,Iuhf)

C SCRNC_DEBUG_LOG.md step 2/3: independent ringcc-side Hbar
C diagonalization self-check. This is the actually-exercised path for
C excite=EOM(Sf)-DRCCD runs (RCL_DRIVER's own call site only fires for
C Solve_4lambda=.True., which this pathway never sets) -- Hbar lists
C 54-59 are fully built at this point, so diagonalize them now to get
C ringcc's own independent EOM(Sf)-DRCCD eigenvalues.
C
C Rcl_heffxc (called inside Rcl_check_rcchbar) destructively overwrites
C the same same-spin list(s) it reads from -- 54/55 for UHF, 56 for RHF
C -- adding the Fock diagonal on top, same pattern as the scrnc Hbar
C double-build bug fixed 2026-08-14. That's fine for Rcl_drvrpa's own
C diagonalization, but would corrupt those lists for whatever downstream
C production EOM step consumes rcl_hbar's raw output next. Save/restore
C around the diagnostic call so it stays side-effect-free.
      If (Iuhf .Ne. 0) Then
         Isz54 = Isymsz(Isytyp(1,54),Isytyp(1,54))
         Isz55 = Isymsz(Isytyp(1,55),Isytyp(1,55))
         I900  = 1
         I910  = I900 + Isz54*Iintfp
         I920  = I910 + Isz55*Iintfp
         If (I920 .Gt. Maxcor) Call Insmem("Rcl_hbar-chk",I920,Maxcor)
         Call Getall(Work(I900),Isz54,1,54)
         Call Getall(Work(I910),Isz55,1,55)
      Else
         Isz56 = Isymsz(Isytyp(1,56),Isytyp(1,56))
         I900  = 1
         I920  = I900 + Isz56*Iintfp
         If (I920 .Gt. Maxcor) Call Insmem("Rcl_hbar-chk",I920,Maxcor)
         Call Getall(Work(I900),Isz56,1,56)
      Endif

      Call Rcl_check_rcchbar(Work(I920),Maxcor-I920+1,Iuhf)

      If (Iuhf .Ne. 0) Then
         Call Putall(Work(I900),Isz54,1,54)
         Call Putall(Work(I910),Isz55,1,55)
      Else
         Call Putall(Work(I900),Isz56,1,56)
         Call Checksum('rcl_hbar-restored-List-56:',Work(I900),Isz56)
      Endif

      Return
      End
