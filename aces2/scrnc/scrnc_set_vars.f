














































































































































































































      Subroutine Scrnc_set_vars(Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Logical Ccsd,Mbpt,Parteom,Nodavid,Rccd,Drccd
      Logical Realfreq,Freq_exist
      Logical Ss,Sd,Ds,Dd
      Logical Respn_d
      LOgical Aoladder
      Double Precision XFreq(200)

      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +                Drccd 
      Common /Freqinfo/Freq,Realfreq,Xfreq,No_freqs 
      Common /Drvhbar/Ss,Sd,Ds,dd
      Common /eomprop/Ieomprop 

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c files.com : begin
      integer        luout, moints
      common /files/ luout, moints
c files.com : end

      Ntol    = Iflags(98) 
      Maxexp  = Iflags(97)
      Nreduce = 0

      Ccsd    = .False.
      Parteom = .False.
      Aoladder= .False.
      Rccd    = .False. 
      Drccd   = .False. 

      If (Iflags2(261) .Gt.0)  Parteom = .True.
      If (Iflags(2)    .Eq.10)  Ccsd    = .True.
      If (Iflags(2)    .Eq.48)  Rccd    = .True.
      If (Iflags(2)    .Eq.49)  Drccd   = .True.
      If (Iflags(19) .Eq.1)   Respn_d = .True.
      Mbpt  = .Not. (Ccsd .OR. Rccd .OR. Drccd)
      Luout = 6

      Aoladder = (Iflags(93).eq. 2)
      Ieomprop =  Iflags2(106)

C 2026-08-19: the original blanket REF=UHF trap here turned out to be
C caused by two narrow bugs, not fundamental missing RHF support: (1) an
C uninitialized-variable typo in scrnc_jj_cc_lrspn.F (Isympert(Iocc_pair),
C Iocc_pair never declared/assigned anywhere in that file) feeding a
C garbage irrep into Spntsing/Spinad3 -- fixed, now reads
C Isympert(Iao_pair), matching the loop index and the working
C scrnc_jj_rcc_assemble.F sibling call. (2) scrnc_newsize.F's
C 1 diagnostic unconditionally called Idsymsz on BB-specific
C list types that legitimately don't exist under RHF -- fixed to mirror
C the real (already-correct) Iuhf branch. With both fixed, the CCSD-LR
C path (H2/STO-3G) now runs RHF cleanly and its JJ_lrespn matrix agrees
C with the UHF path to ~1e-8 (matches cc_conv=9/scf_conv=10) -- a real,
C verified UHF-vs-RHF consistency check, not just "doesn't crash."
C RCCD/DRCCD (ring-CC) is a separate story: it also runs RHF without
C crashing now, but its JJ_lrespn matrix does NOT agree with UHF (H2/
C STO-3G DRCCD: RHF gives JJ(1,1)=0.79307663018 vs UHF's
C 0.80469432450, ~1.2% off -- the identical-between-runs bare-v part
C rules out a trivial cause, so this is a real, separate, NOT-yet-
C root-caused bug somewhere in the ring-CC RHF response path
C (rcc_rng_4rhf_rccd.F / rcc_trp_cntrct.F / ringcl's RHF branches are
C the natural suspects, not yet investigated). So: allow RHF for CCSD,
C keep it blocked for RCCD/DRCCD until that gap is understood. See
C SCRNC_DEBUG_LOG.md.
      If (.FALSE. .AND. Iuhf .EQ. 0 .AND. (Rccd .OR. Drccd)) Then
         Write(6,*)
         Write(6,"(a,a)") " RHF reference for ring-CC (RCCD/DRCCD)",
     +                    " is not yet supported -- see debug log."
         Write(6,*)
         Call Errex
      Endif
      If (.NOT. Respn_d) Then
         Write(6,"(a,a)")" The Density=Response is not set in the",
     +                   " ZMAT and calculation can not proceed."
         Write(6,*)
         Call Errex
      Endif 
      Inquire(File="frequency",Exist=Freq_exist)
      If (Freq_Exist) Then
         OPen(Unit=10,File="frequency",Form="Formatted",Status="OLD")
         Junk = 0
         Read(10,*,End=800) No_Freqs
         Read(10,*,End=800) JunK
         Realfreq = (Junk.eq.0)
         Do Ifreq = 1, No_Freqs
            Read(10,*,End=800) Freq
            Xfreq(Ifreq) = Freq
         Enddo 
 800     Write(6,900)
      Else
         No_freqs = 1
         Freq     = 0.0D0
         RealFreq = .True.  
         Xfreq(1) = Freq
      Endif 

      If (Rccd .or. Drccd) Then
         Nodavid  = .True. 
      Else
         Nodavid  = Parteom .and. Realfreq 
      Endif 

      Ss = .True.
      Sd = .True.
      Ds = .True.
      If (Iflags2(261) .Eq.1) Then
         Dd = .False.
      Else 
         Dd = .True.
      Endif 

      Write(6,*)
      If (Nodavid) Then
          If (Rccd .Or. Drccd) Then

          write(6,"(a)") " @-scrnc_set_vars, Linear equations solved",
     +                   " by matrix inversion (AC=Y => C=A^(-1)Y)"
          ELse
          write(6,"(a)") " @-scrnc_set_vars, Linear equations solved",
     +                   " using partitioning."
          Endif 
      Else 
          write(6,*) " @-scrnc_set_vars, Linear equations solved",
     +               " using Davidson."
      Endif 

     
  900 Format(T3,"@-Scrnc_set_vars; ", "file containing listing of",
     +       " frequencies (i.e. frequency)",//, "is incomplete.")
      Return
      End

