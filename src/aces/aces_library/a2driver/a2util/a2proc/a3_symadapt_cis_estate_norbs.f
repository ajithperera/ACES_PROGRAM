      subroutine a3_symadapt_cis_estate_norbs(ener,iocc,orb,dens,nlorb,
     &                                    edens,Oed2AScal, Ioed2Aord,
     &                                    tmp1,scr,nao,nbas,maxcor,
     &                                    iuhf,iroot,Nbfirr, Nirrep,
     &                                    O_MOS, T_MOS, Eneg, Cis_exes,
     &                                    Smult,Nroots,Root_count,Ipick)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)



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



c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)


c ***NOTE*** This is a genuine (though not serious) limit on what Aces3 can do.
c     12 => s,p,d,f,g,h,i,j,k,l,m,n
      integer maxangshell
      parameter (maxangshell=12)


C
      parameter (one=1.0D0)
      parameter (zilch=0.0D0)
      parameter (DENS_THRESH=1.0D-08)
C
      double precision ener(nao),orb(nao*nao),scr(maxcor)
      double precision dens(nao*nao),edens(nao,nao),nlorb(nao*nao)
      double precision nelec, tmp1(Nao,Nao)
      
      Dimension Oed2AScale(Nao), Ioed2Aorder(Nao)
      Dimension Nprim_shell(Maxangshell*Mxatms)
      Dimension Orig_nprim_shell(Maxangshell*Mxatms)
      Integer   Reorder_Shell(Maxangshell*Mxatms)
C
      character*2 iroot
      character*6 denstype(2)
      character*8 denstyper_aa(24)
      character*8 denstypel_aa(24)
      character*8 denstyper_bb(24)
      character*8 denstypel_bb(24)
      character*8 cscforb(2), Transform(4), Scfvecs(2)
      character*8 cexxcoef(2)
      character*8 String
      character*6 Rhf_string 
      character*6 cnumdrop(2)
      Character*8 Dump_string
      character*5 sptype(2)

      double precision  iocc(nao), Eneg(Nao)
      double precision  Cis_exes(24),Smult(24)
      integer  idrppop(8),idrpvrt(8)
      integer  nocc(8,2), Nbfirr(8) 
      integer  Root_count 
      integer  O_MOS(Nao), T_MOS(Nao)
C
      data denstype /'REOMDN','LEOMDN'/
      data denstyper_aa /'REOMDN1A','REOMDN2A','REOMDN3A','REOMDN4A',
     +                   'REOMDN5A','REOMDN6A','REOMDN7A','REOMDN8A',
     +                   'REOMDN9A','REOMD10A','REOMD11A','REOMD12A',
     +                   'REOMD13A','REOMD14A','REOMD15A','REOMD16A',
     +                   'REOMD17A','REOMD18A','REOMD19A','REOMD20A',
     +                   'REOMD21A','REOMD22A','REOMD23A','REOMD24A'/
      data denstyper_bb /'REOMDN1B','REOMDN2B','REOMDN3B','REOMDN4B',
     +                   'REOMDN5B','REOMDN6B','REOMDN7B','REOMDN8B',
     +                   'REOMDN9B','REOMD10B','REOMD11B','REOMD12B',
     +                   'REOMD13B','REOMD14B','REOMD15B','REOMD16B',
     +                   'REOMD17B','REOMD18B','REOMD19B','REOMD20B',
     +                   'REOMD21B','REOMD22B','REOMD23B','REOMD24B'/
      data denstypel_aa /'LEOMDN1A','LEOMDN2A','LEOMDN3A','LEOMDN4A',
     +                   'LEOMDN5A','LEOMDN6A','LEOMDN7A','LEOMDN8A',
     +                   'LEOMDN9A','LEOMD10A','LEOMD11A','LEOMD12A',
     +                   'LEOMD13A','LEOMD14A','LEOMD15A','LEOMD16A',
     +                   'LEOMD17A','LEOMD18A','LEOMD19A','LEOMD20A',
     +                   'LEOMD21A','LEOMD22A','LEOMD23A','LEOMD24A'/
      data denstypel_bb /'LEOMDN1B','LEOMDN2B','LEOMDN3B','LEOMDN4B',
     +                   'LEOMDN5B','LEOMDN6B','LEOMDN7B','LEOMDN8B',
     +                   'LEOMDN9B','LEOMD10B','LEOMD11B','LEOMD12B',
     +                   'LEOMD13B','LEOMD14B','LEOMD15B','LEOMD16B',
     +                   'LEOMD17B','LEOMD18B','LEOMD19B','LEOMD20B',
     +                   'LEOMD21B','LEOMD22B','LEOMD23B','LEOMD24B'/

      data sptype   /'Alpha','Beta '/
      data Transform /"OOTRANSA", "VVTRNASA", "OOTRANSB", "VVTRANSB"/
      data Scfvecs /"SCFEVCA0", "SCFEVCB0"/
C
      call aces_com_info
      call aces_com_syminf
      call aces_com_sym
C
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c the mo ordering is important for getting the occupations right
c
c depending on the place in the program (if the user is not using xaces2)
c it may be either scf or correlated
c
c    so            mo                        mo
c ao CMP2ZMAT   so cscforb(ispin)    =    ao AOBASMOS
c
c the so ordering is zmat
c the mo ordering is correlated if calc is greater than scf 
c or if this is a vibrational calculation
c
      Call Getrec(20, "JOBARC", "NSHELLS" , 1, nshells)
      Call Getrec(20, "JOBARC", "NPRMSHEL", nshells, Nprim_shell)
      Call Getrec(20, "JOBARC", "BNPAKORD", nshells, Reorder_Shell)
C
      Call Getrec(20, "JOBARC", "ERD2A2CS", Nbas*Iintfp, Oed2AScale)
      Call Getrec(20, "JOBARC", "ERDORDER", Nbas, Ioed2Aorder)
C
      i000=1
      i010=i000+nbas*nbas
      i020=i010+nbas*nbas
      i030=i020+nbas*nbas
      i040=i030+nbas
      iend=i040+nbas

      if(Iend.gt.maxcor)

     &  call insmem('a3_symadapt_estate_norbs',iend,maxcor)

      If (Nroot .Gt. 12) Then
        Write(6,"(a,a)") "More than 12 roots is"
     &                   " not allowed for this analysis."
        Call Errex
      Endif 

      N_o_mos = 0
      N_t_mos = 0
      Nrootsp1= Nroots+1 

C Track the occ-occ and virtual orbitals that at least have 0.5
C contribution to the transition. Note that cutoff value is simply a choice
C that I made and also those  above the cutoff are the only orbitals that are 
C printed. 

      If (Iuhf .EQ. 0) Cutoff = 0.20D0
      If (Iuhf .EQ. 1) Cutoff = 0.20D0

      If (Ipick .EQ. 1) Then
         Do i=1, Nroots
            Smult(i) = 1.0
         Enddo
      else if (Ipick .EQ. Nrootsp1) Then 
         Do i=1,Nroots 
            Ioff = i+Nroots 
            Smult(Ioff) = 2.0
         Enddo
         Call Dcopy(Nroots,CIS_exes,1, CIS_exes(Nroots+1),1) 
      Endif 

      If (Ipick .eq. 1) then
      Write(6,"(a)") " The singlet CIS excitation energies"
      Write(6,"(6(1x,F12.6))") (CIS_exes(i),i=1,Nroots)
      Write(6,"(a)") " The spin mutiplicities of the CIS states"
      Write(6,"(6(1x,F12.6))") (Smult(i),i=1,Nroots)
      Else 
      Write(6,"(a)") " The triplet CIS excitation energies"
      Write(6,"(6(1x,F12.6))") (CIS_exes(i+Nroots),i=1,Nroots)
      Write(6,"(a)") " The spin mutiplicities of the CIS states"
      Write(6,"(6(1x,F12.6))") (Smult(i+Nroots),i=1,Nroots)
      Endif 
      do ispin=1,iuhf+1

         Do Idens = 1, 2

C Generate the occupation numbers for each irrep based on eigen
C values (SCF) and the number of basis functions per irrep.
C
        If (Ispin .EQ. 1) Call Getrec(20, "JOBARC", "SCFEVLA0", 
     &                                Nbas*Iintfp, Ener)
        If (Ispin .EQ. 2) Call Getrec(20, "JOBARC",
     &                                             "SCFEVLB0",
     &                                              Nbas*Iintfp, 
     &                                              Ener)
C
        If (ispin .EQ. 1) Then
            Call Occupy(Nirrep, Nbfirr, Nbas, ener, scr, Nocc(1,1),
     &                  1)
        Else 
            Call Occupy(Nirrep, Nbfirr, Nbas, ener, scr, Nocc(1,2),
     &                  1)
        Endif 
        If (Iuhf .EQ. 0) Call Icopy(8, Nocc(1, 1), 1, Nocc(1,2), 1)


      Noccs = 0
      Do irrep = 1, Nirrep
         Noccs = Noccs + Nocc(irrep,Ispin)
      Enddo
      Nvrts = Nbas - Noccs
C 
C Set the Norb depending on whether we do occ-occ or vrt-vrt block of
C the transition density. The current setup is to handle occ-occ block
C first (idens=1).
C
      If (Idens .Eq. 1) Then
         Norbs = Noccs
         Ioff  = 0
      Else
         Norbs = Nvrts
         Ioff  = Noccs
      Endif
C
C Get the occ-occ and vrt-vrt MO density matrices 
C
      If (Iuhf .eq. 0) Then
          Length = Norbs * Norbs * iintfp 
          RHF_String = Denstype(Idens)
          call getrec(20,'JOBARC',RHF_string//iroot,
     +                    Length,Dens)
         
      Else
          If (Ispin .Eq. 1) Then
              Length = Norbs * Norbs * iintfp 
              If (Idens .EQ. 1) String = Denstyper_aa(Root_count)
              If (Idens .EQ. 2) String = Denstypel_aa(Root_count)
              call getrec(20,'JOBARC',String,Length,Dens)
          Else if (Ispin .eq. 2) Then
              Length = Norbs * Norbs * iintfp 
              If (Idens .EQ. 1) String = Denstyper_bb(Root_count)
              If (Idens .EQ. 2) String = Denstypel_bb(Root_count)
              call getrec(20,'JOBARC',String,Length,Dens)
          Endif 
      Endif 
C
C     Symmetrize and Diagonalize MO basis density
C
        call symmet2(Dens, Norbs)
C
C Eig would have been fine except that it reorders eigenvalues and
C vectors.
CSSS        call eig (dens,Nlorb,1,Norbs,-1)

        I050 = Iend
        I060 = I050 + Norbs
        I070 = I060 + Norbs
        I080 = I070 + Norbs*Norbs
        I090 = I080 + Norbs*Norbs
        Iend = I090 + 4*Norbs 
        If (iend .gt. Maxcor) 
     +     call insmem('a3_symadapt_estate_norbs',iend,maxcor)
 
        call dgeev("N","V",Norbs,Dens,Norbs,Scr(I050),Scr(I060), 
     +              Scr(I070),Norbs,Scr(I080),Norbs,Scr(I090),
     +              4*Norbs,Ierror) 

        If (Ierr .Ne. 0) Then
            Write(6,"(a)") "Eigenvector solver failed" 
            Call Errex
        Endif 

        Call Dcopy(Norbs*Norbs,scr(I080),1,Nlorb,1)

        If (Idens .EQ. 1) call dcopy(Norbs,Scr(I050),1,iocc,1)
        If (Idens .EQ. 2) call dcopy(Norbs,Scr(I050),1,
     &                               iocc(Noccs+1),1)

        nelec = 0.0D0
        I     = 0
        If (Idens .EQ. 1) Call IZero(O_mos, Nao)
        If (Idens .EQ. 2) Call IZero(T_mos, Nao)

        DO imo = 1, Norbs
           If (Idens .eq. 1) nelec = nelec + (iocc(imo))
           If (Idens .eq. 2) nelec = nelec + (iocc(Noccs+imo))
C
           If (Idens .EQ. 1) Then
              If (Abs(iocc(imo)) .GT. Cutoff) Then
                 I = I + 1
                 O_mos(I) = Noccs - imo + 1
              Endif 
           Else if (Idens .EQ. 2) then
              If (Abs(iocc(Noccs+ imo)) .GT. Cutoff) Then
                 I = I + 1
                 T_mos(I) = imo + noccs 
              Endif
           Endif

        ENDDO

        If (Idens .EQ. 1) then
           N_o_mos = 0
           Do Imo = 1, Nao
              If (O_MOS(Imo) .Gt. 0) N_o_mos = N_o_mos + 1
           Enddo
        Elseif (Idens .EQ. 2) then
           N_t_mos = 0
           Do Imo = 1, Nao
              If (T_MOS(Imo) .Gt. 0) N_t_mos = N_t_mos + 1
           Enddo
        Endif 

        Write(6,*)
        WRITE (6,"(a)") ' Trace of diagonalized MO density matrix: '
        WRITE (6,"(F10.5)")  nelec
        Write(6,"(a)") "Eigenvalues"
        If (Idens .EQ. 1) Write(6,"(6(1x,F10.5))") 
     &                         (iocc(imo), imo=1, Norbs)
        If (Idens .EQ. 2) Write(6,"(6(1x,F10.5))") 
     &                         (iocc(imo+noccs), imo=1, Norbs)
        If (Idens .EQ. 1) Write(6,"(6(1x,I4))") (O_mos(I),I=1,Norbs)
        If (Idens .EQ. 2) Write(6,"(6(1x,I4))") (T_mos(I),I=1,Norbs)
C
C Write the OO and VV transformation matrices to JOBARC 
C
C        Call putrec(20, "JOBARC", TRANSFORM(Idens + (ISpin-1)), 
C       &            Norbs*Norbs*Iintfp, Nlorb)
C
C Do the OO and VV transformation. Retrive the SCF vectors of correct
C spin type during first iteration of the inner loop (OO block)
C
        call getrec(20,'JOBARC', SCFVECS(Ispin), Nbas*Nbas*iintfp,dens)


      Enddo 


        call getrec(20,'JOBARC', SCFVECS(Ispin), Nbas*Nbas*iintfp,
     +              scr(i020))
        call getrec(20,'JOBARC','CMP2CART',nao*nbas*iintfp,scr(i000))

        call xgemm('n','n',nao,nbas,nbas,one,scr(i000),nao,scr(i020),
     &              nbas,zilch,orb,nao)

C
C Lets write these corrsponding orbitals to the JOBARC file so 
C that the MOLDEN file can be written.
C
        If (Ispin .EQ. 1) Dump_String = "CRORBA"//iroot
        If (Ispin .EQ. 2) Dump_String = "CRORBB"//iroot
        Call Putrec(20,'JOBARC', Dump_String, nao*nbas*iintfp, Orb)

        Write(6,*)
        WRITE (6,"(a)") "The active occupied orbitals"
        Write(6,"(6(i3))") (O_MOS(i), i=1, N_t_mos)
        WRITE (6,"(a)") "The active virtuals orbitals"
        Write(6,"(6(i3))") (T_MOS(i), i=1, N_t_mos)
C 
        Do Imo = 1,  N_o_mos
           Eneg(Imo) = Ener(O_MOS(Imo))
        Enddo 
        Do Imo = 1,  N_t_mos
           Eneg(N_o_mos + Imo) = Ener(T_MOS(Imo))
        Enddo 
C
        N_ot_total = N_o_mos + N_t_mos

        Call Dcopy(Nbas*Nao, Orb, 1, Nlorb, 1)
CSSS
CSSS        Call Dcopy(N_t_mos*Nao, Orb(Noccs*Nao+1), 1, 
CSSS     &             Nlorb(N_o_mos*Nao+1), 1)
        Write(6, "(1x,a,F12.6,a)") " The Excitation Energy :" ,
     &                               CIS_Exes(Root_count), " eV"
        If (Iuhf .EQ. 0) Then
           Write(6,"(1x,a,a)") " The spin-state :"," Singlet"
        Else
           If (Int(Smult(Root_count)) .EQ. 1) Write(6,"(1x,a,a)")
     &                                  " The spin-state :"," Singlet"
           If (Int(Smult(Root_count)) .EQ. 2) Write(6,"(1x,a,a)")
     &                                  " The spin-state :"," Triplet"
        Endif
        Write(6,*) 
        Call Get_irreps_ex(Nlorb, Eneg, Scr, Imemleft*Iintfp, 
     &                     N_ot_total, Nbas, 1, Nocc, Iuhf, 
     &                     Ispin, O_MOS, T_MOS, N_o_mos, 
     &                     N_t_mos, "EXCITED")

      Enddo
C
      return
      end




