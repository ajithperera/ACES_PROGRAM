













































































































































































































      Subroutine Tdcc_driver(Work,Maxmem,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxmem)
      Logical MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      Integer AAAA_LENGTH_IJKA
      Integer Isympert(3)
      Character*8 Label(3),Work_label
      Double Precision Mu_0, Mu_tilde_0, Mu_0_dot

      Logical MU, MU_tilde, MU_dot, MU_tilde_dot



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

       Common /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2

       Data Ione /1/

       Sing   = .False.
       Icontl = Iflags(4) 
       Iconvg = 1
       Ncycle = 0
       Sing   = (Iflags(2) .gt. 9)

       MU           = .True.
       MU_tilde     = .True.
       MU_dot       = .True.
       MU_tilde_dot = .True.
C
C Unfortunately I have to use the METH common block. At the moment 
C We will set to CCSD using calc key-word. 
C
       CCSD = (Iflags(2) .eq. 10)

       If (CCSD) Then
          Coulomb = .False.
          Call Parread(iuhf)
         If (ispar) Then
           write(6,2010) paralpha
           write(6,2011) parbeta
           write(6,2012) pargamma
           write(6,2013) pardelta
           write(6,2014) parepsilon
 2010      format(' PCCSD   alpha parameter : ', F14.6)
 2011      format(' PCCSD    beta parameter : ', F14.6)
 2012      format(' PCCSD   gamma parameter : ', F14.6)
 2013      format(' PCCSD   delta parameter : ', F14.6)
 2014      format(' PCCSD epsilon parameter : ', F14.6)
          write(6,"(a,a)") ' Perform a parameterized CC HBAR',
     &                    ' calculations'
          Write(6,*)
          if (coulomb) write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
          write(6,*)
          Fae_scale    = (Paralpha - 1.0D0)
          Fmi_scale    = (Parbeta  - 1.0D0)
          Wmnij_scale  = Pargamma
          Wmbej_scale  = Pardelta
          Gae_scale    = Paralpha  
          Gmi_scale    = Parbeta  
        Else
          write(6,*) '   Perform a regular CC HBAR calculations'
          write(6,*)
          Fae_scale    = 0.0D0
          Fmi_scale    = 0.0D0
          Wmnij_scale  = 1.0D0
          Wmbej_scale  = 1.0D0
          Gae_scale    = 1.0D0
          Gmi_scale    = 1.0D0
         Endif 
      Endif

      Write(6,*) "Checking T, L and Hbar lists"
      Call check_t(work,Memleft,Iuhf)
      Call check_l(work,Memleft,Iuhf)
      Call checkhbar(Work,Memleft,Iuhf)
C
C Obtain the symmetry of dipole moment integrals.
C
  
      Nbfns   = Nocco(1) + Nvrto(1)
      Call Getrec(20, "JOBARC","NBASTOT ", Ione, Naobfns)

      Call Tdcc_get_pert_type(Work,Memleft,Nbfns,Label,Isympert,
     +                        Iuhf)

      write(6,*)
      write(6,"(a,1x,I4,1x,I4)") "The Nbfns and Naobfns:",
     +                            Nbfns,Naobfns

      Call Tdcc_form_lists(Work,Memleft,Iuhf)

C Loop over the 3-components of the diople integrals 

      Do Ixyz = 1, 3

         Irrepx     = Isympert(Ixyz)
         Ipert      = Ixyz
         Work_label = Label(Ixyz)

         Call Tdcc_reset_lists(Work,Memleft,Iuhf,Irrepx)

         Length   = Naobfns * (Naobfns+1)/2
         Length21 = Nbfns * Nbfns
         Length22 = Nbfns * Naobfns
         Lenoo    = Irpdpd(Irrepx,21) + Iuhf * Irpdpd(Irrepx,22)
         Lenvv    = Irpdpd(Irrepx,19) + Iuhf * Irpdpd(Irrepx,20)
         Lenvo    = Irpdpd(Irrepx,9)  + Iuhf * Irpdpd(Irrepx,10)
         Lenvo2_aa= Irpdpd(Irrepx,9)  * Irpdpd(Irrepx,9) 
         Lenvo2_bb= Irpdpd(Irrepx,10) * Irpdpd(Irrepx,10)
         Lenvo2_ab= Irpdpd(Irrepx,11) * Irpdpd(Irrepx,11)

        Write(6,"(a,a)") "Length,Length21,Length22,Lenoo,Lenvv,Lenv0"
     +                   ",Lenvo2_aa,Lenvo2_bb,Lenvo2_ab: "
        Write(6,"(9(1x,I5))")Length,Length21,Length22,Lenoo,Lenvv,
     +                       Lenv0,Lenvo2_aa, Lenvo2_bb, Lenvo2_ab
        Write(6,*)
         Ibgn = Ione
         I000 = Ibgn
         I010 = I000 + Lenoo 
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Length 
         I050 = I040 + Length21
         I060 = I050 + Length22
         I070 = I060 + Length22
         Iend = I070

         Memleft =  Maxmem - Iend
         If (Iend .Gt. Maxmem) Call Insmem("@-Tdvee_driver",Iend,
     +                                      Length)
         Call Tdcc_prep_dipole_ints(Work(I030),Work(I040),Work(I050),
     +                              Work(I060),Work(I000),Work(I010),
     +                              Work(I020),Work(Iend),
     +                              Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
     +                              Naobfns,Iuhf,Ipert,Work_label,
     +                              Irrepx)
C Everything is setup to do the 9 Equations given in Nascimento and 
C Deprince III, JCTC, 12, 5834-5840, (2016). In my opinion this is 
C the kind of work that our group should have been doing instead of
C applied statistics (i.e. Fitting DFT parameters....)
 
         If (MU) Then
            Ibgn   = I030
            I040   = I030 + Lenvo
            Iend   = I040
            Memleft= Maxmem - Iend

            Call Tdcc_form_mu_0_0(Work(Iend),Memleft,Work(I000),
     +                            Work(I010),Work(I020),Irrepx,
     +                            Lenoo,Lenvv,Lenvo,Iuhf,Mu_0)

            Call Tdcc_form_mu_0_S(Work(Iend),Memleft,Work(I000),
     +                            Work(I010),Work(I020),Work(I030),
     +                            Irrepx,Lenoo,Lenvv,Lenvo,Iuhf)
         Endif 
           
         If (MU_tilde) Then 
        
            Ibgn   = I030
            I040   = I030 + Lenvo
            Iend   = I040

            Memleft= Maxmem - Iend
            If (Iend .Gt. Memleft) Call Insmem("@-Tdvee_driver",Iend,
     +                                         Memleft)

            Call Tdcc_form_mutilde_0_0(Work(Iend),Memleft,
     +                                 Work(I000),Work(I010),
     +                                 Work(I020),Irrepx,Lenoo,
     +                                 Lenvv,Lenvo,Iuhf,Mu_tilde_0)

            Call Tdcc_form_mutilde_0_s(Work(Iend),Memleft,
     +                                 Work(I000),Work(I010),
     +                                 Work(I020),Work(I030),
     +                                 Irrepx,Lenoo,Lenvv,Lenvo,
     +                                 Iuhf,Mu_0)
      call tdcc_mutilde_0_s_debug(Work(I000),Work(I010),
     +                            Work(I020),Work(I030),
     +                            Work(Iend),
     +                            Memleft,Iuhf,
     +                            Irrepx,2,Mu_0)
           
            If (Iuhf .EQ. 0) Then
               Call Tdcc_form_mutilde_0_d_rhf(Work(Iend),Memleft,
     +                                        Work(I000),Work(I010),
     +                                        Work(I020),
     +                                        Irrepx,Lenoo,Lenvv,Lenvo,
     +                                        Iuhf,Mu_0)
      call tdcc_mutilde_0_d_rhf_debug(Work(I000),Work(I010),
     +                                     Work(I020),
     +                                     Work(Iend),
     +                                     Memleft,Iuhf,
     +                                     Irrepx,2,Mu_0)
            Else
               call Tdcc_form_mutilde_0_d_uhf(Work(Iend),Memleft,
     +                                        Work(I000),Work(I010), 
     +                                        Work(I020),
     +                                        Irrepx,Lenoo,Lenvv,Lenvo,
     +                                        Iuhf,Mu_0)
      call tdcc_mutilde_0_d_uhf_debug(Work(I000),Work(I010),
     +                                Work(I020),
     +                                Work(Iend),
     +                                Memleft,Iuhf,
     +                                Irrepx,2,Mu_0)
            Endif 

         Endif 

         If (MU_dot) Then 

            Call Tdcc_form_mu_dot_0(Work,Memleft,Irrepx,Iuhf,Mu_0_dot)

            Call Tdcc_form_mu_dot_s(Work,Memleft,Irrepx,Iuhf)

            Call Tdcc_form_mu_dot_d(Work,Memleft,Irrepx,Iuhf)

         Endif 

         If (MU_tilde_dot) Then   
   
            Call Tdcc_form_mutilde_dot_s(Work,Memleft,Irrepx,Iuhf)

            Call Tdcc_form_mutilde_dot_d(Work,Memleft,Irrepx,Iuhf)
   
         Endif

C Built the overlaps (otherwise known as the aut-correlation function)

  40     Continue 

         Call Tdcc_built_ac_func(Work,Memleft,Irrepx,Iuhf)

      Enddo 

      Return
      End
