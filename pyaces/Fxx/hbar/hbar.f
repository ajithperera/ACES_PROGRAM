







































































































































































































      Subroutine Hbar_(Work, Maxcor, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Maxcor)

      Logical  LCCD, LCCSD,  CCSD, MBPT2, CC2, HBAR_4LCCSD
      Logical  HBAR_4DCCSD, HBAR_4DCCD
      Logical  CCD, MBPT3,M4DQ,M4SDQ,M4SDTQ,QCISD
      Logical  DCCSD, DCCD
      Logical  Do_tau,Nonhf,Adc2,Adcc2
      Logical  Term1,Term2,Term3,Term4,Term5,Term6

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


      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2

C Here, Formwl,wtwtw manupilate Hbar(mb,ej).
C Formw4 and formw5 manupilate Hbar(ia,jk) and Hbar(ab,ci).
C Hbarijka form the Hbar(ij,ka).  Hbar(mn,ij), Hbar(a,e), Hbar(m,i)
C and Hbar(m,e) is formed in the genint call from the CC code. 
c The fixfbar is needed to CCSD, CC2, CCD, MBPT(2) (to remove T*F 
C extra contrbution). 

C This tdee_init list call is to test time dependent CC without 
C forming Hbar (by setting T1 and T2 to zero). It should be used
C only in that context.
    
CSSS      call tdee_init_lists(Work,Maxcor,Iuhf)

      Call Getrec(0, "JOBARC","LAMBDA  ",Length,Junk)
      Call Putrec(20,"JOBARC","HBAR    ",1,Junk)

      If (iflags(91) .GT. 0 .And. Length .GT. 0) Then
         Write(6,*) 
         Write(6,"(a)")   " !!!    Warning   !!!"
         Write(6,"(a,a)") " Excited state derivative or response",
     &                    " properties requested and therefore",
     &                    " lambda vectors are needed!"
         Write(6,"(a,a)") " We can not run xlambda and xhbar instead",
     &                    " of that choose xhbar and xlhbar"
         Write(6,"(a)")   " sequence."
         Write(6,*)
         Call Errex
      Endif 
      
      Coulomb = .False. 
      Call Parread(iuhf)
      If (ispar) Then
         Write(6,*)
         write(6,"(a,a)") ' Perform a parameterized CC HBAR',
     &                    ' calculations'
         Write(6,*)
         write(6,2010) paralpha
         write(6,2011) parbeta
         write(6,2012) pargamma
         write(6,2013) pardelta
         write(6,2014) parepsilon
 2010    format(' PCCSD   alpha parameter : ', F14.6)
 2011    format(' PCCSD    beta parameter : ', F14.6)
 2012    format(' PCCSD   gamma parameter : ', F14.6)
 2013    format(' PCCSD   delta parameter : ', F14.6)
 2014    format(' PCCSD epsilon parameter : ', F14.6)
         if (coulomb) Write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
         write(6,*)
         Fae_scale    = (Paralpha - 1.0D0)
         Fmi_scale    = (Parbeta  - 1.0D0)
         Wmnij_scale  = Pargamma
         Wmbej_scale  = Pardelta
         Gae_scale    = Paralpha
         Gmi_scale    = Parbeta  
      Else
         write(6,"(a)") ' Perform a canonical CC HBAR calculations'
         Fae_scale    = 0.0D0
         Fmi_scale    = 0.0D0
         Wmnij_scale  = 1.0D0
         Wmbej_scale  = 1.0D0
         Gae_scale    = 1.0D0
         Gmi_scale    = 1.0D0
      Endif
C
C It is much cleaner to built Hbar based on individual theoretical
C model from the outset. 
C
      CCSD  = .FALSE.
       CCD  = .FALSE.
       CC2  = .FALSE.
      MBPT2 = .FALSE.
      LCCSD = .FALSE.
      LCCD  = .FALSE.
      DCCSD = .FALSE.
      DCCD  = .FALSE.

      If (iflags(2) .EQ. 5) Then
         LCCD        = .TRUE.
         HBAR_4LCCSD = .FALSE.
      Elseif (iflags(2) .EQ. 8) Then
         CCD         = .TRUE.
         HBAR_4LCCSD = .FALSE.
CSSS#ifdef _DCC_FLAG
         If (Ispar) Then
            CCD         = .FALSE.
            DCCD        = .TRUE.
            HBAR_4LCCSD = .FALSE.
            HBAR_4DCCD  = .TRUE. 
         Endif 
CSSS#endif
      Elseif (iflags(2) .EQ. 6) Then
         LCCSD       = .TRUE.
         HBAR_4LCCSD = .TRUE.
      Elseif (iflags(2) .EQ. 10) Then
         CCSD        = .TRUE.
         HBAR_4LCCSD = .FALSE.
CSSS#ifdef _DCC_FLAG         
         If (Ispar) Then
            CCSD        = .FALSE.
            DCCSD       = .TRUE.
            HBAR_4LCCSD = .FALSE.
            HBAR_4DCCSD = .TRUE.
         Endif 
CSSSS#endif 
          
      Elseif (iflags(2) .EQ. 1) Then
         MBPT2       = .TRUE.
         HBAR_4LCCSD = .FALSE.
      Elseif (iflags(2) .EQ. 47) Then
         CC2         = .TRUE.
         HBAR_4LCCSD = .FALSE.
      Endif 
 
      Nonhf  = (iflags(38) .ne. 0)

      If ((iflags2(117) .EQ. 12) .AND. MBPT2)
     *   ADC2  = .TRUE.
      If ((iflags2(117) .EQ. 13) .AND. CC2) 
     *   ADCC2 = .TRUE.
 
      If (CCSD) Then

         Write(*,"(a)") " Calculation of CCSD Hbar"

         Call Inilam(Iuhf)
         Call Stllst(Work,Maxcor,Iuhf)
        Write(6,"(a)") "Hbar content at the entry"
        call checkhbar(Work,Macor,Iuhf)

         Call Formwl(work,Maxcor,Iuhf)
         Call Wtwtw(Work,Maxcor,Iuhf)
        
         Term1  = .TRUE.
         Term2  = .TRUE.
         Term3  = .TRUE.
         Term4  = .TRUE.
         Term5  = .TRUE.
         Term6  = .TRUE.
         DO_TAU = .TRUE.

         Call formw4(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)

         Call formw5(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)

         Call Hbrijka0(Work,Maxcor,Iuhf)
         Call Fixfbar(Work,Maxcor,Iuhf)
         Call Modaibc(Work,Maxcor,Iuhf,1.0D0)
         
         Write(*,"(a)") " Calculation of CCSD Hbar is completed"

         Write(*,*) 
         Write(*,"(a)") " Statistics of Hbar elements:" 
         Call Checkhbar(Work,Maxcor,Iuhf)

      Elseif  (CC2) then
            
         CCSD = .TRUE.
         Call Inilam(Iuhf)
         Call Stllst(Work,Maxcor,Iuhf)
         Call Formwl(work,Maxcor,Iuhf)
         Call Wtwtw(Work,Maxcor,Iuhf)

         Term1  = .TRUE.
         Term2  = .TRUE.
         Term3  = .FALSE.
         If (Nonhf) Term3 = .True.
         Term4  = .TRUE.
         Term5  = .FALSE.
         Term6  = .TRUE.
         DO_TAU = .TRUE.

         Call formw4(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)
         Term5   = .TRUE.

         Call formw5(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)

         Call Hbrijka0(Work,Maxcor,Iuhf)
         Call Fixfbar(Work,Maxcor,Iuhf)
         Call Modaibc(Work,Maxcor,Iuhf,1.0D0)

         Write(*,"(a)") " Calculation of CC2 Hbar is completed"

         Write(*,*) 
         Write(*,"(a)") " Statistics of Hbar elements:" 
         Write(*,*) 
         Call Checkhbar(Work,Maxcor,Iuhf)

      Elseif (LCCSD) then

         CCSD = .TRUE.
         Call Inilam(Iuhf)
         Call Stllst(Work,Maxcor,Iuhf)
         Call Formwl(work,Maxcor,Iuhf)
         Call Wtwtw(Work,Maxcor,Iuhf)

         Term1  = .TRUE.
         Term2  = .TRUE.
         Term3  = .FALSE.
         Term4  = .TRUE.
         Term5  = .FALSE.
         Term6  = .TRUE.
         DO_TAU = .FALSE.

         Call formw4(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)
         Term5   = .TRUE.

         Call formw5(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)

         Call Hbrijka0(Work,Maxcor,Iuhf)
         Call Modaibc(Work,Maxcor,Iuhf,1.0D0)

         Write(*,"(a)") " Calculation of LCCSD Hbar is completed"

         Write(*,*) 
         Write(*,"(a)") " Statistics of Hbar elements:" 
         Write(*,*) 
         Call Checkhbar(Work,Maxcor,Iuhf)

      Elseif ((MBPT2 .OR. CCD .OR. LCCD) .AND. (.NOT. ADC2) ) then

C Notice that we are taking advantage of the fact that T1=0.0 (ie. We
C go over the contractions with T1). This may not be ideal in the
C computaional efficiency sense.
C 
         CCSD = .TRUE.
         Call Inilam(Iuhf)
         Call Stllst(Work,Maxcor,Iuhf)

         Call Formwl(work,Maxcor,Iuhf)
         Call Wtwtw(Work,Maxcor,Iuhf)

         Term1  = .TRUE.
         Term2  = .TRUE.
         Term3  = .TRUE.
         Term4  = .TRUE.
         Term5  = .TRUE.
         Term6  = .TRUE.
         DO_TAU = .TRUE.

         Call formw4(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)
         Term5   = .TRUE.

         Call formw5(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)

         Call Hbrijka0(Work,Maxcor,Iuhf)
         Call Fixfbar(Work,Maxcor,Iuhf)
         Call Modaibc(Work,Maxcor,Iuhf,1.0D0)
         
         Write(*,"(a)") " Statistics of Hbar elements:" 
         Write(*,*) 
         Call Checkhbar(Work,Maxcor,Iuhf)
         If (CCD) then
            Write(6,*)
            Write(*,"(a)") " Calculation of CCD Hbar is completed"
         Elseif (MBPT2) then
            Write(6,*)
            Write(*,"(a)")" Calculation of MBPT(2) Hbar is completed"
         Elseif (LCCD) then
            Write(6,*)
            Write(*,"(a)")" Calculation of LCCD Hbar is completed"
         Endif 

      Elseif (MBPT2 .AND. ADC2) then

         CCSD = .TRUE. 
         Call Inilam(Iuhf)
         Call Stllst(Work,Maxcor,Iuhf)
         Call Formwl(work,Maxcor,Iuhf)
         Call Wtwtw(Work,Maxcor,Iuhf)

         Term1  = .TRUE.
         Term2  = .TRUE.
         Term3  = .TRUE.
         Term4  = .TRUE.
         Term5  = .TRUE.
         Term6  = .TRUE.
         DO_TAU = .TRUE.
   
         IRREPX = 1
         If (IUHF .NE. 0) Then
         IAAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
         IBBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         IABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         ELse
         IABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         Endif 

         If (IUHF .NE. 0) Then
            Call Getall(Work, IAAAA_LENGTH_IJAB, IRREPX, 44) 
            Call Dzero(Work,  IAAAA_LENGTH_IJAB)
            Call Putall(Work, IAAAA_LENGTH_IJAB, IRREPX, 44) 

            Call Getall(Work, IBBBB_LENGTH_IJAB, IRREPX, 45) 
            Call Dzero(Work,  IBBBB_LENGTH_IJAB)
            Call Putall(Work, IBBBB_LENGTH_IJAB, IRREPX, 45) 

            Call Getall(Work, IABAB_LENGTH_IJAB, IRREPX, 46) 
            Call Dzero(Work,  IABAB_LENGTH_IJAB)
            Call Putall(Work, IABAB_LENGTH_IJAB, IRREPX, 46) 
         Else
            Call Getall(Work, IABAB_LENGTH_IJAB, IRREPX, 46) 
            Call Dzero(Work,  IABAB_LENGTH_IJAB)
            Call Putall(Work, IABAB_LENGTH_IJAB, IRREPX, 46) 
         Endif 
 
         Call Rnabij(Work, Maxcor, Iuhf, "T")

         Call formw4(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)

         Term5   = .TRUE.

         Call formw5(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)

         Call Hbrijka0(Work,Maxcor,Iuhf)
         Call Fixfbar(Work,Maxcor,Iuhf)
         Call Modaibc(Work,Maxcor,Iuhf,1.0D0)

         Call symmetrize_hbar(Work,Maxcor,Iuhf)

         Write(*,"(a)") " Calculation of ADC(2) Hbar is completed"

         Write(*,*) 
         Write(*,"(a)") " Statistics of Hbar elements:" 
         Write(*,*) 
         Call Checkhbar(Work,Maxcor,Iuhf)
      Elseif (DCCSD) then

         Write(6,*) 
         Write(6,"(a)") " Calculation of pCCSD Hbar"
         Write(6,*) 
         CCSD = .TRUE.
         Call Inilam(Iuhf)
         Call Stllst(Work,Maxcor,Iuhf)
C
C Note that the W(MNIJ) intm. come in as Hbar counterpart. The important 
C term is the  last one (Hbar(MNIJ)) calls for +1/4 Tau(ij,ef)<kl||ef>). 
C Similarly, instead of the Hbar(abcd), we have just the <ab||cd> integral. 
C 
C Construct the Hbar(IMBEJ)
C
         Call dcc_hbar_mbej(Work,Maxcor,IUhf)

         Term1  = .TRUE.
         Term2  = .TRUE.
         Term3  = .TRUE.
         Term4  = .TRUE.
         Term5  = .TRUE.
         Term6  = .TRUE.
         DO_TAU = .TRUE.
C
C Construct the Hbar(IAJK)
C
         Call formw4(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)
C
C Construct the Hbar(ABCI)
C
         Call formw5(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)
C
C Construct the Hbar(IJKA)
C 
         Call Hbrijka0(Work,Maxcor,Iuhf)
C
C There are two verison of Hbar(AE) and Hbar(MI) and The Hbar(ME) is the 
c standard CCSD F(ME). The lists 91(1,2) and 92(1,2) contains the Hbar
C (MI) and (AE) as it is for the standard CCSD calculation. The list 
C 91(10,11) and 92(10,11) are scalled according to the pCC parameters. 
C Both copies are needed in EOM calculations. 
C
         Call dcc_hbar_fae(Work,Maxcor,Iuhf)
         Call dcc_hbar_fmi(Work,Maxcor,Iuhf)
C 
C Construct the Hbar(AIBC)
C 
         Call Modaibc(Work,Maxcor,Iuhf,1.0D0)
C
C Add the (+-)12F(EM)*T(E,I) (or F(EM)*TM(A,M)) to the F(MI) and
C F(EA) intermediates. We need to it for both copies (one is scalled
C and the other is not.
c

         Call Fixfbar(Work,Maxcor,Iuhf)
         Call Dcc_fixfbar(Work,Maxcor,Iuhf)

         Call checkhbar(Work,Length,Iuhf)
         Write(*,"(a)") " Calculation of pCCSD Hbar is completed"

      Elseif (DCCD) then

         CCSD = .TRUE.
         Call Inilam(Iuhf)
         Call Stllst(Work,Maxcor,Iuhf)

         Call dcc_hbar_mbej(Work,Maxcor,IUhf)

         Term1  = .TRUE.
         Term2  = .TRUE.
         Term3  = .TRUE.
         Term4  = .TRUE.
         Term5  = .TRUE.
         Term6  = .TRUE.
         DO_TAU = .TRUE.
C
C Construct the Hbar(IAJK)
C
         Call formw4(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)
C
C Construct the Hbar(ABCI)
C
         Call formw5(Work,Maxcor,Iuhf,Term1,Term2,Term3,Term4,Term5,
     &              Term6,Do_Tau,.True.)
C
C Construct the Hbar(IJKA)
C
         Call Hbrijka0(Work,Maxcor,Iuhf)
C
C There are two verison of Hbar(AE) and Hbar(MI) and The Hbar(ME) is the 
c standard CCSD F(ME). The lists 91(1,2) and 92(1,2) contains the Hbar
C (MI) and (AE) as it is for the standard CCSD calculation. The list 
C 91(10,11) and 92(10,11) are scalled according to the pCC parameters. 
C Both copies are needed in EOM calculations. 

         Call dcc_hbar_fae(Work,Maxcor,Iuhf)
         Call dcc_hbar_fmi(Work,Maxcor,Iuhf)
C 
         Call Modaibc(Work,Maxcor,Iuhf,1.0D0)
C
C Add the (+-)12F(EM)*T(E,I) (or F(EM)*TM(A,M)) to the F(MI) and
C F(EA) intermediates. We need to it for both copies (one is scalled
C and the other is not.
c
         Call Fixfbar(Work,Maxcor,Iuhf)
         Call Dcc_fixfbar(Work,Maxcor,Iuhf)

         Call checkhbar(Work,Length,Iuhf)

         Write(*,"(a)") " Calculation of pCCD Hbar is completed"

      Endif 

      Return
      End


  
  
