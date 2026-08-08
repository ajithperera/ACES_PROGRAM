










      SUBROUTINE FORMX(ICORE,MAXCOR,IUHF)
C
C   THIS ROUTINE CALCULATES THE X(IJ,AB) CONTRIBUTION REQUIRED IN 
C   METHODS WHICH INVOLVE QUADRUPLE EXCITATION TO COMPUTE GRADIENTS
C
C   MBPT(4) :   X(IJ,AB) =   1/4  T(MN,EF) T(IJ,EF) T(MN,AB)
C
C                          - 1/4 P(AB) T(MN,EF) T(MN,BE) T(IJ,AF)
C
C                          - 1/4 P(IJ) T(MN,EF) T(IM,EF) T(JN,AB)
C
C                          + 1/8 P(IJ)P(AB) T(MN,EF) T(JN,BF) T(IN,AF)
C
C  CCD  :      X(IJ,AB) =   1/4  L(MN,EF) T(IJ,EF) T(MN,AB)
C    
C                         - 1/4 P(AB) L(MN,EF) T(MN,BE) T(IJ,AF)
C
C                         - 1/4 P(IJ) L(MN,EF) T(IM,EF) T(JN,AB)
C 
C                         + 1/8 P(IJ) P(AB) L(MN,EF) T(JN,BF) T(MN,EF)
C
C  QCISD :     THE SAME AS CCD
C
C  CCSD :      ADDITIONAL TAU-CONTRIBUTIONS AS GIVEN BELOW
C
CEND
C
C CODED AUGUST/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      INTEGER POP,VRT,DIRPRD
      INTEGER AAAA_LENGTH_IJAB,BBBB_LENGTH_IJAB,ABAB_LENGTH_IJAB
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

C
      DATA HALF/0.5D0/
C
      IF(IFLAGS(3).EQ.2) THEN
       CALL SETLST2(ICORE,MAXCOR,IUHF)
      ENDIF

CSSS#ifdef 1
C We need to construct the P(D)CCSD V(mn,ij) intermediate.
      If (Ispar) Then
         CALL PDCC_FORMV1(ICORE,MAXCOR,IUHF)
      Endif 
CSSS#endif 
C
C CALCULATE THE INTERMEDIATE  H = - T(IM,AE) T(JN,BF)
C 
C For PCC and DCC the work done in the above block has been 
C moved to gamdrv.

      IF (.NOT. (ISPAR .AND. COULOMB)) CALL FORMH4(ICORE,MAXCOR,IUHF)
C
C CONTRACT THE INTERMDIATES TO GIVE THE X-CONTRIBUTION
C
      CALL V1INX2(ICORE,MAXCOR,IUHF)

      Write(6,*) "After V1 in X2"
C Gamma(IJ,AB)
      Write(*,*)
      Irrepx=1
      IF (IUHF .Gt.0) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,114),ISYTYP(2,114))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,115),ISYTYP(2,115))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,116),ISYTYP(2,116))
         Call Getall(ICORE, AAAA_LENGTH_IJAB, IRREPX, 114)
         Call Checksum("AAAA_IJAB", ICORE, AAAA_LENGTH_IJAB)
         Call Getall(ICORE, BBBB_LENGTH_IJAB, IRREPX, 115)
         Call Checksum("BBBB_IJAB", ICORE, BBBB_LENGTH_IJAB)
         Call Getall(ICORE, ABAB_LENGTH_MNIJ, IRREPX, 116)
         Call Checksum("ABAB_IJAB", ICORE, ABAB_LENGTH_IJAB)
      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,116),ISYTYP(2,116))
         Call Getall(ICORE, ABAB_LENGTH_IJAB, IRREPX, 116)
         Call Checksum("ABAB_IJAB", ICORE, ABAB_LENGTH_IJAB)
      Endif
C
C For PCC, the work done in the above block has been moved to gamdrv.
C Also, notes the changes to h4inx2. For PCC and DCC, the contributions
C are broken down into pieces. 

      IF (.NOT. (ISPAR .AND. COULOMB)) CALL H4INX2(ICORE,MAXCOR,IUHF)
C
      If (Ispar) Then
C
C Note that in above PDCC_FORMV1, the V(MNIJ) is formed with
C paralpha scalling for T2 contributions. Since the above
c contraction is (T2 + T1T1) (gamma*T2+T1T1)*L2, we need add a term
C (1-gamma)*T1*T1*T2*L2. Lets, add it here.

        CALL CONSTRUCT_T2L2V1(ICORE,MAXCOR,IUHF)
        CALL CONSTRUCT_T1T1T2L2(ICORE,MAXCOR,IUHF)

C Substract the extra piece that was done in H4INX2. This is
C the scalling of the W(mb,ej) intermediate. The last argument
C is to false to indicate that there is no Tau formation.
C
         If (Ispar .AND. .NOT. Coulomb) Then
             IBOT=1
             IF(IUHF.EQ.0) IBOT=3

             DO ISPIN=IBOT,3
                CALL H4X2_4PDCC(ICORE,MAXCOR,ISPIN,IUHF,.FALSE.)
             ENDDO 
         Endif 
C
      Write(6,*) "DCC: After V1 and H4 in X2"
      CALL Checkgamma1(Icore,Maxcor,Iuhf)
C     
C for pCC the  G(mi) and G(ae) intermediates are scaled 
C See post_dccl_mods.F in lambda code. 01/2016, Ajith Perera
C
         cALL G1INX2(ICORE,MAXCOR,IUHF,8)
      Write(6,*) "DCC: After scalled G1 in X2"
      CALL Checkgamma1(Icore,Maxcor,Iuhf)
         CALL G2INX2(ICORE,MAXCOR,IUHF,8)

      Write(6,*) "DCC: After scalled G1 and G2 in X2"
      CALL Checkgamma1(Icore,Maxcor,Iuhf)
C
C Note that the G(mi) and G(ai) stored in 191,192 (9,10) locations
C are scalled with the parameters Gae_scale and Gmi_scale. Since
C the above contraction is (T2 +  T1T1)(C * T2*L2), it introduce
C and errorneous contribution C*T1T1*T2*L2 instead of T1T1*T2*L2
C In order to correct this lets add (1-C)*T1T1*T2*L2 contribution.
C (C=Fae_scale,Fmi_scale)
C
         If (.NOT. Coulomb) Then
            CALL CONSTRUCT_T1T1G1(ICORE,MAXCOR,IUHF,0,1.0D0)
            CALL CONSTRUCT_T1T1G1(ICORE,MAXCOR,IUHF,0,-Gae_scale)
         Endif 
      Write(6,*) "DCC: After correcting for G1 scales in X2"
      CALL Checkgamma1(Icore,Maxcor,Iuhf)
         CALL CONSTRUCT_T1T1G2(ICORE,MAXCOR,IUHF,0,1.0D0)
         CALL CONSTRUCT_T1T1G2(ICORE,MAXCOR,IUHF,0,-Gmi_scale)

      Write(6,*) "DCC: After correcting for G1 and G2i scales in X2"
      CALL Checkgamma1(Icore,Maxcor,Iuhf)
      Else
CSSS#else
         CALL G1INX2(ICORE,MAXCOR,IUHF,0)
         CALL G2INX2(ICORE,MAXCOR,IUHF,0)
      Endif 
CSSS#endif 

      Write(6,*)
      Write(6,*) "After G1 and G2 in X2"
      CALL Checkgamma1(Icore,Maxcor,Iuhf)

CSSS#ifdef 1
C We need to reconstruct the V(mn,ij) CCSD intermediate. 
      If (Ispar) Then
         CALL FORMV1(ICORE,MAXCOR,IUHF)
      Endif 
CSSS#endif 
C
C ALL DONE 
C
      RETURN
      END
