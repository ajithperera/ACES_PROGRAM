










C
      SUBROUTINE Scrnc_form_mubar_s(MIJ,MAB,MAI,MIJBAR,MABBAR,MAIBAR,
     &                             MIABAR,SCR,MAXCOR,IUHF,IRREPX,
     &                             ISIDE,MBAR00,IEOMPROP)
C
C THIS ROUTINE CALCULATES ALL ONE-PARTICLE PARTS OF THE TRANSFORMED
C ONE-PARTICLE OPERATOR:
C
C                _
C                M = exp(-T) M exp(T)
C
C TO WIT:
C    _ 
C    M (ij) = M(ij) + T(je) M(ei) 
C    _
C    M (ab) = M(ab) - T(ma) M(mb)
C    _                                       _             _
C    M (ai) = M(ai) + T(im,ae) M(me) + T(ie) M(ae) - T(ma) M(mi)
C    _                                        
C    M (ia) = M(ia) (If antisymmetric M(ia) = -M(ai))
C
C NOTE THAT THE M(ae) ELEMENTS USED IN THE THIRD EQUATION CORRESPOND
C TO THE *TRANSFORMED* DIPOLE OPERATOR, WHICH IS HOW THE T1*T1 TERMS
C ARE INCLUDED. 
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL PRINT
      LOGICAL Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      DOUBLE PRECISION ONE,ONEM,ZILCH,SNRM2,TWO,SCR,SDOT,
     &                 SIGN1, SIGN2, HALF, HALFM
      DOUBLE PRECISION MIJ(*),MAB(*),MAI(*),MEXPECT,MBAR00
      DOUBLE PRECISION MIJBAR(*),MABBAR(*),MAIBAR(*),MIABAR(*)
      DIMENSION SCR(MAXCOR)
      COMMON /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +                Drccd
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ONEM,ZILCH,TWO /1.0D0,-1.0D0,0.0D0,2.0D0/
      DATA HALF,HALFM /0.5D0,-0.5D0/
C
      LISTT1=90+100*(ISIDE-1)
C
      SIGN1 = ONE
      SIGN2 = ONE

C DO OCCUPIED-OCCUPIED PART
C      
C Iside = 1      
C    _ 
C    M (ij) = M(ij) + M(ie) * T(ej) = M(ij) + M(ei) * T(ej)
C    _
C    M (ij) = M(ij) - M(ei) * T(ej) if the perturbation is antisymmetric
C
C Iside = 2
C    _                _                       _
C    L (ij) = M(ji) + M(ei) * L(je) = M(ij) + M(ei) * L(je)
C    _                  _
C    L (ij) = - M(ij) + M(ei) * L(ej) if the perturbation is antisymmetric
C           _                                                         _
C Note that L(ji) which is what actually being evaluated is stored in L(ij).
C           _
C Moreover, L(ij) has no use in the rest of the program. 

      DO 10 ISPIN=1,1+IUHF
       LENOO=IRPDPD(IRREPX,20+ISPIN)
       LENVO=IRPDPD(IRREPX,8+ISPIN)
       I000=1
       CALL GETLST(SCR(I000),1,1,1,ISPIN,LISTT1)
       IOFFT0=I000
       IOFFOO0=1+(ISPIN-1)*IRPDPD(IRREPX,21)
       IOFFVO0=1+(ISPIN-1)*IRPDPD(IRREPX,9)
       DO 30 IRREPJ=1,NIRREP
        IRREPI=DIRPRD(IRREPJ,IRREPX)
        IRREPE=IRREPJ
        NUMJ=POP(IRREPJ,ISPIN)
        NUMI=POP(IRREPI,ISPIN)
        NUME=VRT(IRREPE,ISPIN)
        NROW=NUMI
        NCOL=NUMJ
        NSUM=NUME
        IOFFT=IOFFT0+(ISYMOFF(IRREPJ,1,8+ISPIN)-1)
        IOFFOO=IOFFOO0+(ISYMOFF(IRREPJ,IRREPX,20+ISPIN)-1)
        IOFFVO=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE*SIGN1,MAI(IOFFVO),
     &              NSUM,SCR(IOFFT),NSUM,ONE*SIGN2,MIJBAR(IOFFOO),NROW)
30     CONTINUE
10    CONTINUE
C
C DO VIRTUAL-VIRTUAL PART
C
C Iside = 1      
C    _
C    M (ab) = M(ba) - M(ma) * T(bm) = M(ab) - M(am) * T(bm) 
C    _
C    M (ab) = - M(ab) - (-) M(am) * T(bm) if the perturbation is antisymmetric
C 
C This is little tricky when the perturbation is not symmetric. The actual
C                           _                          _
C contraction evaluates the M(ba), but it is stored in M(ab) which is consistent
C with rest of the ACESII ordering of (virt,virt) quantities. But bare M is 
C ordered as M(ab) and note the negative sign in front of M(ab) when the
C perturbation is not symmetric.
C
C Iside = 2
C    _                _
C    L (ab) = M(ab) - M(am) * L(mb) 
C    _                _
C    L (ab) = M(ab) - M(am) * L(mb) if the perturbation is antisymmetric

      DO 110 ISPIN=1,1+IUHF
       LENVV=IRPDPD(IRREPX,18+ISPIN)
       LENVO=IRPDPD(IRREPX,8+ISPIN)
       I000=1
       I010=I000+NT(ISPIN)
       CALL GETLST(SCR(I000),1,1,1,ISPIN,LISTT1)
       IOFFT0=I000
       IOFFVV0=1+(ISPIN-1)*IRPDPD(IRREPX,19)
       IOFFVO0=1+(ISPIN-1)*IRPDPD(IRREPX,9)
       DO 130 IRREPB=1,NIRREP
        IRREPA=DIRPRD(IRREPB,IRREPX)
        IRREPM=IRREPB
        NUMA=VRT(IRREPA,ISPIN)
        NUMB=VRT(IRREPB,ISPIN)
        NUMM=POP(IRREPM,ISPIN)
        NROW=NUMA
        NCOL=NUMB
        NSUM=NUMM
        IOFFT=IOFFT0+(ISYMOFF(IRREPM,1,8+ISPIN)-1)
        IOFFVV=IOFFVV0+(ISYMOFF(IRREPB,IRREPX,18+ISPIN)-1)
        IOFFVO=IOFFVO0+(ISYMOFF(IRREPM,IRREPX,8+ISPIN)-1)
        CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM*SIGN1,MAI(IOFFVO),
     &             NROW,SCR(IOFFT),NCOL,ONE*SIGN1,MABBAR(IOFFVV),
     &             NROW)
130    CONTINUE
110   CONTINUE
C
C DO VIRTUAL-OCCUPIED PART
C
C Iside = 1
C    _                                       _             _
C    M (ai) = M(ai) + T(ae,im) M(me) + T(ei) M(ae) - T(am) M(mi)
C
C Iside = 2       
C    _                                       _             _
C    L (ia) = M(ia) + L(im,ae) M(em) + L(ie) M(ea) - L(ma) M(im)
C 
      DO 210 ISPIN=1,1+IUHF
       LENVV=IRPDPD(IRREPX,18+ISPIN)
       LENVO=IRPDPD(IRREPX,8+ISPIN)
       LENOO=IRPDPD(IRREPX,20+ISPIN)
       LENT1=NT(ISPIN)
       LISTT2AB=100*(ISIDE-1)+38-ISPIN
       LISTT2AA=100*(ISIDE-1)+33+ISPIN
       NT2DSZAB=IRPDPD(IRREPX,ISYTYP(1,LISTT2AB))
       NT2DISAB=IRPDPD(IRREPX,ISYTYP(2,LISTT2AB))
       NT2DSZAA=IRPDPD(IRREPX,ISYTYP(1,LISTT2AA))
       NT2DISAA=IRPDPD(IRREPX,ISYTYP(2,LISTT2AA))
       I000=1
       I0GV=I000+LENT1
       I0GO=I0GV+NFEA(ISPIN)
       I010=I0GO+NFMI(ISPIN)

C Bug fix 2026-08-14: was NT2DSZAB*NT2DSZAB (squared -- copy-paste typo
C off the AA term above, which correctly uses SZ*DIS). Numerically inert
C for closed-shell-like UHF references where alpha/beta occ/virt counts
C match (NT2DSZAB happens to equal NT2DISAB there), which is why this
C survived undetected -- but wrong in general, and a real corruption risk
C (I020 landing too close, overlapping still-needed T2AB data) for any
C case where the AB block truly is non-square.
       I020=I010+MAX(NT2DSZAA*NT2DISAA,NT2DSZAB*NT2DISAB)
       IF(IUHF.EQ.0)I030=I020+NT2DSZAA*NT2DISAA
       CALL GETLST(SCR(I000),1,1,1,ISPIN,LISTT1)
       IF(ISIDE.EQ.2)THEN
        CALL GETLST(SCR(I0GO),1,1,1,ISPIN,191)
        CALL GETLST(SCR(I0GV),1,1,1,ISPIN,192)
       ENDIF
       IOFFT0=I000
       IOFFVV0=1+(ISPIN-1)*IRPDPD(IRREPX,19)
       IOFFOO0=1+(ISPIN-1)*IRPDPD(IRREPX,21)
       IOFFVO0=1+(ISPIN-1)*IRPDPD(IRREPX,9)
       DO 230 IRREPI=1,NIRREP
        IRREPA=DIRPRD(IRREPI,IRREPX)
        NUMA=VRT(IRREPA,ISPIN)
        NUMI=POP(IRREPI,ISPIN)
        NROW=NUMA
        NCOL=NUMI
        IOFFVO=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
C
C Iside = 1
C  
C    RHS(ai) = M(ae)  * T(ei) + M(ae) [Iside = 1]
C
C Iside = 2
C
C     LHS(ai) = M(ea) * L(ie) +  M(ia) = M(ea) * L(ie) + M(ai)
C     LHS(ai) = M(ea) * L(ie) - M(ai) if the perturbation is antisymmetric
C
C Refer to the earlier comments. 
C
        IRREPE=IRREPI
        NUME=VRT(IRREPE,ISPIN)
        NSUM=NUME
        IOFFT=IOFFT0+(ISYMOFF(IRREPI,1,8+ISPIN)-1)
        IOFFVV=IOFFVV0+(ISYMOFF(IRREPE,IRREPX,18+ISPIN)-1)
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,MAB(IOFFVV),
     &             NROW,SCR(IOFFT),NSUM,ONE*SIGN2,MAIBAR(IOFFVO),NROW)
C 
C Iside = 1
C                        _  
C    RHS(ai) = - T(am) * M(mi) + RHS(ai)
C          
C Iside = 2
C                        _                      
C    LHS(ai) = - L(ma) * M(im) + LHS(ai)
C
C  Refer to the earlier comments
C
        IRREPM=IRREPA
        NUMM=POP(IRREPM,ISPIN)
        NSUM=NUMM
        IOFFT=IOFFT0+(ISYMOFF(IRREPM,1,8+ISPIN)-1)
        IF(ISIDE.EQ.1)THEN
         IOFFOO=IOFFOO0+(ISYMOFF(IRREPI,IRREPX,20+ISPIN)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,SCR(IOFFT),
     &              NROW,MIJBAR(IOFFOO),NSUM,ONE,MAIBAR(IOFFVO),NROW)
        ELSE
         IOFFOO=IOFFOO0+(ISYMOFF(IRREPM,IRREPX,20+ISPIN)-1)
         CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,SCR(IOFFT),
     &              NROW,MIJ(IOFFOO),NCOL,ONE,MAIBAR(IOFFVO),NROW)
        ENDIF
C
        IF(ISIDE.EQ.2)THEN
C    
C  LHS(ai) = G(ea) * M(ie) + LHS(ai) [Iside = 2]
C  LHS(ai) = M(ma) * G(im) + LHS(ai) [Iside = 2]
C
         IRREPE=DIRPRD(IRREPI,IRREPX)
         NUME=VRT(IRREPE,ISPIN)
         IOFFG=I0GV+(ISYMOFF(IRREPA,1,18+ISPIN)-1)
         IOFFM=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
         IOFFTAR=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
         NSUM=NUME
         CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,SCR(IOFFG),
     &              NSUM,MIABAR(IOFFVO),NSUM,ONE,MAIBAR(IOFFTAR),NROW)
C
C  LHS(ai) = M(ma) * G(im) + LHS(ai) [Iside = 2]
C
         IRREPM=DIRPRD(IRREPA,IRREPX)
         NUMM=POP(IRREPM,ISPIN)
         IOFFG=I0GO+(ISYMOFF(IRREPI,1,20+ISPIN)-1)
         IOFFM=IOFFVO0+(ISYMOFF(IRREPM,IRREPX,8+ISPIN)-1)
         IOFFTAR=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
         NSUM=NUMM
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,MIABAR(IOFFM),NROW,
     &              SCR(IOFFG),NSUM,ONE,MAIBAR(IOFFTAR),NROW)
        ENDIF 
C
230    CONTINUE
C 
C [Iside = 1]
C
C  RHS(ai) = T(ai,em)*M(me) + RHS(ai) = T(ai,em)*M(em) + RHS(ai)
C
C  RHS(ai) = - T(ai,em)*M(em) + RHS(ai) If the perturbation is antisymmetric
C 
C [Iside = 2]
C                     _
C  LHS(ai) = L(ai,em)*M(em) + LHS(ai)
C
C If ISIDE = 1 or (ISIDE = 2 and IEOMPROP not equal to 2) include the
C following contributions. 

      IF (.NOT. (ISIDE .EQ. 2 .AND. IEOMPROP .EQ. 2)) THEN
C
          IF(ISPIN.EQ.1)THEN
             IOFFVO1=1
             IOFFVO2=1+IUHF*IRPDPD(IRREPX,9)
          ELSE
             IOFFVO1=1+IUHF*IRPDPD(IRREPX,9)
             IOFFVO2=1
          ENDIF
C
          CALL GETLST(SCR(I010),1,NT2DISAB,1,IRREPX,LISTT2AB)
C 2026-08-22: raw, pre-combination T(list37-family) checksum -- settles
C whether LISTT2AB (Iside-dependent, =37 R-side/=137 L-side for RCCD
C ISPIN=1) actually equals the converged singlet channel amplitude
C T_S, independent of the 2*T37-T34 combination done below. See
C SCRNC_DEBUG_LOG.md, 2026-08-22 (Lead B empirical checksum test).
      Write(6,"(a,I2,a,I5,a,20F15.10)")
     +     " LEADB_RAW37 Iside=",Iside," LISTT2AB=",LISTT2AB," val=",
     +     (SCR(I010-1+I),I=1,MIN(NT2DSZAB*NT2DISAB,20))
          IF(IUHF.EQ.0)THEN
C 2026-08-19: RCCD's own T2(Ij,Ab) (list 37, read above as LISTT2AB)
C is never antisymmetrized (rcc_dwmbej.F's own comment) -- unlike
C CCSD's, so list 39 (the generic RNABIJ 'Ab<->bA'-type resort of list
C 46, valid for CCSD's antisymmetrized amplitudes) turns out IDENTICAL
C to list 37 for RCCD, element for element (confirmed empirically,
C H2/STO-3G and H2/DZ both). That silently collapses this routine's
C "2*T37-T39" CCSD-style dressing combination to just "T37" for RCCD,
C losing the intended second (exchange-type) contribution entirely --
C the real source of the still-open RHF/UHF response mismatch (see
C SCRNC_DEBUG_LOG.md). rcc_dwmbej.F's OWN Hbar construction already
C solves this correctly using list 34 (not 39) as list 37's exchange-
C type partner -- reuse that same list here for RCCD specifically.
C CCSD (and DRCCD, deferred separately, not fixed here) keep list 39.
             LISTT2B = 100*(ISIDE-1) + 39
             IF (RCCD) LISTT2B = 100*(ISIDE-1) + 34
             CALL GETLST(SCR(I020),1,NT2DISAA,1,IRREPX,LISTT2B)
C 2026-08-22: raw list34 (RCCD's "exchange partner", flagged in
C SCRNC_DEBUG_LOG.md 2026-08-20 "continued 3" as actually being ring-CC's
C TRIPLET channel's own primary amplitude, REFLST=44/SPCAS2='AAAA') --
C print BEFORE the SAXPY combination below so it can be compared
C directly against T_T=+0.163485120840 (H2/STO-3G RHF RCCD, list 61).
      Write(6,"(a,I2,a,I5,a,20F15.10)")
     +     " LEADB_RAW34 Iside=",Iside," LISTT2B =",LISTT2B," val=",
     +     (SCR(I020-1+I),I=1,MIN(NT2DSZAA*NT2DISAA,20))
C 2026-08-22 FIX (see SCRNC_DEBUG_LOG.md "Lead B" test): the empirical
C checksum test above confirms list34/list134 EXACTLY equals the
C converged triplet-channel amplitude T_T (RCC_RHF_CMPENEG list 61) and
C list37/list137 EXACTLY equals the converged singlet-channel amplitude
C T_S (list 63) -- both to ~8 printed digits. Finding 1's proven
C spin-adaptation identity T(ABAB) = (T_S - T_T)/2 (see log -- this is
C the SAME quantity UHF's raw, un-dressed list37 IS) therefore predicts
C the correct RHF reconstruction of that spin-orbital ABAB amplitude is
C (T37-T34)/2, NOT the CCSD-style 2*T37-T34 previously shipped here
C (empirically tuned, not derived, in 498542f8 -- see 2026-08-20
C "continued 3", which already cast doubt on its theoretical grounding).
C VERIFIED: H2/STO-3G RCCD response-only ratio (UHF/RHF) improved
C 5.18x->1.88x; H2O/STO-3G RCCD improved on ALL 49 diagonal JJ elements,
C no regressions. Guarded to RCCD only; CCSD's own list-39
C antisymmetrized-amplitude formula (2*T37-T39) is untouched; DRCCD
C (separate flag) also untouched, reconfirmed bit-identical.
             IF (RCCD) THEN
                CALL SSCAL(NT2DSZAB*NT2DISAB,HALF,SCR(I010),1)
                CALL SAXPY(NT2DSZAB*NT2DISAB,HALFM,SCR(I020),1,
     &                     SCR(I010),1)
             ELSE
                CALL SSCAL(NT2DSZAB*NT2DISAB,TWO,SCR(I010),1)
                CALL SAXPY(NT2DSZAB*NT2DISAB,ONEM,SCR(I020),1,
     &                     SCR(I010),1)
             ENDIF
          ENDIF
      call checksum("T2AB(pre-xgemm)              :",SCR(I010),
     +               NT2DSZAB*NT2DISAB)
      call checksum("MAI(IOFFVO2)(pre-xgemm)      :",MAI(IOFFVO2),
     +               NT2DISAB)
      Write(6,"(a,I2,a,I2,a,I2,a,I5,a,I5)")
     +     " RAWTRACE Iside=",Iside," Irrepx=",Irrepx," Ispin=",Ispin,
     +     " NT2DSZAB=",NT2DSZAB," NT2DISAB=",NT2DISAB
      Write(6,"(a,20F15.10)") " RAWTRACE T2AB(1:NT2DSZAB,dis=1)  :",
     +     (SCR(I010-1+I),I=1,MIN(NT2DSZAB,20))
      Write(6,"(a,20F15.10)") " RAWTRACE MAI(IOFFVO2)            :",
     +     (MAI(IOFFVO2-1+I),I=1,MIN(NT2DISAB,20))
          CALL XGEMM('N','N',NT2DSZAB,1,NT2DISAB,ONE*SIGN1,SCR(I010),
     &                NT2DSZAB,MAI(IOFFVO2),NT2DISAB,ONE,
     &                MAIBAR(IOFFVO1),NT2DSZAB)
      call checksum("MAIBAR(post-AB-xgemm)        :",MAIBAR(IOFFVO1),
     +               NT2DSZAB)
      Write(6,"(a,20F15.10)") " RAWTRACE MAIBAR(post-AB-xgemm)   :",
     +     (MAIBAR(IOFFVO1-1+I),I=1,MIN(NT2DSZAB,20))
          IF(IUHF.NE.0)THEN
             CALL GETLST(SCR(I010),1,NT2DISAA,1,IRREPX,LISTT2AA)
      call checksum("T2AA(pre-xgemm)              :",SCR(I010),
     +               NT2DSZAA*NT2DISAA)
      call checksum("MAI(IOFFVO1)(pre-xgemm)      :",MAI(IOFFVO1),
     +               NT2DISAA)
      Write(6,"(a,I2,a,I2,a,I2,a,I5,a,I5)")
     +     " RAWTRACE Iside=",Iside," Irrepx=",Irrepx," Ispin=",Ispin,
     +     " NT2DSZAA=",NT2DSZAA," NT2DISAA=",NT2DISAA
      Write(6,"(a,20F15.10)") " RAWTRACE T2AA(1:NT2DSZAA,dis=1)  :",
     +     (SCR(I010-1+I),I=1,MIN(NT2DSZAA,20))
      Write(6,"(a,20F15.10)") " RAWTRACE MAI(IOFFVO1)            :",
     +     (MAI(IOFFVO1-1+I),I=1,MIN(NT2DISAA,20))
             CALL XGEMM('N','N',NT2DSZAA,1,NT2DISAA,ONEM*SIGN1,
     &                   SCR(I010),NT2DSZAA,MAI(IOFFVO1),NT2DISAA,ONE,
     &                   MAIBAR(IOFFVO1),NT2DSZAA)
      call checksum("MAIBAR(post-AA-xgemm)        :",MAIBAR(IOFFVO1),
     +               NT2DSZAA)
      Write(6,"(a,20F15.10)") " RAWTRACE MAIBAR(post-AA-xgemm)   :",
     +     (MAIBAR(IOFFVO1-1+I),I=1,MIN(NT2DSZAA,20))
          ENDIF
       ENDIF
C
       IF (ISIDE .EQ. 2 .AND. IRREPX .EQ. 1 .AND. IEOMPROP .EQ. 0) THEN
          CALL GETLST(SCR(I000),1,1,1,ISPIN,190)
          CALL SAXPY (NT(ISPIN),MBAR00,SCR(I000),1,MAIBAR(IOFFVO0),1)
       ENDIF
       If (iside .eq.1)  call Checksum("scrnc_form_mubar_r(f)        :", 
     +                                  MAIBAR(Ioffvo0),Lenvo)

       If (iside .eq.2)  call Checksum("scrnc_form_mubar_l(f)        :", 
     +                                 MAIBAR(Ioffvo0),Lenvo)
       IF (ISIDE .EQ. 1) THEN
          CALL PUTLST(MAIBAR(IOFFVO0),1,1,1,ISPIN,480)
       ELSE IF (ISIDE .EQ. 2) THEN
          CALL PUTLST(MAIBAR(IOFFVO0),1,1,1,ISPIN,482)
       ENDIF 
C
210   CONTINUE
C
      RETURN
      END
