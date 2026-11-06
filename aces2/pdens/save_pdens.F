      SUBROUTINE SAVE_PDENS(DOOA,DOOB,DVVA,DVVB,DVOA,DVOB,DOVA,
     &                      DOVB,IRREPX,IUHF)
C
C THIS SUBROUTINE SAVES THE 1ST ORDER PERTURB DENSITY MATRIX ON THE
C GAMLAM FILE. THIS IS ONLY REQUIRED FOR SECOND DERIVATIVE CALCULATIONS.
C
C THE LISTS WRITTEN ARE:
C
C  DOOA    1,495
C  DOOB    2,495    UHF AND ROHF ONLY
C  DVVA    3,495
C  DVVB    4,495    UHF AND ROHF ONLY
C  DVOA    5,495  
C  DVOB    6,495    UHF AND ROHF ONLY
C  DOVA    7,495
C  DOVB    8,495    UHF AND ROHF ONLY
C

      INTEGER POP,VRT
      LOGICAL ANTI_SYMM
      DIMENSION DOOA(1),DOOB(1),DVVA(1),DVVB(1),DVOA(1),
     &          DVOB(1),DOVA(1),DOVB(1)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      LOGICAL MADE
      SAVE MADE
      DATA MADE /.FALSE./
C
      IENTER = 0
      IOFF = 0
C
C THE DENSITY HAS THE SYMMETRY OF THE PERTURBATION (IRREPX), SO THE
C LISTS ARE MADE BIG ENOUGH FOR ANY IRREP AND THEN RESIZED TO IRREPX.
C
      DO 100 ISPIN=1,IUHF+1
C
       MAXOO = 0
       MAXVV = 0
       MAXVO = 0
       DO IRREP = 1, NIRREP
          MAXOO = MAX(MAXOO, IRPDPD(IRREP,20+ISPIN))
          MAXVV = MAX(MAXVV, IRPDPD(IRREP,18+ISPIN))
          MAXVO = MAX(MAXVO, IRPDPD(IRREP, 8+ISPIN))
       ENDDO
       IF (.NOT. MADE) THEN
          CALL UPDMOI(3,MAXOO,ISPIN,495,IENTER,IOFF)
          CALL UPDMOI(3,MAXVV,2+ISPIN,495,IENTER,IOFF)
          CALL UPDMOI(3,MAXVO,4+ISPIN,495,IENTER,IOFF)
          CALL UPDMOI(3,MAXVO,6+ISPIN,495,IENTER,IOFF)
       ENDIF
       CALL ACES_LIST_RESIZE(ISPIN,  495,IRPDPD(IRREPX,20+ISPIN))
       CALL ACES_LIST_RESIZE(2+ISPIN,495,IRPDPD(IRREPX,18+ISPIN))
       CALL ACES_LIST_RESIZE(4+ISPIN,495,IRPDPD(IRREPX, 8+ISPIN))
       CALL ACES_LIST_RESIZE(6+ISPIN,495,IRPDPD(IRREPX, 8+ISPIN))
C
100   CONTINUE
      MADE = .TRUE.
C
#if _DEBUG_LVLM
      Print*, "-------The A and B Density Mats in SAVE_PDENS ------"
      Print*, "occ-occ A and B"
      call output(dooa, 1, nf1(1), 1, 1, nf1(1), 1, 1)
      Write(6,*)
      call output(doob, 1, nf1(2), 1, 1, nf1(2), 1, 1)
      Print*, "vrt-vrt A and B"
      call output(dvva, 1, nf2(1), 1, 1, nf2(1), 1, 1)
      Write(6,*)
      call output(dvvb, 1, nf2(2), 1, 1, nf2(2), 1, 1)
      Write(6,*)
      Print*, "occ-vrt A and B"
      call output(dova, 1, nt(1), 1, 1, nt(1), 1, 1)
      Write(6,*)
      call output(dovb, 1, nt(2), 1, 1, nt(2), 1, 1)
      Write(6,*)
      Print*, "vrt-occ A and B"
      call output(dvoa, 1, nt(1), 1, 1, nt(1), 1, 1)
      Write(6,*) 
      call output(dvob, 1, nt(2), 1, 1, nt(2), 1, 1)

#endif
      CALL PUTLST(DOOA,1,1,1,1,495)
      CALL PUTLST(DVVA,1,1,1,3,495)
      CALL PUTLST(DVOA,1,1,1,5,495)
      CALL PUTLST(DOVA,1,1,1,7,495)
C
      IF (IUHF.NE.0) THEN
C
       CALL PUTLST(DOOB,1,1,1,2,495)
       CALL PUTLST(DVVB,1,1,1,4,495)
       CALL PUTLST(DVOB,1,1,1,6,495)
       CALL PUTLST(DOVB,1,1,1,8,495)
C
      ENDIF
C
      RETURN
      END
