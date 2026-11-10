C
      SUBROUTINE R2rmom(VLIST,NOC,CENTR1,EXPS,CONT,LABEL,NPRIM,NCONTR,
     &                  ITRM,ITRN,CTRN,NCONT,BUF,IBUF,AOINT,LBUF,INDX,
     &                  NGROUP,VS,SCR,MXPRM2)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Evaluates the AO integrals of -r^2 r_c (c = x,y,z), the trace  &
C &  of the third moment, for the electric dipole/octupole term of  &
C &  the beyond-dipole oscillator strength (the OCTUP* integrals    &
C &  are traceless). Same structure as magqmom.F.                  &
C &  Ajith Perera, 11/2026                                          &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND
      CHARACTER*8 DUML1,LABELD(3),JUNKLAB
C
      COMMON /PRTBUG/ IPRTFS
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      COMMON /NUCINDX/ ICT
      COMMON /PARM/ISL2,ISL3,ISL4,ISLX,THRESH,MTYP,ILNMCM,INT,INUI,
     &             NCENTS,ID20,NSYM,NAO(12),NMO(12),JUNK,PTOT(3,10), 
     &             INTERF,JOBIPH,RUNMOS
      COMMON /HIGHL/ LMNVAL(3,84),ANORM(84)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      COMMON /FILES/LUOUT,MOINTS
C
      DIMENSION VLIST(NOC,10),EXPS(1),CONT(1),LABEL(1),ITRM(1),
     &          ITRN(256,1),CTRN(256,1),NCONTR(1),NPRIM(1),INDX(1)
      DIMENSION IHOLD(10),IADH(10),VS(MXPRM2,10),TEMP(4,3)
      DIMENSION BUF(LBUF,10),IBUF(IINTFP*LBUF,10),AOINT(1)
      DIMENSION CENTR1(3),SCR(1),GBUF(600),IGNDX(600)
      DIMENSION LMN(200),JMN(200)
      DIMENSION A(3),B(3),EP(3),DISTQD(3),AAA(200),BBB(200),
     &          V(10),DV(9)
C     
      DATA LABELD /'R2R_X   ','R2R_Y   ','R2R_Z   '/

      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
C
      LXQ = 0
      THR2 = THRESH/100000.D+00
      IWPB = (LBUF-2)/2
      LENBUF = IINTFP*LBUF
      ZERO = 0.0D+00
      ONE  = 1.00D+00
      TWO  = 2.00D+00
      EIGHT = 8.00D+00
C
C Do some initializations (No nuclear contribution for this term)
C
      CALL SETRHF
      NT = 3
      DO 100 I = 1, NT
         IHOLD(I) = 0
         IBUF(1,I) = 0
         IBUF(2,I) = -1
100   CONTINUE
C
      IEXH  = 0
      ICONT = 0
      ICOFF = 0
C  
C Loop over the distinct group of exponents
C      
      DO 50 I = 1, NGROUP
C     
        ICNT = AND(IBTSHR(LABEL(I), 10), 2**10-1)
        ITYP = IBTAND(LABEL(I), 2**10-1)
C        
        DO 410 II = 1, 3
           A(II) = VLIST(ICNT, II)
 410    CONTINUE

        IE = NPRIM(I)
        ICE = NCONTR(I)
        JEXH = 0
        JCONT = 0
        JCOFF = 0
C     
        DO 49 J = 1, I
C     
           JCNT = IBTAND(IBTSHR(LABEL(J),10), 2**10 - 1)
           JTYP = IBTAND(LABEL(J), 2**10-1)
C
           DO 411 JJ = 1, 3
              B(JJ) = VLIST(JCNT, JJ)
 411       CONTINUE
C 
           JE = NPRIM(J)
           JCE = NCONTR(J)
C
           DO 2 M = 1, 10
              DO 28 MM = 1, MXPRM2
                 VS(MM, M) = ZERO
 28           CONTINUE
 2         CONTINUE
C     
           IEX = IEXH
C
C Loop over the prmitives
C 
           DO 40 II = 1, IE
              IEX = IEX + 1 
              JEX = JEXH
              DO 41 JJ = 1, JE
                 JEX = JEX + 1
                 IJIND = IE*(JJ-1) + II
C     
C Get the common center of gaussian on A and on B. The resulting
C center remains the same for new gaussians resulting from the
C differentiation
C
                 CALL CIVPT(A, EXPS(IEX), B, EXPS(JEX), EP, GAMAP,
     &                      EFACTP)
C     
C Get the angular momentum quantum numbers of the gaussian
C centered on B
C
                 LB = LMNVAL(1, JTYP)
                 MB = LMNVAL(2, JTYP)
                 NB = LMNVAL(3, JTYP)
C
C Raw third moments of the primitive pair (B undifferentiated);
C DV(c) = -(r^2 r_c) (electron charge included)
C
                 DO 25 K = 1, 3
                    DISTQD(K) = EP(K) - CENTR1(K)
 25              CONTINUE
                 CALL RHFTCE (AAA, A, EP, ITYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JTYP, JTM, ONE, JMN)
                 CALL DZERO(DV,3)
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 211
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 212
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
                       CALL THR_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C V: xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
                       DV(1) = DV(1) - EFACTP*DD*(V(1)+V(6)+V(8))
                       DV(2) = DV(2) - EFACTP*DD*(V(4)+V(2)+V(9))
                       DV(3) = DV(3) - EFACTP*DD*(V(5)+V(7)+V(3))
  212                  CONTINUE
                    ENDDO
  211               CONTINUE
                 ENDDO
C
                 DO 30 MM = 1, 3
                    VS(IJIND, MM) = DV(MM) + VS(IJIND, MM)
 30              CONTINUE
C
#ifdef _DEBUG_LVLM1
                 write(luout,*)
                 write(luout,*) 'Center1  ', i ,'     ',ITYP
                 write(luout,*) 'Center2  ', j ,'     ',JTYP
                 Write (luout, 999) (VS(IJIND, IIII), IIII = 1 ,10)
#endif
C
 41           CONTINUE
 40        CONTINUE
C     
           DO 101 M = 1, NT
C
C PERFORM GENERAL CONTRACTION, SPECIAL CODE FOR 1X1 CASE
C
              IF (IE*JE .EQ. 1) THEN
                 VS(1,M) = VS(1,M)*CONT(ICONT+1)*CONT(JCONT+1)
              ELSE
                 CALL MXM (VS(1,M), IE, CONT(JCONT+1), JE, SCR, JCE)
                 CALL MXMT (CONT(ICONT+1), ICE, SCR, IE, VS(1,M), JCE)
              ENDIF
C&&&&&&&&&DEBUG
C              Write (luout, 999) (VS(1,IIII), IIII = 1 ,3)
C&&&&&&&&&&&&          
              DO 111 II = 1, ICE
                 IIND = INDX(ICOFF + II)
                 JTOP = JCE
                 IF (I .EQ. J) JTOP = II
C
                 DO 112 JJ = 1, JTOP
                    JJND = INDX(JCOFF + JJ)
                    VSINT = VS(ICE*(JJ-1)+II, M)
                    IF(ABS(VSINT) .LT. THRESH) GO TO 112
                    IHOLD(M) = IHOLD(M) + 1
                    IPUT = IHOLD(M) + 2
                    BUF(IPUT, M) = VSINT
                    IOFF = IPUT + IWPB
                    IBUF(IOFF*IINTFP, M) = IBTOR(IBTSHL(MAX0(IIND,JJND)
     &                                    ,16), MIN0(IIND,JJND))
                    IF(IHOLD(M).LT.IWPB) GO TO 112
                    IBUF(1, M) = IWPB
                    CALL OUTAPD(IBUF(1,M), LENBUF, ID20, LXQ)
                    IBUF(2, M) = LXQ
                    IHOLD(M) = 0
 112             CONTINUE
 111          CONTINUE
C     
 101       CONTINUE
           JEXH = JEXH + JE
           JCONT = JCONT + JE*JCE
           JCOFF = JCOFF + JCE
C
 49     CONTINUE
        IEXH = IEXH + IE
        ICONT = ICONT + IE*ICE
        ICOFF = ICOFF + ICE
 50   CONTINUE
C     
C     CLOSE BINS
C     
C
      DO 102 I=1, NT
         IF(IHOLD(I).EQ.0) GO TO 103
         IBUF(1,I) = IHOLD(I)
         CALL OUTAPD(IBUF(1,I), LENBUF, ID20, LXQ)
         IADH(I) = LXQ
         GO TO 102
 103     CONTINUE
C
         IADH(I) = IBUF(2,I)
 102  CONTINUE
C
C     NOW TRANSFORM TO SYMMETRY INTS
C
      NU2 = NSABF*(NSABF+1)/2
      NU = NU2 + 3
      DO 600 NX = 1, NT
C
         DO 70 I = 1, NU2
            AOINT(I) = ZERO
 70      CONTINUE
C
         LL=IADH(NX)
 75      CONTINUE
         IF(LL.EQ.-1) GO TO 76
         CALL INTAPD(IBUF, LENBUF, ID20, LL)
         LOOP = IBUF(1, 1)
         LL = IBUF(2, 1)
         DO 74 IQ = 1, LOOP
            X = BUF(IQ + 2, 1)
            IOFF = IQ + 2 + IWPB
            I = IBTAND(IBTSHR(IBUF(IOFF*IINTFP,1), 16), 2**16 - 1)
            J = IBTAND(IBUF(IOFF*IINTFP, 1), 2**16 - 1)
            IE = ITRM(I)
            DO 702 II =1,IE
               CC = CTRN(II,I)
               IS = ITRN(II,I)
               JE = ITRM(J)
               DO 704 JJ = 1, JE
                  JS = ITRN(JJ,J)
                  IF (IS .LT. JS) GO TO 704
                  XP = CC*CTRN(JJ,J)*X
                  IJF = IS*(IS-1)/2+JS
                  AOINT(IJF) = XP + AOINT(IJF)
 704           CONTINUE
 702        CONTINUE
C     
            IF(I .EQ. J) GO TO 74
            DO 722 II = 1, JE
               CC = CTRN(II, J)
               IS = ITRN(II, J)
               DO 724 JJ = 1, IE
                  JS = ITRN(JJ, I)
                  IF(IS .LT. JS) GO TO 724
                  XP = CC*CTRN(JJ, I)*X
                  IJF = IS*(IS-1)/2+JS
                  AOINT(IJF)=XP+AOINT(IJF)
 724           CONTINUE
 722        CONTINUE
 74      CONTINUE
C     
         GO TO 75
 76      CONTINUE
C     
C No VPOUT ("CRAPS interface") copy; JOBARC only.
C
         CALL PUTREC(20, 'JOBARC', LABELD(NX), NU2*IINTFP, AOINT)
C
C      
         IF(ISL2.EQ.0) GO TO 600
C     CALL PRINTM(AOINT,NCONT)
C     
 600  CONTINUE
C
      TI = T2 - T1
      TS = T3 - T2
      TR = T3 - T1
C
C     WRITE(6,530) TI, TS, TR
 530  FORMAT(' INTS ',F8.3,' SYMMETRY ',F8.3,' TOTAL ',F8.3)
C     CLOSE(26)
      RETURN
      END
