










      SUBROUTINE SCRNC_XXXX_PQXX(ICORE,MAXCOR,IUHF,LSTAOINC,LSTMOINC,
     &                           IRREPX,IFLAG)
C
C This routine drives the transformation of (XX,XX) Coulomb integral
C to (PQ,XX) and write it as (pq,XX)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION X,SNRM2,SDOT,ABLAD,ZILCH,ONE
      LOGICAL TAU,NOINC
      CHARACTER*2 SPLAB(3)
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)

      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
      DATA SPLAB /'AA','BB','AB'/
C
      NNM1O2(I)=(I*(I-1))/2
C
      Write(6,*)
      Write(6,"(a)") " Entered - scrnc_xxxx_pqxx"
      IONE=1
      LUINT=10
      NMO=NOCCO(1)+NVRTO(1)
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
C
C TRANSFORM TARGET VECTOR BACK TO MO REPRESENTATION
C
C (XX,XX) -> (PQ,XX)

      DO ISPIN=1,(IUHF+1)
       ABLAD=ZILCH
       DO IRREPXX=1,NIRREP
        IRREPIJ=DIRPRD(IRREPXX,IRREPX)

        IF (IFLAG .EQ.3) NUMIJ=IRPDPD(IRREPIJ,20+ISPIN)
        IF (IFLAG .EQ.4) NUMIJ=IRPDPD(IRREPIJ,18+ISPIN)
        IF (IFLAG .EQ.5) NUMIJ=IRPDPD(IRREPIJ, 8+ISPIN)
        NUMXX=IRPDPDAO(IRREPXX)
        write(6,*) 
        Write(6,"(a,i2,i2,i4)") " Irrepl,irrepr,numij",IRREPIJ,
     +                           IRREPXX,numij

        I000=1
        I010=I000+NUMXX*NUMXX*IINTFP
        I020=I010+NUMIJ*NUMXX*IINTFP
        I030=I020+MAX(NT(1),NT(2),NAO*NMO)*IINTFP
        I040=I030+MAX(NT(1),NT(2),NAO*NMO)*IINTFP
        I050=I040+MAX(NT(1),NT(2),NAO*NMO)*IINTFP
        IF(I050.GE.MAXCOR) CALL INSMEM('SCRNC_XX_PQXX',I050,MAXCOR)

        CALL GETLST(ICORE(I000),1,NUMXX,1,IRREPXX,LSTAOINC)

        Write(6,*)
CSSS        Call Checksum("I(XX|XX):",ICORE(I000),NUMXX*NUMXX)
CSSS        call output(icore(i000),1,numxx,1,numxx,numxx,numxx,1)
C Diagnostic (2026-08-18, sym=on/off investigation, see
C SCRNC_DEBUG_LOG.md "continued 19") -- plain dump of the raw
C I(XX,XX) block exactly as read here, for Iflag=4 (AB/VV), Irrepxx=1,
C to check directly whether record 17's own row/column is populated in
C THIS routine's own input, before any transform touches it. Routed
C through a properly Double-Precision-typed helper (this routine's own
C "Implicit Integer (A-Z)" would silently bit-reinterpret Icore's
C elements as integers if read directly here -- same trap as the
C already-fixed Jj_lrespn bug, "continued 8").
        If (Iflag .Eq. 4 .And. Irrepxx .Eq. 1 .And. Ispin .Eq. 1) Then
           Call Scrnc_dump_dble_array(Icore(I000),Numxx*Numxx,
     +          "ixxxx_iflag4_irrepxx1.dat",Numxx)
        Endif
        CALL SCRNC_I1TRAN(IFLAG,ICORE(I010),ICORE(I000),ICORE(I020),
     &                   ICORE(I030),ICORE(I040),NAO,NMO,ISPIN,IUHF,
     &                   IRREPXX,IRREPIJ,.FALSE.,"EMPTY")

        CALL PUTLST(ICORE(I010),1,NUMXX,1,IRREPXX,LSTMOINC+ISPIN)

        If (Iflag .eq. 3) then

        write(6,*) "@-scrnc_xxx_pqxx, I(IJ|XX)"
        call checksum("I(IJ,XX):",icore(i010),numij*numxx)

        write(6,"(a)") " I(IJ,XX)"
        call output(icore(i010),1,numij,1,numxx,numij,numxx,1)
C        call transp(icore(i010),icore(i020),Numxx,numij)
C        write(6,"(a)") " I(XX,IJ)"
C        Ioff = 0
C        do irrepr=1,Nirrep
C        irrepl=Dirprd(irrepr,1)
C        numij=irpdpd(irrepr,20+ispin)
C        nrow=iaopop(irrepl)
C        ncol=iaopop(irrepl)
C        do i=1,numij
C        call output(icore(i020+ioff),1,nrow,1,ncol,nrow,ncol,1)
C        ioff = Ioff + nrow*ncol
C        enddo
C        enddo 

        elseif (iflag.eq.4) then

        write(6,*) "@-scrnc_xxx_pqxx, I(AB|XX)"
        call checksum("I(AB,XX):",icore(i010),numij*numxx)
        write(6,"(a)") " I(AB,XX)"
        call output(icore(i010),1,numij,1,numxx,numij,numxx,1)
C        call transp(icore(i010),icore(i020),Numxx,numij)
C        Ioff = 0
C        write(6,"(a)") " I(XX,AB)"
C        do i=1,numij
C        call output(icore(i020+ioff),1,nao,1,nao,nao,nao,1)
C        ioff = Ioff + numxx
C        enddo

        elseif (iflag.eq.5) then
        write(6,*) "@-scrnc_xxx_pqxx, I(AI|XX)"
        call checksum("I(AI,XX):",icore(i010),numij*numxx)

        write(6,"(a)") " I(AI,XX)"
        call output(icore(i010),1,numij,1,numxx,numij,numxx,1)
c        call transp(icore(i010),icore(i020),Numxx,numij)
c        Ioff = 0
c        write(6,"(a)") " I(XX,AI)"
c        do i=1,numij
c        call output(icore(i020+ioff),1,nao,1,nao,nao,nao,1)
c        ioff = Ioff + numxx
c        enddo

        endif 
       ENDDO
      ENDDO

      RETURN
      END
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_dump_dble_array(X,N,Fname,Numxx)
C
C Diagnostic helper (2026-08-18, see SCRNC_DEBUG_LOG.md "continued
C 19") -- plain one-value-per-line dump of a Double Precision array,
C called from routines (like Scrnc_xxxx_pqxx.F) that themselves see
C their Work-array argument as Implicit Integer -- passing it here
C forces the correct Double Precision reinterpretation of the same
C underlying memory, since Fortran argument passing is by address.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension X(N)
      Character*(*) Fname
C
      Open(Unit=91,File=Fname,Form="Formatted",Status="UNKNOWN")
      Write(91,*) Numxx
      Do I = 1, N
         Write(91,*) X(I)
      Enddo
      Close(Unit=91)
C
      Return
      End
