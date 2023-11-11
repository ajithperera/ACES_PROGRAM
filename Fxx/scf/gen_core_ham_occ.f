











      SUBROUTINE GEN_CORE_HAM_OCC(SCR1, SCRN, ONEH, EVAL, SCRA,
     &                            SCR2, SCRB, XFORM, SCRTMP,
     &                            NOCC, LDIM1, LDIM2, IUHF,
     &                            NBAS)

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
c molcas.com : begin
      logical seward, petite_list
      character*8 fnord
      integer luord
      integer ipmat(8,2), lbbt, lbbs
      common /molcas_com/ seward, petite_list, fnord, luord,
     &                    ipmat, lbbt, lbbs
c molcas.com : end
c symm2.com : begin

c This is initialized in vscf/symsiz.

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end


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



C
      DIMENSION SCR1(NBAS*NBAS),SCRN(NBAS*NBAS),ONEH(LDIM1),
     &           XFORM(LDIM1),SCR2(LDIM2),SCRTMP(LDIM2),
     &           SCRA(LDIM2),SCRB(LDIM2), NOCC(16),
     &           EVAL((IUHF+1)*NBAS)
C
      INDX2(I,J,N)=I+(J-1)*N
C
      CALL ZERO(SCR1,NBAS*NBAS)
      CALL ZERO(SCRN,NBAS*NBAS)
C
      ICNT = 0
      DO I=1,NIRREP
         IF(NBFIRR(I).NE. 0) THEN
C
           CALL EXPND2(ONEH(ITRIOF(I)),SCR2,NBFIRR(I))
           CALL EXPND2(XFORM(ITRIOF(I)),SCRA,NBFIRR(I))
           CALL TRANSP(SCRA,SCRB,NBFIRR(I),NBFIRR(I))
           CALL MXM(SCRB,NBFIRR(I),SCR2,NBFIRR(I),SCRTMP,NBFIRR(I))
           CALL MXM(SCRTMP,NBFIRR(I),SCRA,NBFIRR(I),SCR2,NBFIRR(I))
           CALL ZERO(SCR1,NBAS*NBAS)
           CALL EIG(SCR2,SCR1,NBFIRR(I),NBFIRR(I),0)
           DO  J=1,NBFIRR(I)
               ICNT=ICNT+1
               EVAL(ICNT)=SCR2(INDX2(J,J,NBFIRR(I)))
           ENDDO
         ENDIF
C
      ENDDO
C
      IF (IUHF.EQ.1) THEN
         CALL DCOPY(NBAS,EVAL(1),1,EVAL(NBAS+1),1)
      ENDIF
C
      CALL OCCUPY(NIRREP,NBFIRR,NBAS,EVAL,SCR1,IUHF)
C
      WRITE(6,5000)
      WRITE(6,5001)(NOCC(I),I=1,NIRREP)
      WRITE(6,5002)(NOCC(8+I),I=1,NIRREP)
      WRITE(6,5010)
 5000 FORMAT(T3,'Occupancies from core Hamiltonian:',/)
 5001 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5002 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
 5010 FORMAT(/)
C
      RETURN
      END

