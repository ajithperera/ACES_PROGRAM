










      SUBROUTINE EXPORT_AODENS(WORK, DENS, SCR2, LDIM1, LDIM2, MAXCOR,
     &                         NBASX, NBAS, IUHF)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)



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
C
      DIMENSION WORK(MAXCOR), DENS((IUHF+1)*LDIM1), SCR2(LDIM2)
C
      COMMON /FLAGS/  IFLAGS(100)
C
C  Read AO density from AODENS file.
C
      OPEN(UNIT=61,FILE='AODENS',STATUS='OLD',ACCESS='SEQUENTIAL',
     &     FORM='FORMATTED')
C
      I000 = 1
      I010 = I000 + NBASX*NBASX
      I020 = 1010 + NBASX*NBASX
      I030 = I020 + NBASX*NBAS 
      IEND = I030 + NBASX*NBASX
C
      CALL GETAODENS(WORK(I000), 61, NBASX)
      IF (IUHF .EQ. 1)  CALL GETAODENS(WORK(I010), 61, NBASX)
        
      CLOSE(61,STATUS='KEEP')
C
      Write(6,*) "The AO density matrix in EXPORT_GUESS_MO_ORBS"
      call output(work(i000), 1, nbasx, 1, nbasx, nbasx, nbasx, 1)
      if (iuhf .eq. 1) call output(work(i010), 1, nbasx, 1, nbasx,
     &                             nbasx, nbasx, 1)
C
      CALL GETREC(20,'JOBARC','ZMAT2CMP',NBASX*NBAS*IINTFP,
     &             WORK(I020))
      CALL XGEMM('N','N',NBAS,NBASX,NBASX,1.0D+00,WORK(I020),
     &             NBAS,WORK(I000),NBASX,0.0D+00,WORK(I030),NBAS)
      CALL XGEMM('N','T',NBAS,NBAS,NBASX,1.0D+00,WORK(I020),NBAS,
     &            WORK(I030),NBAS,0.0D+00,WORK(I000),NBAS)
      IF (IUHF .EQ. 1) THEN
         CALL XGEMM('N','N',NBAS,NBASX,NBASX,1.0D+00,WORK(I020),
     &               NBAS,WORK(I010),NBASX,0.0D+00,WORK(I030),NBAS)
         CALL XGEMM('N','T',NBAS,NBAS,NBASX,1.0D+00,WORK(I020),NBAS,
     &               WORK(I030),NBAS,0.0D+00,WORK(I010),NBAS)
      ENDIF
C
      Write(6,*) "The sym ada. density matrix in EXPORT_AODENS"
      call output(work(i000), 1, nbas, 1, nbas, nbas, nbas, 1)
      if (iuhf. eq. 1) call output(work(i010), 1, nbas, 1, nbas, nbas, 
     &                             nbas, 1)

      DO I = 1, NIRREP
         CALL GETBLK(WORK(I000), SCR2, NBFIRR(I), NBAS, IREPS(I))
         CALL SQUEZ2(SCR2, DENS(ITRIOF(I)), NBFIRR(I))
C
        Write(6,*) "The symmetry packed density in EXPORT_AODENS"
        call output(scr2, 1, nbfirr(i), 1, nbfirr(i), nbfirr(i),
     &              nbfirr(i), 1)
C
         IF (IUHF .EQ. 1) THEN
            CALL GETBLK(WORK(I010), SCR2, NBFIRR(I), NBAS, IREPS(I))
        Write(6,*) "The symmetry packed density in EXPORT_AODENS"
        call output(scr2, 1, nbfirr(i), 1, nbfirr(i), nbfirr(i),
     &              nbfirr(i), 1)
            CALL SQUEZ2(SCR2, DENS(ITRIOF(I)+LDIM1), NBFIRR(I))
         ENDIF
      ENDDO
c
      RETURN
      END

