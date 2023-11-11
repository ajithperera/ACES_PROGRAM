c THIS ROUTINE LOADS THE AO INTEGRALS FROM THE CORRESPONDING
c INTEGRAL FILE AND CONSTRUCTS THE FOCK MATRIX (ARRAY F) USING
c THE DENSITY MATRIX SUPPLIED IN D.
c
c COLOUMB TERMS:
c
c    F(I,J) = F(I,J) + D(K,L) (IJ|KL)
c    F(K,L) = F(K,L) + D(I,J) (IJ|KL)
c
c EXCHANGE TERMS:
c
c    F(I,K) = F(I,K) - 1/4 D(J,L) (IJ|KL)
c    F(J,L) = F(J,L) - 1/4 D(I,K) (IJ|KL)
c    F(I,L) = F(I,L) - 1/4 D(J,K) (IJ|KL)
c    F(J,K) = F(J,K) - 1/4 D(I,L) (IJ|KL)
c
c  J. GAUSS AND J.F. STANTON, QTP 1993
c  A. YAU AND A. PERERA, 2000 - MOLCAS integral file support

      SUBROUTINE MKRHFF(F,D,H,BUF,IBUF,NTOTAL,NBAST,NBAS,
     &                  IMAP,ILNBUF,LUINT,ADDH,
     &                  naobasfn,iuhf,scfks,scfksexact,
     &                  scfkslastiter,iter,
     &                  V,z1,ksa,
     &                  screxc,scr,scr2,
     &                  screval,valao,valgradao,totwt,
     &                  max_angpts,natoms,intnumradpts,ncount,kshf)

c   F - SYMMETRY BLOCKED FOCK MATRIX
c   D - SYMMETRY BLOCKED DENSITY MATRIX
c   H - SYMMETRY BLOCKED ONE-ELECTRON INTEGRALS
c
c   SEWARD:
c      BUF  - square-packed, unscaled density matrix
c             (vs. scaled, triangular-packed in D)
c      IBUF - work buffer for doubles, used like the array Work in
c             MOLCAS/ftwo.f
c             BEWARE: The compiler addresses this array as integers, so
c                     all the offset counters must be scaled by iintfp.
c                     In addition, you cannot (meaningfully) address any
c                     one element in the array, you must pass the
c                     address to a subroutine.
c      IMAP - NOT USED!!!
c      ILNBUF - maximum number of doubles that will fit in IBUF
c               (plug-in for GetMem(...,'MAX',...))
c   .NOT.SEWARD:
c      BUF  - BUFFER OF SIZE ILNBUF*IINTFP FOR INTEGRALS
c      IBUF - BUFFER OF SIZE ILNBUF FOR INDICES
c      IMAP - SYMMETRY VECTOR (FULL ARRAY -> SYMMETRY BLOCKED)
c      ILNBUF - BUFFER LENGTH
c
c   NTOTAL  - TOTAL LENGTH OF D AND F
c   NBAST   - TOTAL NUMBER OF BASIS FUNCTIONS
c   NBAS(I) - NUMBER OF BASIS FUNCTIONS IN IRREP I
c
c   LUINT - UNIT NUMBER FOR INTEGRAL FILE
c
c   ADDH - LOGICAL SWITCH FOR ADDING H TO F

c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
C
C-----------------------------------------------------------------     
c ARGUMENT LIST
      INTEGER NTOTAL, NBAST, NBAS(8), ILNBUF, LUINT
      DOUBLE PRECISION F(NTOTAL), D(NTOTAL), H(NTOTAL), BUF(ILNBUF)
      INTEGER IBUF(ILNBUF), IMAP(NBAST,NBAST)
      LOGICAL ADDH

C--------------------------------------------------------------
C KS related decalartions, prakash 23/04
#include<aces.h>
#include<intgrtflags.com>
C
       integer iuhf,iter,naobasfn,ncount,natoms,max_angpts,
     &         intnumradpts
       logical scfks,scfksexact,scfkslastiter,kshf
C
       integer z1(naobasfn,2)
C
       double precision V(naobasfn,naobasfn,iuhf+1),
     &                  ksa(NTOTAL),
     &                  screxc(naobasfn,naobasfn),
     &                  scr(naobasfn,naobasfn,iuhf+1),
     &                  scr2(nbast,nbast),
     &                  screval(nbast),coef_nonloc, 
     &                  valao(naobasfn,max_angpts,
     &                        intnumradpts,ncount),
     &                  valgradao(naobasfn,max_angpts,
     &                            intnumradpts,ncount,3),
     &                  totwt(ncount,intnumradpts,max_angpts)
C-----------------------------------------------------------------
c INTERNAL VARIABLES FOR/FROM VMOL
      DOUBLE PRECISION X
      INTEGER I, J, K, L, IOFF, IOFFT
      INTEGER IJ, IK, IL, JK, JL, KL, INDI, INDJ, INDK, INDL
      INTEGER IRREP, INDS, ILENGTH
      INTEGER NAOBUF, NUMINT, NUT
      CHARACTER*80 FNAME
c INTERNAL VARIABLES FOR/FROM MOLCAS
      double precision DDot, Temp
      integer nSym,
     &        iSym, kSym, nBs, iBs, jBs, kBs, mBs,
     &        iRc, iOpt, lbuf, nSub,
     &        ijPairs, klPairs, nInts, lw2,
     &        ipInts, ipTri, ipSqr, jpTri, jpSqr, kpTri, kpSqr

c COMMON BLOCKS
#include "molcas.com"
#include "machsp.com"
#include "syminf.com"
#include "flags.com"
cwc1
#include "icore.com"
      common /GMSdirect/ dirscf,fdiff,schwrz,pople,i011,i012,i013
      logical dirscf,fdiff,schwrz,pople
      integer i011,i012,i013
cwc0

c PARAMETERS
      DOUBLE PRECISION HALFM, FOURTHM, ZERO, HALF, ONE, TWO
      PARAMETER (HALFM=-0.5D0, FOURTHM=-0.25D0)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)

c FUNCTION DEFINITIONS
      INTEGER IAND,IOR,ISHFT
      INTEGER INT, IUPKI, IUPKJ, IUPKK, IUPKL, IPACK, INDX
      IUPKI(INT)=IAND(INT,IALONE)
      IUPKJ(INT)=IAND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=IAND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=IAND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=IOR(IOR(IOR(I,ISHFT(J,IBITWD)),
     &                             ISHFT(K,2*IBITWD)),
     &                             ISHFT(L,3*IBITWD))
      INDX(I,J)=J+(I*(I-1))/2
C
C---------------------------------------------------------------------
c prakash 23/04
      if ( scfks) then
         call dzero(ksa,ntotal)
      end if
c ----------------------------------------------------------------------
c initialize IMAP
       CALL IZERO(IMAP,NBAST*NBAST)
       IOFF=0
       IOFFT=0
       DO IRREP=1,NIRREP
          DO I=1,NBAS(IRREP)
             DO J=1,NBAS(IRREP)
                INDS=INDX(MAX(I,J),MIN(I,J))+IOFFT
                IMAP(I+IOFF,J+IOFF)=INDS
             END DO
          END DO
         IOFF=IOFF+NBAS(IRREP)
          IOFFT=IOFFT+NBAS(IRREP)*(NBAS(IRREP)+1)/2
       END DO

c initialize the Fock matrix

       DO I=1,NTOTAL
          F(I)=ZERO
       END DO
c ----------------------------------------------------------------------
ccccccccccccccccccccccccc
c START VMOL PROCESSING c
ccccccccccccccccccccccccc

c IIII CONTRIBUTES TO BOTH COLOUMB AND EXCHANGE INTEGRALS
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      NAOBUF=0
      NUMINT=0
      CALL LOCATE(LUINT,'TWOELSUP')
 1    CONTINUE
         READ(LUINT) BUF, IBUF, NUT
         NAOBUF=NAOBUF+1
         DO INT=1,NUT
            X=BUF(INT)
            INDI=IUPKI(IBUF(INT))
            INDJ=IUPKJ(IBUF(INT))
            INDK=IUPKK(IBUF(INT))
            INDL=IUPKL(IBUF(INT))
            IJ=IMAP(INDJ,INDI)
            KL=IMAP(INDL,INDK)
            IK=IMAP(INDK,INDI)
            JL=IMAP(INDL,INDJ)
            IL=IMAP(INDL,INDI)
            JK=IMAP(INDJ,INDK)
            IF (INDI.EQ.INDJ) X=X*HALF
            IF (INDK.EQ.INDL) X=X*HALF
            IF (IJ.EQ.KL)     X=X*HALF
            F(IJ)=F(IJ)+D(KL)*X
            F(KL)=F(KL)+D(IJ)*X

C prakash 23/04
C For B3LYP type DFT we need a paramter for non local exhange 
C coef_pot_nonlocal<=coef_pot_exch(fun_exch_hf)
C fun_exch_hf<=b3lypa
C b3lypa< dftfunc.com as a parameter
C
               if (scfks.AND..not.scfkslastiter) then
                  if(kshf) then
                     coef_pot_nonlocal=0.d0
                  end if      
                  if(coef_pot_nonlocal .ne.0.d0)then
                     coef_nonloc=coef_pot_nonlocal  
                  else
                     goto 50
                  end if
               else
                  coef_nonloc=1.D0
               end if
                                
            F(IK)=F(IK)+coef_nonloc*FOURTHM*D(JL)*X
            F(JL)=F(JL)+coef_nonloc*FOURTHM*D(IK)*X
            F(IL)=F(IL)+coef_nonloc*FOURTHM*D(JK)*X
            F(JK)=F(JK)+coef_nonloc*FOURTHM*D(IL)*X

  50     continue

      END DO
      NUMINT=NUMINT+NUT
      IF (NUT.NE.-1) GOTO 1
      CLOSE(UNIT=LUINT,STATUS='KEEP')

c FOR NIRREP > 1, READ ALSO IIJJ AND IJIJ INTEGRAL FILES
      IF (NIRREP.GT.1) THEN

c IIJJ CONTRIBUTES ONLY TO COLOUMB INTEGRALS
         CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
         OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &        ACCESS='SEQUENTIAL')
         NAOBUF=0
         NUMINT=0
         CALL LOCATE(LUINT,'TWOELSUP')
 101     CONTINUE
            READ(LUINT) BUF, IBUF, NUT
            NAOBUF=NAOBUF+1
            DO INT=1,NUT
               X=BUF(INT)
               INDI=IUPKI(IBUF(INT))
               INDJ=IUPKJ(IBUF(INT))
               INDK=IUPKK(IBUF(INT))
               INDL=IUPKL(IBUF(INT))
               IJ=IMAP(INDJ,INDI)
               KL=IMAP(INDL,INDK)
               IF (INDI.EQ.INDJ) X=X*HALF
               IF (INDK.EQ.INDL) X=X*HALF
               F(IJ)=F(IJ)+D(KL)*X
               F(KL)=F(KL)+D(IJ)*X
            END DO
            NUMINT=NUMINT+NUT
         IF (NUT.NE.-1) GOTO 101
         CLOSE(UNIT=LUINT,STATUS='KEEP')

c IJIJ CONTRIBUTES ONLY TO EXCHANGE
         CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
         OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &        ACCESS='SEQUENTIAL')
         NAOBUF=0
         NUMINT=0
         CALL LOCATE(LUINT,'TWOELSUP')
 201     CONTINUE
            READ(LUINT) BUF, IBUF, NUT
            NAOBUF=NAOBUF+1
            DO INT=1,NUT
               X=BUF(INT)
               INDI=IUPKI(IBUF(INT))
               INDJ=IUPKJ(IBUF(INT))
               INDK=IUPKK(IBUF(INT))
               INDL=IUPKL(IBUF(INT))
               IJ=INDX(INDJ,INDI)
               KL=INDX(INDL,INDK)
               IK=IMAP(INDK,INDI)
               JL=IMAP(INDL,INDJ)
               IF (IJ.EQ.KL) X=X*HALF
C____________________________________________________________________

                if (scfks.AND. .not. scfkslastiter) then
                    if(kshf) then
                       coef_pot_nonlocal=0.d0
                    end if
                    if (coef_pot_nonlocal .ne.0d0) then
                       coef_nonloc=coef_pot_nonlocal
                    else 
                       goto 60
                    end if
                else
                   coef_nonloc=1.d0
                end if
               F(IK)=F(IK)+coef_nonloc*FOURTHM*D(JL)*X
               F(JL)=F(JL)+coef_nonloc*FOURTHM*D(IK)*X

   60       continue
C___________________________________________________________________
         END DO
         NUMINT=NUMINT+NUT
         IF (NUT.NE.-1) GOTO 201
         CLOSE(UNIT=LUINT,STATUS='KEEP')

c     END IF (NIRREP.GT.1)
      END IF

c SCALE DIAGONAL PARTS OF F BY A FACTOR OF TWO
      CALL SCALEF(F,NTOTAL,NBAS)
      CALL XSCAL(NTOTAL,TWO,F,1)
C-----------------------------------------------------------------
C ks terms
      if(scfks.AND. .not. scfkslastiter) then
         if(iter.gt.1)then
           call dzero(V,naobasfn*naobasfn*(1+iuhf))
           call integxc(V,z1,screxc,scr,scfkslastiter,
     &                       valao,valgradao,totwt,
     &        intnumradpts,max_ANGPTS,ncount,kshf)
           call mat_trans(MAT_FULL,MAT_TRI,v(1,1,1),ksa,0)
           call B_AXPY(NTOTAL,one,ksa,1,F,1)
         end if
      end if
C end prakash
C________________________________________________________________
C 
c ADD THE 1-ELECTRON CONTRIBUTION
      IF (ADDH) CALL XAXPY(NTOTAL,ONE,H,1,F,1)
#ifdef _PRINTCRAP
C      write(*,*) "The fock matrix:"
C      call realprt(F,ntotal)
CSS      Write(6, "(6(1x,F10.5))") (H(i),i=1,ntotal)
      call checksum("FA", F, ntotal)
      call checksum("DA", D, ntotal)
#endif /* _PRINTCRAP */

ccccccccccccccccccccccc
c END VMOL PROCESSING c
ccccccccccccccccccccccc

      RETURN
      END

