










      SUBROUTINE PRIGRD(GRAD,ISAVE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)

c IMPORTANT: The ECP library contains three routines (doitgr, ecpabl, ecpder)
c that use this parameter but do not include this file. Changing this value
c requires updating those files as well. [Yau: If we linked them together,
c then compiling ecp would require the vdint directory. Alternatively, we
c could put the mxcent.par file into the top-level include directory.]

C MXCENT : Maximum number of atoms currently allowed
C MXCOOR : The number of Cartesian that correspond to MXCENT

      INTEGER MXCENT, MXCOOR
      PARAMETER (MXCENT=200, MXCOOR=3*MXCENT)

c These parameters are gathered from vmol and vdint and are used by ecp
c as well. It just so happens that the vmol parameters do not exist in
c vdint and vice versa. LET'S TRY TO KEEP IT THAT WAY!

c VMOL PARAMETERS ------------------------------------------------------

C     MAXPRIM - Maximum number of primitives for a given shell.
      INTEGER    MAXPRIM
      PARAMETER (MAXPRIM=72)

C     MAXFNC  - Maximum number of contracted functions for a given shell.
C               (vmol/readin requires this to be the same as MAXPRIM)
      INTEGER    MAXFNC
      PARAMETER (MAXFNC=MAXPRIM)

C     NHT     - Maximum angular momentum
      INTEGER    NHT
      PARAMETER (NHT=7)

C     MAXATM  - Maximum number of atoms
      INTEGER    MAXATM
      PARAMETER (MAXATM=100)

C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers.
      INTEGER    MXTNPR
      PARAMETER (MXTNPR=MAXPRIM*MAXPRIM)

C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers.
      INTEGER    MXTNCC
      PARAMETER (MXTNCC=180000)

C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers.
      INTEGER    MXTNSH
      PARAMETER (MXTNSH=200)

C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions).
c mxcbf.par : begin

c MXCBF := the maximum number of Cartesian basis functions (limited by vmol)

c This parameter is the same as MAXBASFN. Do NOT change this without changing
c maxbasfn.par as well.

      INTEGER MXCBF
      PARAMETER (MXCBF=1000)
c mxcbf.par : end

c VDINT PARAMETERS -----------------------------------------------------

C     MXPRIM - Maximum number of primitives for all symmetry
C              inequivalent centers.
      INTEGER    MXPRIM
      PARAMETER (MXPRIM=MXTNPR)

C     MXSHEL - Maximum number of shells for all symmetry inequivalent centers.
      INTEGER    MXSHEL
      PARAMETER (MXSHEL=MXTNSH)

C     MXCORB - Maximum number of contracted basis functions.
      INTEGER    MXCORB
      PARAMETER (MXCORB=MXCBF)

C     MXORBT - Length of the upper or lower triangle length of MXCORB.
      INTEGER    MXORBT
      PARAMETER (MXORBT=MXCORB*(MXCORB+1)/2)

C     MXAOVC - Maximum number of subshells per center.
      INTEGER    MXAOVC,    MXAOSQ
      PARAMETER (MXAOVC=32, MXAOSQ=MXAOVC*MXAOVC)

c     MXCONT - ???
      INTEGER    MXCONT
      PARAMETER (MXCONT=MXAOVC)

      CHARACTER NAMEX*6
      CHARACTER*6 NAMES(MXCENT),NAME
      CHARACTER*6 NAMES_REORD(MXCENT)
      LOGICAL DCORD,DCORGD,NOORBT,DOPERT
      LOGICAL FIRST
C
      DIMENSION GRAD(MXCOOR),CGRAD(MXCOOR)
      DIMENSION WMAT(3,3),CGRAD_ROT(MXCOOR)
      DIMENSION IMAP(MXCENT)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON /NUCLEC/ NAMEX(MXCOOR)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
C
      SAVE FIRST, NAMES
C
      DATA FIRST/.TRUE./

      CALL GETREC(20,"JOBARC","REORIENT",9*IINTFP,WMAT)

C
      IF (FIRST) THEN
       FIRST = .FALSE.
       DO 100 IATOM = 1, NUCDEP   
        NAME = NAMEX(3*IATOM)
        NAMES(IATOM) = NAME(1:4)//'  '
100    CONTINUE
      END IF
      IF (MAXLOP .EQ. 0) THEN
       IOFF = 0
       DO 200 IATOM = 1, NUCDEP
        WRITE (LUPRI,1000) ' ', NAMES(IATOM), (GRAD(IOFF+J),J=1,3)
        IF(ISAVE.GT.0) THEN
         WRITE(81,3000) CHARGE(IATOM),(GRAD(IOFF+J),J=1,3)
        ENDIF
        IOFF = IOFF + 3
  200  CONTINUE
       CALL PUTREC(20,'JOBARC','GRADIENT',3*NATOMS*IINTFP,GRAD)

       WRITE (LUPRI,'(/)')
       CALL HEADER('Molecular energy gradient',-1)
       CALL TRAGRD(GRAD,CGRAD,3*NUCDEP)

C Reorder the Cartesian gradients to ZMAT orientation before
C write to the output. This is useful for vibrational correction
C of properties related work. Only the standard out
C reflect this new reordered values. The rest still works in the
C computational order (in order to avoid any internal conflicts)
C Since this block has no symmetry, these transformations and 
C reordering more often than not do anything. Leaving them for
C consistency with the block that use symmetry. 
C Ajith Perera, Jan. 2021.

       CALL GETREC(20,'JOBARC','MAP2ZMAT',NATOMS,IMAP)
       CALL GETREC(20,'JOBARC',"REORIENT",9*IINTFP,WMAT)
       CALL CMP2ZMAT(CGRAD,GRAD,IMAP,NAMES_REORD,NATOMS)
       CALL DGEMM("N","N",3,NATOMS,3,1.0D0,WMAT,3,GRAD,3,0.0D0,
     &            CGRAD_ROT,3)
       KOFF = 0
       DO IATOM = 1, NATOMS
          WRITE (LUPRI,1010) ' ',NAMES_REORD(IATOM)(1:4),IATOM,
     &                          (CGRAD_ROT(KOFF + K), K = 1, 3)
          KOFF = KOFF+3
       ENDDO

      ELSE

       DO 300 I = 1, NCRREP(0)
        WRITE (LUPRI,'(25X,A6,F17.10)') NAMEX(IPTCOR(I)), GRAD(I)
  300  CONTINUE

       WRITE (LUPRI,'(/)')
       CALL HEADER('Molecular energy gradient',-1)
       CALL TRAGRD(GRAD,CGRAD,3*NUCDEP)

C Reorder the Cartesian gradients to ZMAT orientation before 
C write to the output. This is useful for vibrational correction
C of properties related work. Only the standard out
C reflect this new reordered values. The rest still works in the
C computational order (in order to avoid any internal conflicts)
C Ajith Perera, Dec. 2020.

       CALL GETREC(20,'JOBARC','MAP2ZMAT',NATOMS,IMAP)
       CALL GETREC(20,'JOBARC',"REORIENT",9*IINTFP,WMAT)
       CALL CMP2ZMAT(CGRAD,GRAD,IMAP,NAMES_REORD,NATOMS)
       CALL DGEMM("N","N",3,NATOMS,3,1.0D0,WMAT,3,GRAD,3,0.0D0,
     &            CGRAD_ROT,3)

       KOFF = 0
       DO IATOM = 1, NATOMS
          WRITE (LUPRI,1010) ' ',NAMES_REORD(IATOM)(1:4),IATOM,
     &                          (CGRAD_ROT(KOFF + K), K = 1, 3)
          KOFF = KOFF+3
       ENDDO 

       IOFF = 0
       DO 310 IATOM = 1, NUCIND
        MULTI = MULT(ISTBNU(IATOM))
        IF (MULTI .EQ. 1) THEN
CSSS         WRITE (LUPRI,1000) ' ',NAMES_REORD(IATOM),
CSSS     &                          (CGRAD_ROT(IOFF+K),K=1,3)
         IF(ISAVE.GT.0) THEN
          WRITE(81,3000) CHARGE(IATOM),(CGRAD(IOFF+K),K=1,3)
         ENDIF
         IOFF = IOFF + 3
        ELSE
         DO 320 J = 1, MULTI
CSSS          WRITE (LUPRI,1010) ' ',NAMES_REORD(IATOM)(1:4),J,
CSSS     &    (CGRAD_ROT(IOFF + K), K = 1, 3)
          IF(ISAVE.GT.0) THEN
           WRITE(81,3000) CHARGE(IATOM),(CGRAD(IOFF+K),K=1,3)
          ENDIF
          IOFF = IOFF + 3
          JOFF = JOFF + 1
320      CONTINUE
        END IF
310    CONTINUE
       CALL PUTREC(20,'JOBARC','GRADIENT',3*NATOMS*IINTFP,CGRAD)
      END IF
      WRITE (LUPRI,'(/)')
c      IF(ISAVE.GT.0) THEN
c       CLOSE(UNIT=81,STATUS='KEEP')
c      ENDIF
C
CJDW 5/12/95. Write coordinates, energy and gradient to POLYRATE.
C
      IF(IFLAGS2(113).GT.0)THEN
C
      CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &              .TRUE.,.TRUE.,.TRUE.,.FALSE.)
C
      ENDIF
C
      RETURN
1000  FORMAT (A1,A6,F17.10,2F24.10)
1010  FORMAT (A1,A4,I2,F17.10,2F24.10)
3000  FORMAT(4F20.10)
      END
