










C
C This does a very simple thing. It reads alpha spin Fock matrices
C (occ-occ, occ-vrt and vrt-vrt) to find out whether there are non-zero
C off-diagonal elements (note that the list 91, 92 and 93 are without
C diagonal elements. This is required when a user specified NONHF flags
C for references other than ROHF or QRHF. There is the possibility that
C the user truly has non-HF terms or he is just fooling around. The code
C needs to know this to function properly. I am using the integral
C tolerence as the threshold for being non-zero (that make certain sense).
C 09/04 Ajith Perera.
C
      SUBROUTINE ANALYZE_FOCK(FOO, FVV, FVO, NTOCC, NTVRT, NTOCCVRT,
     &                        NONHF_TERMS_EXIST)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD, POP, VRT
      DIMENSION FOO(NTOCC), FVV(NTVRT), FVO(NTOCCVRT)
      LOGICAL NONHF_TERMS_EXIST
C




























































































































































































c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
C
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /INFO/NOCCO(2), NVRTO(2)
C
      CALL GETLST(FOO, 1, 1, 1, 5, 91)
      CALL GETLST(FVV, 1, 1, 1, 5, 92)
      CALL GETLST(FVO, 1, 1, 1, 5, 93)
C
C
      SQRSUM_OCCVRT = DDOT(NTOCCVRT, FVO, 1, FVO, 1)
      SQRSUM_OCCOCC = DDOT(NTOCC, FOO, 1, FOO, 1)
      SQRSUM_VRTVRT = DDOT(NTVRT, FVV, 1, FVV, 1)
C
      IEXP = IFlags2(108)
      TOLRNCE = 10.0D0**(-IEXP)
C
      NONHF_TERMS_EXIST = (SQRSUM_OCCVRT .GT. TOLRNCE .OR.
     &                     SQRSUM_OCCOCC .GT. TOLRNCE .OR.
     &                     SQRSUM_VRTVRT .GT. TOLRNCE)
C
      Write(6,"(a,l)") "NONHF_TERMS_EXIST :", NONHF_TERMS_EXIST
      RETURN
      END
