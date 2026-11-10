      SUBROUTINE SETRES(IRREPX,LISTRS,IUHF)
C
C Shape the resorted-doubles lists written by
C RESORT(...,IRREPX,LIST0,LISTRS) for symmetry IRREPX: LISTRS+3 and
C LISTRS+5 (RHF), plus LISTRS, LISTRS+1, LISTRS+2, LISTRS+4 (UHF).
C Same as the "RESORTED LISTS" block of newlst.F (434-439, 454-459).
C
      IMPLICIT INTEGER (A-Z)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DO I=0,5
         IF (IUHF.NE.0 .OR. I.EQ.0 .OR. I.EQ.3 .OR. I.EQ.5) THEN
            L=LISTRS+I
            CALL NEWTYP2(IRREPX,L,ISYTYP(1,L),ISYTYP(2,L),.TRUE.)
         ENDIF
      ENDDO
      RETURN
      END
