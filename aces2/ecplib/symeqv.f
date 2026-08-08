










      Subroutine SymEqv (Ntotatoms, IGenBy)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Prepares IGENBY list in a much less cumbersome way
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IMPLICIT INTEGER (A-Z)

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C
      DIMENSION MEMBER(MXATMS),ORBPOP(MXATMS),IGENBY(Ntotatoms)
C
C INITIALIZE GENBY TO 999.  THIS WILL TAKE CARE OF DUMMY ATOMS.
C
      DO 20 I=1,NTOTATOMS
         IGENBY(I)=999
20    CONTINUE
C
C GET MEMBER AND ORBPOP VECTORS FROM JOBARC.
C
      ISIZE=1
      CALL GETREC(20,'JOBARC','COMPNORB',ISIZE,IORBIT)
      CALL GETREC(20,'JOBARC','COMPMEMB',NTOTATOMS,MEMBER)
C
      CALL GETREC(20,'JOBARC','COMPPOPV',IORBIT,ORBPOP)
C
C LOOP OVER ORBITS AND ZERO OUT POSITION OF FIRST MEMBER OF EACH ORBIT
C   IN GENBY LIST.
C
      IOFF=1
      DO 50 I=1,IORBIT
         IGENBY(MEMBER(IOFF))= I 
         DO 51 J=1,ORBPOP(I)-1
            IGENBY(MEMBER(IOFF+J))=I 
51       CONTINUE
         IOFF=IOFF+ORBPOP(I)
50    CONTINUE
      
      RETURN
      END
