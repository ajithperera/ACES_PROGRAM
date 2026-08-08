










      SUBROUTINE DO_ISOTOPES_JJ(IUNITIS,JJFC,NDIM,NATOMS,IEOMPROP)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DOUBLE PRECISION JJFC(NDIM)

      DO ISPCC = 1, 4

         IF (ISPCC.EQ. 1) THEN
            LENGTH = NATOMS*NATOMS
            CALL GETREC(20,"JOBARC","JCC_FC  ",LENGTH,JJFC) 
            IFERMI = 1
            ISDIP  = 0
            IPSO   = 0
            IDSO   = 0
            NTPERT = NATOMS
            ICALL  = 0
         ELSEIF (ISPCC.EQ. 2) THEN
            LENGTH = 9*NATOMS*NATOMS
            CALL GETREC(20,"JOBARC","JCC_PSO ",LENGTH,JJFC)
            IFERMI = 0
            ISDIP  = 0
            IPSO   = 1
            IDSO   = 0
            NTPERT = 3*NATOMS
            ICALL  = 1
         ELSEIF (ISPCC.EQ. 3) THEN
            LENGTH = 36*NATOMS*NATOMS
            CALL GETREC(20,"JOBARC","JCC_SDIP" ,LENGTH,JJFC) 
            IFERMI = 0
            ISDIP  = 1
            IPSO   = 0
            IDSO   = 0
            NTPERT = 6*NATOMS
            ICALL  = 1
         ELSEIF (ISPCC.EQ. 4) THEN
            LENGTH = 9*NATOMS*NATOMS
            CALL GETREC(20,"JOBARC","JCC_DSO " ,LENGTH,JJFC) 
            IFERMI = 0
            ISDIP  = 0
            IPSO   = 0
            IDSO   = 1 
            NTPERT = 3*NATOMS
            ICALL  = 3
         ENDIF 
         REWIND(IUNITIS)
         CALL FACTOR_IS(JJFC, NTPERT, IFERMI, ISDIP, IPSO, IDSO,
     &                  IEOMPROP,ICALL,IUNITIS)
        
      ENDDO
 
      RETURN
      END 

     
