










      SUBROUTINE T2MATCH(IF, IFR, JF, JFR, AF, AFR, BF, BFR, T)

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T
C

      PARAMETER (MAX_ACT_SPACE = 50)



      DOUBLE PRECISION T1AS(MAX_ACT_SPACE**2), T2AS(MAX_ACT_SPACE**4)
C
      DIMENSION I1(MAX_ACT_SPACE**2), A1(MAX_ACT_SPACE**2),
     &          I2(MAX_ACT_SPACE**4), A2(MAX_ACT_SPACE**4),
     &          J2(MAX_ACT_SPACE**4), B2(MAX_ACT_SPACE**4)
      DIMENSION IR1(MAX_ACT_SPACE**2), AR1(MAX_ACT_SPACE**2),
     &          IR2(MAX_ACT_SPACE**4), JR2(MAX_ACT_SPACE**4),
     &          AR2(MAX_ACT_SPACE**4), BR2(MAX_ACT_SPACE**4)
C
      COMMON /T1ASPACE/ I1, IR1, A1, AR1, T1AS, NT1SIZE
      COMMON /T2ASPACE/ I2, IR2, J2, JR2, A2, AR2, B2, BR2, T2AS,
     &                  NT2SIZE
    


      DO I = 1, NT2SIZE
            
         IA  = I2 (I)
         IAR = IR2(I)
         JA  = J2 (I)
         JAR = JR2(I)
         AA  = A2 (I)
         AAR = AR2 (I)
         BA  = B2 (I)
         BAR = BR2 (I)
         IF (IA .EQ. IF .AND. JA .EQ. JF .AND. AA .EQ. AF .AND.
     &       BA .EQ. BF .AND. IAR .EQ. IFR .AND. JAR .EQ. JFR
     &      .AND. AAR .EQ. AFR .AND. BAR .EQ. BFR) THEN


             T = T2AS(I)
             RETURN 
         ENDIF
         
      ENDDO
   
      RETURN
      END
      
