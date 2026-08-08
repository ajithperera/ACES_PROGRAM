










      SUBROUTINE AUGCNT(W1,W2,N)
C
C....    INCREMENT COUNT OF CONTRACTED BASIS FUNCTIONS
C....    THESE ARE STORED IN THE 16 CHARACTER POSITIONS OF W1//W2,
C....    AS I2,1X,I2,1X,I2,1X,I1,1X,I1,1X,I1,1X,I1.
C
C....    IT IS NOT ASSUMED THAT THESE WORDS ARE ADJACENT IN MEMORY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      character*8 w1,w2,buf
      DIMENSION BUF(2), IANG(7)
      MAXAQN = MAX(MAXAQN,N)
c      BUF(1) = W1
c      BUF(2) = W2
c      DECODE(16,1000,BUF) IANG
C
      read(w1,'(i2,1x,i2,1x,i2)')(iang(j),j=1,3)
      read(w2,'(1x,i2,1x,i1,1x,i1,i1)')(iang(j),j=4,7)

      IANG(N) = IANG(N) + 1

      read(w1,'(i2,1x,i2,1x,i2)')(iang(j),j=1,3)
      read(w2,'(1x,i2,1x,i1,1x,i1,i1)')(iang(j),j=4,7)
c      ENCODE(16,1000,BUF) IANG
c      write(buf,1000)iang
c      W1 = BUF(1)
c      W2 = BUF(2)
      RETURN
1000  FORMAT(I2,1X,I2,1X,I2,1X,I1,1X,I1,1X,I1,1X,I1)
      END
