










      Subroutine Track_vrt_ntos(CNtos,Ntos,Cmos,Dvv,Numb,Nbas,
     +                          ICOLL,IMAP_A,IMAP_B,Lunitn)

      Implicit Double Precision (A-H, O-Z)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Double Precision Ntos
      Integer Track

      Dimension CNtos(Nbas,Numb), Cmos(Nbas,Numb),Dvv(Numb,Numb)
      Dimension Ntos(Numb,Numb)
      Dimension IMAP_A(*),IMAP_B(*)
C     Dimension Ovlp(Maxbasfn,Maxbasfn), Track(Maxbasfn)
      Dimension Ovlp(Numb,Numb), Track(Numb)

      Do i = 1, Numb
         Do j = 1, Numb
C           Ovlp(j,i) = Ddot(Nbas, CNtos(1,i), 1, Cmos(1,j), 1)
            Ovlp(j,i) = Ntos(j,i)
         Enddo 
      Enddo 
      
C
C For a given i (NTO) find the MO that has the largest overlap
C
      Do i = 1, Numb
         Test = 0.0D0
         Do j = 1, Numb 
            If (Dabs(Ovlp(j,i)) .Gt. Test) Then
                Test = Dabs(Ovlp(j,i))
C               Track(i) = j + Nocco(1)
                Track(i) = IMAP_A(ICOLL+j-1)
            Endif 
         Enddo
      Enddo 
C     Do i = 1, Numb
C        Test = 0.0D0
C        Do j = 1, Numb
C            If (Dabs(Ntos(j,i)) .Gt. Test) Then
C                Test = Dabs(Ntos(j,i))
C                Track(i) = j + Nocco(1)
C                Track(i) = IMAP_A(ICOLL+j-1) 
C            Endif
C        Enddo
C     Enddo

      Do i = 1, Numb
         WRITE(LUNITN,'(A,I7,1X,F15.10)') "VRT ",Track(i),
     +                                     2.0D0*Dvv(i,i)
      Enddo
C
      Return
      End
