










      Subroutine Track_occ_ntos(CNtos,Ntos,Cmos,Doo,Numi,Nbas,
     +                          ICOLL,IMAP_A,IMAP_B,Lunitn)

      Implicit Double Precision (A-H, O-Z)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end

      Double Precision Ntos
      Integer Track

      Dimension CNtos(Nbas,Numi), Cmos(Nbas,Numi),Doo(Numi,Numi)
      Dimension Ntos(Numi,Numi)
      Dimension IMAP_A(*),IMAP_B(*)
C     Dimension Ovlp(Maxbasfn,Maxbasfn), Track(Maxbasfn)
      Dimension Ovlp(Numi,Numi), Track(Numi)

      Do i = 1, Numi
         Do j = 1, Numi
C           Ovlp(j,i) = Ddot(Nbas, CNtos(1,i), 1, Cmos(1,j), 1)
            Ovlp(j,i) = Ntos(j,i)
         Enddo 
      Enddo 
      
C
C For a given i (NTO) find the MO that has the largest overlap
C
      Do i = 1, Numi
         Test = 0.0D0
         Do j = 1, Numi 
             If (Dabs(Ovlp(j,i)) .Gt. Test) Then
                 Test = Dabs(Ovlp(j,i))
C                Track(i) = j
                 Track(i) = IMAP_A(ICOLL+j-1)
             Endif 
         Enddo
      Enddo 
      Do i = 1, Numi
         WRITE(LUNITN,'(A,I7,1X,F15.10)') "OCC ",Track(i),
     +                                     2.0D0*Doo(i,i)
      Enddo 
C
      Return
      End
