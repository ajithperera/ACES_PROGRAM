










      Subroutine Dbls_correctns_2rpa(Work,Length,Iuhf,Rcc_vecl,
     &                               Len_ph_aa,Len_ph_bb,Len_ph_pq,
     &                               E2_rcc_aa,E2_rcc_bb,Listz)
   
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Length),Rcc_vecl(Len_ph_pq)

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
        
C Do the dot product with the left EOM(SF)-CCD vectors to
C compute the second-order correcion.

      I000 = 1
      Ioff = 1
      Do Ispin = 1, Iuhf+1
         If (Ispin .Eq. 1) Then
            I010 = I000 + Len_ph_aa
         Elseif (Ispin .Eq. 2) Then 
            I010 = I000 + Len_ph_bb
         Endif 

         Call Getlst(Work,1,1,1,2+Ispin,Listz)
         Ioff = Ioff + Len_ph_aa * (Ispin - 1)
         If (Ispin .eq. 1) call checksum("CR(a,i)",Work,Len_ph_aa,s)
         If (Ispin .eq. 2) call checksum("CR(b,j)",Work,Len_ph_bb,s)
         If (Ispin .Eq. 1) Then

C This was done to compare with the UHF for debugging purposes.
C            If (Iuhf .EQ. 0) Call Dscal(Len_ph_aa,1.00/Dsqrt(2.0D0),
C     +                                  Rcc_vecl(Ioff),1)

            call checksum("CL(a,i)",Rcc_vecl(Ioff),Len_ph_aa,s)
            E2_rcc_aa = Ddot(Len_ph_aa,Work(I000),1,Rcc_vecl(Ioff),1)

         Else if (Ispin .Eq. 2) Then
       call checksum("CL(b,j)",Rcc_vecl(Ioff),Len_ph_bb,s)
            E2_rcc_bb = Ddot(Len_ph_bb,Work(I000),1,Rcc_vecl(Ioff),1)

         Endif
      Enddo 

      Return
      End
