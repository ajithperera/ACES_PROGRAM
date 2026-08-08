      Subroutine Get_external_occnums(Occ_numsa,Occ_numsb,Nbfns_int,
     +                                Frac_Occnum_present,Me)

      Implicit None
      Integer Nalpha,Nbeta
      Logical Frac_Occnum_present
      Integer Iunit, I,Me
      Double Precision Docc_num,Nbfns_ext,Nbfns_int
      Double Precision Occ_numsa(Nbfns_int),Occ_numsb(Nbfns_int)
      Character*80 Blank
      
      Data Iunit /10/

      Inquire(File="FRAC_OCCNUMS",EXIST=FRAC_OCCNUM_PRESENT)

      If (FRAC_OCCNUM_PRESENT) Then
          If (me .eq. 0) Then
          Write(6,*)
          Write(6,"(a,a)") " The occupation numbers are read",
     +                     " from the FRC_OCCNUMS file."
          Write(6,*)
          Endif 
          Open(Unit=IUNIT, FILE="FRAC_OCCNUMS", FORM="FORMATTED")

          Read (Iunit,"(I4)",End=19) Nbfns_ext
          Read (Iunit,"(I4)",End=19) Nalpha 
          Read (Iunit,"(I4)",End=19) Nbeta 
          Read (Iunit,"(80a)") Blank
  
          If (Nbfns_ext .NE.  Nbfns_int) Then
             If (Me .eq. 0) Then
             Write(6,"(3a,I3,a,I3)") " Inconsistent value for the",
     +                               " number of basis functions!", 
     +                               " External: ",Nbfns_ext,
     +                               " and internali: ",Nbfns_int
             Endif 
             Call abort_job() 
          Endif 

          Do I = 1, Nalpha  
             Read(Iunit,"(F10.5)",End=19) Docc_num
             Occ_numsa(I) = Docc_num
          Enddo
          Do I = Nalpha+1, Nbfns_ext
             Occ_numsa(I) = 0.0D0
          Enddo 

          Read (Iunit,"(80a)") Blank

          Do I = 1, Nbeta
             Read(Iunit,"(F10.5)",End=19) Docc_num
             Occ_numsb(I) = Docc_num
          Enddo 
          Do I = Nbeta+1, Nbfns_ext
             Occ_numsa(I) = 0.0D0
          Enddo 

      Endif 

  19  Close(Iunit)

      Return
      End 
