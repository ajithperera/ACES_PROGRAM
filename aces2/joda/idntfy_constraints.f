










      Subroutine Idntfy_constraints(Iredunco, Nredunco, TOTNOFBND, 
     &                              TOTNOFANG, TOTNOFDIH, Not_found)
C
      Implicit Double Precision (A-H, O-Z)
C
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C
      Integer Constrains, Total_constrains, TOTNOFBND, TOTNOFANG,
     &        TOTNOFDIH, Constrained
      Logical Not_found, Sucess
C
      Dimension Iredunco(4, Nredunco), Constrains(4, Maxredunco),
     &          Constrained(Maxredunco)
C 
C Read the constrins from the *CONSTRAINS namelist from the
C ZMAT file.
C
      write(6,"(a)") " Entering read_constrains"
      Call Read_ric_constrains(Constrains,
     &                         Num_bond_Constrains,
     &                         Num_angle_constrains,
     &                         Num_dihs_angle_constrains,
     &                         Not_found)

      If (Not_found) Return

      Write(6,*)
      write(6,"(a,3(1x,I4))") " Constrains:bonds,angles,dihedrals ",
     &                        Num_bond_Constrains, 
     &                        Num_angle_Constrains,
     &                        Num_dihs_angle_Constrains
      Total_constrains = Num_bond_Constrains + Num_angle_Constrains +
     &                   Num_dihs_angle_Constrains

      Write(6,*)
      Print*, "The constrains read from ZMAT"
      Write(6,*)
      Do i = 1, Total_constrains
         Write(6,111) (Constrains(j, i), j=1, 4)
      Enddo
 111  Format(5X, 4(I3, 1X))
C
C If there are no constrains to be found, return and do the full
C optimization. 
C

      Sucess = .False. 
      Ifound = 0
      Call Izero(Constrained, Nredunco)
C
      Do Jc = 1, Num_bond_Constrains
         Do Ic = 1, TOTNOFBND
            
            ICcon = Constrains(1, Jc)
            JCcon = Constrains(2, Jc)

            IRcon = Iredunco(1, Ic)
            JRcon = Iredunco(2, Ic)
            If (ICcon .eq. IRcon .OR. ICcon .eq. JRcon) Then

                If (JCcon .Eq. IRcon .OR. JCcon .eq. JRcon) Then

                    Constrained(Ic) = 1
                    Ifound = Ifound + 1
                    Sucess = .TRUE.

                Endif

            Endif 
            IF (Sucess) GO TO 10

        Enddo
  10    Continue 
        Sucess = .False. 
      Enddo
            
      Ibegin = Num_bond_Constrains + 1
      Jbegin = TOTNOFBND + 1

C A constraint i j k means the angle i-j-k (vertex j). Look for that
C angle first; only if it is not there accept any angle made of the 
C same three atoms (the old rule, which picked the first such angle
C and so could constrain the wrong angle of a three-membered ring).
C 10/2026.
C
      Do Jc = Ibegin, Ibegin + Num_angle_Constrains - 1
         Do Ipass = 1, 2
         Do Ic = Jbegin, Jbegin + TOTNOFANG - 1
            
            ICcon = Constrains(1, Jc) 
            JCcon = Constrains(2, Jc)
            KCcon = Constrains(3, Jc)

            IRcon = Iredunco(1, Ic)
            JRcon = Iredunco(2, Ic)
            KRcon = Iredunco(3, Ic)

            If (Ipass .Eq. 1) Then
               If (JCcon .Eq. JRcon .And. 
     &            ((ICcon .Eq. IRcon .And. KCcon .Eq. KRcon) .Or.
     &             (ICcon .Eq. KRcon .And. KCcon .Eq. IRcon))) Then
                   Constrained(Ic) = 1
                   Ifound = Ifound + 1
                   Sucess = .TRUE.
               Endif
            Else

            If (ICCon .eq. IRcon .OR. ICCon .eq. JRcon .OR. 
     &          ICcon .EQ. KRcon) Then
                If (JCCon .eq. IRcon .OR. JCCon .eq. JRcon .OR.
     &              JCcon .EQ. KRcon) Then
                    If (KCCon .eq. IRcon .OR. KCCon .eq. JRcon .OR.
     &                  KCcon .EQ. KRcon) Then

                        Constrained(Ic) = 1
                        Ifound = Ifound + 1
                        Sucess = .TRUE.

                     Endif
                Endif
            Endif
            Endif

            IF (Sucess) GO TO 20

         Enddo
         Enddo
  20     Continue 
         Sucess = .False. 
      Enddo
C
      Ibegin = Num_bond_constrains + Num_angle_constrains + 1
      Jbegin = TOTNOFBND + TOTNOFANG + 1

C A constraint i j k l means the dihedral i-j-k-l (or l-k-j-i). Look
C for that dihedral first and accept any dihedral of the same four 
C atoms only if it is not there. 10/2026.
C
      Do Jc = Ibegin, Ibegin + Num_dihs_angle_Constrains - 1
         Do Ipass = 1, 2
         Do Ic = Jbegin, Jbegin + TOTNOFDIH - 1

            ICcon = Constrains(1, Jc)
            JCcon = Constrains(2, Jc)
            KCcon = Constrains(3, Jc)
            LCcon = Constrains(4, Jc)

            IRcon = Iredunco(1, Ic)
            JRcon = Iredunco(2, Ic)
            KRcon = Iredunco(3, Ic)
            LRcon = Iredunco(4, Ic)
            If (Ipass .Eq. 1) Then
               If ((ICcon .Eq. IRcon .And. JCcon .Eq. JRcon .And.
     &              KCcon .Eq. KRcon .And. LCcon .Eq. LRcon) .Or.
     &             (ICcon .Eq. LRcon .And. JCcon .Eq. KRcon .And.
     &              KCcon .Eq. JRcon .And. LCcon .Eq. IRcon)) Then
                   Constrained(Ic) = 1
                   Ifound = Ifound + 1
                   Sucess = .TRUE.
               Endif
            Else

            If (ICCon .eq. Ircon .OR. ICcon .Eq. JRcon .OR.
     &          ICcon .eq. Krcon .OR. ICCon .EQ. LRcon) Then
                If (JCCon .eq. Ircon .OR. JCcon .Eq. JRcon .OR.
     &              JCcon .eq. Krcon .OR. JCCon .EQ. LRcon) Then
                    If (KCCon .eq. Ircon .OR. KCcon .Eq. JRcon .OR.
     &                  KCcon .eq. Krcon .OR. KCCon .EQ. LRcon) Then
                        If (LCCon .eq. Ircon .OR. LCcon .Eq. JRcon .OR.
     &                      LCcon .eq. Krcon .OR. LCCon .EQ. LRcon) Then

                            Constrained(Ic) = 1
                            Ifound = Ifound + 1
                            Sucess = .TRUE.

                        Endif 
                    Endif
                Endif
            Endif
            Endif

            IF (Sucess) GO TO 30

         Enddo
         Enddo
  30     Continue
         Sucess = .False.
      Enddo
C
      If (Ifound .Lt. Total_constrains) Then
        Write(6,"(3x,2a)")"At least one constrain specified" ,
     &               " in the input does not correspond to a known"
        Write(6,"(3x,2a)") "RIC. The constrain(s) is(are) ignored",
     &               " and the optimization proceeds."
C
        Write(6,*)
        Write(6,"(3x,2a)") "If this has been a mistake, cancel the job",
     &                " and use the redundent coordinates"
        Write(6,"(3x,a,a)")"printed below to identify the correct ", 
     &                     "degs. of freedoms that can be " 
        Write(6,"(3x,a)")  "constrained."
C
        Write(6,*)
        Write(6,"(3x,a)") "The center lables of redundent coordinates"
        Do i = 1, Nredunco
           Write(6,112) (iredunco(j, i), j=1, 4)
        Enddo
 112    Format(5X, 4(I3, 1X))

      Endif
C
      Write(6,*) 
      Write(6,*) "The CONSTRAINED array"
      Write(6,"(a,I4)") "The total No. of constraints found",Ifound
      Write(6,"(6(1x,I4))") (Constrained(i),i=1,Nredunco)
C
      Call Putrec(20, "JOBARC","CONSTRNS", Nredunco, Constrained)
C
      Return
      End
        
