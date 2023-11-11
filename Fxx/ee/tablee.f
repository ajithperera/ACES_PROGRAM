










      SUBROUTINE TABLEE(EIGVAL,EIGVAL_T,OSCSTR,NATURE,BGN,BGN_IRP,
     &                  END,END_IRP,METHOD)
C
C  THE CALCULATION OF EXCITATION ENERGIES IS SUMMARIZED
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,BGN,END,BGN_IRP,END_IRP
      Logical Triplet
      CHARACTER*1 NATURE(100,8)
      DIMENSION EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),BGN(100,8),
     &          BGN_IRP(100,8),END(100,8),END_IRP(100,8)

      LOGICAL MBPT2, CC,CCD,RCCD,DRCCD,LCCD,LCCSD,EOM_TRPS
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/CALCINFO/NROOT(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD
C
      DATA FACTEV, FACTCM /27.2113957D0,  2.19474625D5/
      DATA ZILCH /0.0D0/
C
C
      Triplet = .FALSE.
      EOM_TRPS = (iflags(87) .EQ. 11)

      WRITE(6,*)
      IF (METHOD .EQ. 3) THEN
        IF (CC .AND. (.NOT. CCD)) THEN
          WRITE(6,*) ' Summary of EOM-CCSD excitation energies'
        ELSE IF (MBPT2) THEN
          WRITE(6,*) ' Summary of EOM-MBPT(2) excitation energies'
        ELSE IF (CCD) THEN
          WRITE(6,*) ' Summary of CCD excitation energies'
        ELSE IF (RCCD) THEN
          WRITE(6,*) ' Summary of rCCD excitation energies'
        ELSE IF (DRCCD) THEN
          WRITE(6,*) ' Summary of drCCD excitation energies'
        ELSE IF (LCCD) THEN
          WRITE(6,*) ' Summary of LCCD excitation energies'
        ELSE IF (LCCSD) THEN
          WRITE(6,*) ' Summary of LCCSD excitation energies'
        ENDIF
      ELSE IF (METHOD .EQ. 7) THEN
        IF (CC) THEN
       WRITE(6,*) ' Summary of P-EOM-CCSD excitation energies'
        ELSE
       WRITE(6,*) ' Summary of P-EOM-MBPT(2) excitation energies'
        ENDIF
      ELSE IF (METHOD .EQ. 8) THEN
        IF (CC) THEN
          WRITE(6,*) ' Summary of BWPT(2)-EOM-CCSD excitation energies'
        ELSE
          WRITE(6,*) ' Summary of BWPT(2)-EOM-MBPT excitation energies'
        ENDIF
      ENDIF

      CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ECC)
      If (EOM_Trps) Then

      Write(6,*)
      Write(6,100)
      Write(6,101)
 100  Format(2x,'Sym',2x,'Origin',4x,'Destination',1x,'EOM-',7x,
     &       'EOM-',8x,'Osc. Str.',4x,'Total Energy')
 101  Format(28x,' CCSD(eV)',3x,'CCSD(T)(eV)')
      Write(6,*)' ----------------------------------------------------',
     &          '-------------------------'
      Write(6,*)
      Write(6,"(2x,a)") "The Singlet Excited States:"
      Write(6,*)
         Do Irrep = 1, Nirrep
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "S") Then
                  Eed   = Eigval(Iroot,Irrep) * Factev
                  Eet   = Eigval_t(Iroot,Irrep) * Factev
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep)

                  Write(6,99) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Eed,Eet,F,Et
               Else
                  Triplet = .True.
               Endif
            Enddo
         Enddo

         If (Triplet) Then
         Write(6,*)
         Write(6,"(2x,a)") "The Triplet Excited States:"
         Write(6,*)
         Do Irrep = 1, Nirrep
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "T") Then
                  Eed   = Eigval(Iroot,Irrep) * Factev
                  Eet   = Eigval_t(Iroot,Irrep) * Factev
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep)

                  Write(6,99) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Eed,Eet,F,Et
               Endif
            Enddo
          Enddo
         Endif

  99   Format(1x,i3,1x,i3,"[",i1,"]",4x,i3,"[",i1,"]",4x,F12.4,1x,F12.4,
     &        3x,E9.4,1x,F15.8)
      Else

      Write(6,*) 
      Write(6,10)
 10   Format(2x,'Sym',2x,'Origin',4x,'Destination',3x,'EE(eV)',5x,
     &       'EE(cm-1)',4x,'Osc. Str.',4x,'Total Energy')
      Write(6,*)' ----------------------------------------------------',
     &          '-------------------------'
      Write(6,*)
      Write(6,"(2x,a)") "The Singlet Excited States:" 
      Write(6,*)
         Do Irrep = 1, Nirrep 
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "S") Then
                  Ee    = Eigval(Iroot,Irrep) * Factev
                  Eecm  = Eigval(Iroot,Irrep) * Factcm
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep) 

                  Write(6,9) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Ee,Eecm,F,Et 
               Else 
                  Triplet = .True.
               Endif 
            Enddo
         Enddo

         If (Triplet) Then
         Write(6,*)
         Write(6,"(2x,a)") "The Triplet Excited States:"
         Write(6,*)
         Do Irrep = 1, Nirrep
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "T") Then
                  Ee    = Eigval(Iroot,Irrep) * Factev
                  Eecm  = Eigval(Iroot,Irrep) * Factcm
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep)

                  Write(6,9) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Ee,Eecm,F,Et
               Endif
            Enddo
          Enddo
          Endif 

      Endif 

      Write(6,*)
      Write(6,*)' ----------------------------------------------------',
     &          '-------------------------'
  9   Format(1x,i3,1x,i3,"[",i1,"]",4x,i3,"[",i1,"]",4x,F12.4,1x,F12.2,
     &       3x,E9.4,1x,F15.8)

      Write(6,*)
      Write(6,"(a,a)") " Comments! If there are empty rows all that ",
     &                 "means assigining singlet or triplet"
      Write(6,"(a,a)") " characted to a state has failed. Most likely",
     &                 " scenario is that they are not" 
      Write(6,"(a,a)") " pure singlets or triplets. Also, the zero",
     &                 " excitation energy indicates that the"
      Write(6,"(a)")   " root is not converged."

      RETURN
 1300    FORMAT(/)
      END
