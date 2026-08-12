













































































































































































































      Subroutine Job_control(Occnums_file) 
      Logical Ca,Fc,Bca,Bfc,Occnums_file,Move
      Logical Nwo
      Logical Scf, Pscf
      Character*80 Fname 
      
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      Data Izz /100/

      Ca  = .False.
      Fc  = .False.
      BCa = .False.
      Bfc = .False.
      Nwo = .False.
      Fcr = .False.

      If (Iflags(39) .EQ. 2 .OR. 
     +    Iflags(39) .EQ. 0)     Then
          Ca = .True. 
      ElseIf (Iflags(39) .EQ. 3) Then
          Fc = .True.
      ElseIf (Iflags(39) .EQ. 4) Then
          Bca = .True.
      ElseIf (Iflags(39) .EQ. 5) Then
          Bfc = .True.
      ElseIf (Iflags(39) .EQ. 6) Then
          Nwo = .True.
      Endif 

      Scf  = .False.
      Pscf = .False.

      If (iflags(2) .Eq. 0) Scf  = .True.
      If (iflags(2) .Eq. 1) PScf = .True.

      If (Nwo) Then
         If (Scf) Then
            Call Runit("xvmol")
            Call Runit("xvmol2ja")
            Call Runit("xvprops")
            Call Runit("xvscf")
            Call Runit("xget_acesinfo")
            Call Runit("python nwchem_mo.py nwchem.out dft")
            Call Runit("xget_acesmo")
            Call Runit("mv NEWMOS.work OLDMOS")
            Iflags(18)      = 1
            Iflags(16) = 0
            Iflags(45)      = 4
            Call A2putrec(20, 'JOBARC', 'IFLAGS  ', Izz, IFLAGS)
            Call Runit("xvscf")
            Call Runit("xprops")
            Call Runit("xa2mix")
         Else 
            Call Runit("xvmol")
            Call Runit("xvmol2ja")
            Call Runit("xvprops")
            Call Runit("xvscf")
            Call Runit("xvtran")
            Call Runit("xintprc")
            Call Runit("xvcc")
            Call Runit("xlambda")
            Call Runit("xdens")
            Call Runit("xprops")
            Call Runit("xa2mix")
         Endif
      Else 
         If (Iflags(22) .EQ. 1) Then
            Call Do_brueckner(BCa,Bfc,Occnums_file)

C After Brueckner is completed, we need to continue with the 
C Brueckner orbitals. Lets setup the flags for that. 
C Set GUESS=READ_SO_MOS and move NEWMOS to OLDMOS and remove the
C BRUECKNER flag. 

            Call Runit("mv NEWMOS OLDMOS")
            Iflags(45)      = 4
            Iflags(22)  = 0
            Iflags(16) = 0
            Call A2putrec(20, 'JOBARC', 'IFLAGS  ', Izz, IFLAGS)

            Call Gfname("OCCNUMS",Fname,Ilength)
            Inquire(File=Fname(1:7),Exist=Occnums_file)
            If (Occnums_file)  Call Runit("rm OCCNUMS")
            Move = (Occnums_file .and. Bfc .and. .not. Bca)
            If (Move) Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xvscf")
            Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xccsd_light")
         Else 
            Call Runit("xvmol")
            Call Runit("xvmol2ja")
            Move = (Occnums_file .and. Fc .and. .not. Ca) 
            If (Move) Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xvscf")
            Write(6,*)
            If (Occnums_file) Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xccsd_light")
         Endif 
      Endif 

      return
      End
