













































































































































































































      Subroutine do_brueckner(Xpccd)

      Implicit Double Precision(A-H,O-Z)
      Logical Converged,Bca,Bfc,Occnums_file,Move
      Character*80 Fname 
      Integer Icycle 
      Logical porbrot_in
      Logical Xpccd 
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



C
      Write(6,"(2a)"), " @-do_brueckner Entering pCCD Brueckner/or",
     +                 " orbital optimization block"
      
      Converged = .FALSE.
      Max_cycle = Iflags2(182)
      No_diis   = Iflags2(183)

      Inquire(file="porbrot.info",exist=porbrot_in)
      Icycle = 1
      Do while (.NOT. Converged)
          If (Icycle .EQ. 1) Then
             Call Runit("xvmol")
             Call Runit("xvmol2ja")
             Call Runit("xvscf")
             If (porbrot_in) Then
                Write(6,"(3a)")" Warning! - porbrot_in file is",
     +                         " present and specified extrenal",
     +                         " orbital rotations are performed!"
                Call Runit("xoprots")
             Endif 
          Else
             Call Runit("xvscf")
          Endif

CSSS          Call Aces_ja_init 
CSSS          Call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
CSSS          Call Aces_ja_fin

          Call Runit("xvtran")
          Call Runit("xintprc")
          If (Xpccd) Then
             Call Runit("xpccd")
CSSS             Call Runit("xpsi4dbg")
          Else
             Call Runit("xvcc")
          Endif
          If (Icycle .Eq. 1) Then
             Call A2getrec(20, 'JOBARC', 'SCFENEG ',Iintfp, Eref)
             Call A2getrec(20, 'JOBARC', 'TOTENERG',Iintfp, Ecor)
             Call A2putrec(20, 'JOBARC', 'SCFENEG0',Iintfp, Eref)
             Call A2putrec(20, 'JOBARC', 'TOTENEG0',Iintfp, Ecor) 
             Call Aces_ja_fin
          Endif 

          Call A2getrec(0, 'JOBARC', 'BRUKTEST', Ibtest, Ijunk)
          Call A2getrec(0, 'JOBARC', 'GRADTEST', Igtest, Ijunk)
          If (Ibtest .Gt. 0) Call A2getrec(20, 'JOBARC', 'BRUKTEST',
     +                                    1, Ibtest)
          If (Igtest .Gt. 0) Call A2getrec(20, 'JOBARC', 'GRADTEST',
     +                                    1, Igtest)
          Icycle = Icycle + 1
          If (Ibtest .Eq. 0) Then
          Write(6,"(a,i4,2a)") "  Starting ",Icycle, " Brueckner",
     +                         " pCCD iteration."
          Elseif (Igtest .Eq. 0) Then
          Write(6,"(a,i4,2a)") "  Starting ",Icycle," orbital ",
     +                         "optimization pCCD iteration."
          Endif
C--------- Eliminate after debugging is completed-----------
          Call Aces_ja_init 
          Call putrec(20,'JOBARC','ORBOPITR',1,Icycle)
          Call Aces_ja_fin
C----------------------------------------------------------

          IF (Ibtest .EQ. 1 .or. Igtest .EQ. 1)  Then
             Converged = .TRUE.
          Else if (Icycle .Gt. Max_cycle) Then
          If (Ibtest .Eq. 0) Then
             If (Max_cycle .Eq. 1) Return
             Write(6,"(a,a,i3,a)") "  The maximum allowed Brueckner",
     +                             " pCCD iterations ", Max_cycle, 
     +                             " has reached and no convergence!."
            Write(6,*)
            Call Errex 
          Elseif (Igtest .Eq. 0) Then
             If (Max_cycle .Eq. 1) Return
             Write(6,"(a,a,i3,a)") "  The maximum allowed orbital",
     +                             " optimization pCCD iterations ", 
     +                                Max_cycle, 
     +                             " has reached and no convergence!" 
            Write(6,*)
            Call Errex 
          Endif 
          Endif 
      Enddo
      
      Return
      End
