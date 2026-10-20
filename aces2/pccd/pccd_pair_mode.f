














































































































































































































      Subroutine Pccd_pair_mode(Ipair)
C Ipair = 1: this xpccd run uses the pair-form solver (Pccd_pair_drv):
C pCCD, OO_METHOD=RFO, not Brueckner (PCCD_REF=2), and none of the
C debug files that need the ACES densities (pccd_frzamp.in,
C pccd_hess.in, pccd_ohess.in). The first word OFF of pccd_pair.in
C selects the CCD-shaped solvers (checks). RHF without dropped MOs (the
C caller).
      Implicit Double Precision(A-H,O-Z)
      Logical pCCD,CCD,LCCD,Exist
      Character*8 Word
      Common /CALC/pCCD,CCD,LCCD
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
      Ipair = 0
      Call Pccd_setmet()
      If (.Not. pCCD) Return
      If (Iflags2(178) .Ne. 6) Return
      If (Iflags2(179) .Eq. 2) Return
      Inquire(File="pccd_frzamp.in",Exist=Exist)
      If (Exist) Return
      Inquire(File="pccd_hess.in",Exist=Exist)
      If (Exist) Return
      Inquire(File="pccd_ohess.in",Exist=Exist)
      If (Exist) Return
      Inquire(File="pccd_pair.in",Exist=Exist)
      If (Exist) Then
         Word = " "
         Open(Unit=44,File="pccd_pair.in",Status="Old")
         Read(44,"(a)",End=5,Err=5) Word
 5       Close(44)
         If (Word(1:3) .Eq. "OFF") Return
      Endif
      Ipair = 1
      Return
      End
