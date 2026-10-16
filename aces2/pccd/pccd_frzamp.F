      Subroutine Pccd_frzamp(Work,Maxcor,Iuhf)
C
C Debug hook for checking orbital Hessians: if the file pccd_frzamp.in
C exists its first line selects
C   SAVE : write the converged T2(AB) (list 46) and L2(AB) (list 146)
C          to pccd_amps.dat
C   LOAD : overwrite lists 46 and 146 with the amplitudes in
C          pccd_amps.dat
C This is called after the left-hand equations are solved and before
C the densities are built, so with LOAD the orbital gradient is that of
C the pCCD Lagrangian at fixed (T,L) in the current orbital basis (the
C frozen-amplitude Hessian is obtained from it by finite differences).
C No-op when the file is absent. RHF only.
C
      Implicit Double Precision(A-H,O-Z)
      Character*8 Mode
      Logical Exist
      Dimension Work(Maxcor),Lists(2)
#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
      Data Lists /46,146/

      Inquire(File="pccd_frzamp.in",Exist=Exist)
      If (.Not. Exist) Return
      If (Iuhf .Ne. 0) Then
         Write(6,"(a)") " @-Pccd_frzamp: RHF only."
         Call Errex
      Endif
      Open(Unit=41,File="pccd_frzamp.in",Status="Old")
      Read(41,"(a)") Mode
      Close(41)
      If (Mode(1:4) .Eq. "SAVE") Then
         Open(Unit=42,File="pccd_amps.dat",Status="Unknown")
      Elseif (Mode(1:4) .Eq. "LOAD") Then
         Open(Unit=42,File="pccd_amps.dat",Status="Old")
      Else
         Write(6,"(2a)") " @-Pccd_frzamp: unknown mode ", Mode
         Call Errex
      Endif

      Do Il = 1, 2
         List = Lists(Il)
         Do Irrep = 1, Nirrep
            Ndis = Irpdpd(Irrep,Isytyp(2,List))
            Ldis = Irpdpd(Irrep,Isytyp(1,List))
            Nsize = Ndis*Ldis
            If (Nsize .Eq. 0) Go to 10
            If (Nsize .Gt. Maxcor) Call Insmem("Pccd_frzamp",Nsize,
     +                                         Maxcor)
            If (Mode(1:4) .Eq. "SAVE") Then
               Call Getlst(Work,1,Ndis,1,Irrep,List)
               Write(42,*) List, Irrep, Nsize
               Write(42,"(4(1x,E25.17))") (Work(I),I=1,Nsize)
            Else
               Read(42,*) Jlist, Jrrep, Jsize
               If (Jlist.Ne.List .Or. Jrrep.Ne.Irrep .Or.
     +             Jsize.Ne.Nsize) Then
                  Write(6,"(2a)") " @-Pccd_frzamp: pccd_amps.dat",
     +                            " mismatch"
                  Call Errex
               Endif
               Read(42,*) (Work(I),I=1,Nsize)
               Call Putlst(Work,1,Ndis,1,Irrep,List)
            Endif
 10         Continue
         Enddo
      Enddo
      Close(42)
      Write(6,"(3a)") " @-Pccd_frzamp: ", Mode(1:4),
     +                " of T2(AB)/L2(AB) done."
      Return
      End
