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
C The ring-ordered copies of T2 and L2 (lists 34-39, 134-139, those that
C exist) are saved and loaded too: the density code reads them (via
C Pccd_h4all), and the solvers leave them holding the amplitudes of the
C current orbitals. Freezing only 46 and 146 changed Gamma(123/125) by
C 1.8e-4 at a 0.02 rotation and made the frozen FD Hessian asymmetric
C by ~1e-3 (H2O STO-3G). After LOAD the one-particle intermediates G1
C (list 192) and G2 (191), which the Lambda driver forms and the density
C code reads (Pccd_g1inx2/g2inx2), are formed again from the loaded T2
C and L2.
C
      Implicit Double Precision(A-H,O-Z)
      Character*8 Mode
      Logical Exist
      Dimension Work(Maxcor),Lists(14)
#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
      Data Lists /46,146,34,35,36,37,38,39,134,135,136,137,138,139/
      Data One /1.0D0/

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

      Do Il = 1, 14
         List = Lists(Il)
         If (Isytyp(1,List) .Eq. 0 .Or. Isytyp(2,List) .Eq. 0) Go to 20
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
 20      Continue
      Enddo
      Close(42)
      If (Mode(1:4) .Eq. "LOAD") Then
         Call Pccd_formg1(Work,Maxcor,Iuhf,One)
         Call Pccd_formg2(Work,Maxcor,Iuhf,One)
      Endif
      Write(6,"(3a)") " @-Pccd_frzamp: ", Mode(1:4),
     +                " of T2(AB)/L2(AB) and ring copies done."
      Return
      End
