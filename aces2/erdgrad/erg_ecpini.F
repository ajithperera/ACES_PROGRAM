      Subroutine Erg_ecpini(Charge)
C
C Set up ecplib for the ECP gradient, as Erd_ecp does for the ECP
C integrals: read the ECP definitions (Chrges, Ecpdef, Ecppar), load
C our shells (Erd_ecpset). The shells are given as Cartesian ones,
C since ecplib's gradient code works with Cartesian components; the
C density it gets is over the Cartesian ZMAT AOs (Erg_ecpgrd). On
C return Zatm holds the effective nuclear charges (core electrons of
C the ECPs removed), which the nuclear attraction and repulsion must
C use. Ecppar assigns its arguments: pass variables.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension Charge(Natm)
      Common /infoa / Natecp
      Character*2 Mtype(Mxatms)
      Character*80 Ecpnam(Mxatms), Cnick(Mxatms)
      Logical Seward

      Do I = 1, Natm
         Charge(I) = Zatm(I)
      Enddo
      Ierr = 0
      Call Chrges(Charge,Mtype,Natm,Ierr)
      If (Ierr .Ne. 0) Then
         Write(6,"(a)") " @Erg_ecpini: error in Chrges"
         Call Errex
      Endif
      Call Ecpdef(Natm,Mtype,Ecpnam,Cnick)
      Ishow  = 0
      Lusew  = 0
      Seward = .False.
      Call Ecppar(Ishow,Cnick,Charge,Natm,Seward,Lusew)
      Do I = 1, Natm
         Zatm(I) = Charge(I)
      Enddo
      Call Ecp_init

C ecplib's gradient code loops over the centers in /infoa/ (xvdint's
C Ecp_basis_init_4grd sets it); its centers are the ZMAT centers.

      Call Getrec(20,'JOBARC','NATOMS  ',1,Natecp)
      Call Erd_ecpset(Nshl,Ishatm,Lshl,Nprsh,Nctsh,Alpha,Iofalp,Ccoef,
     &                Iofcc,.False.)

      Return
      End
