      Subroutine Erd_ecp(Natoms,Iz2a,Coord,Charge)
C
C ECP integrals with ecplib, without VMOL: read the ECP definitions
C (ECPNAM from joda, ECPDATA) as ecplib's Ecp_main does, load our
C shells into ecplib (Erd_ecpset), and let Ecp_int_driver write
C ECP1INTS (ZMAT ordered Cartesian AOs; vscf symmetry adapts it with
C CMP2ZMAT). On return Zatm holds the effective nuclear charges
C (core electrons of the ECPs removed), which the nuclear attraction
C and repulsion must use, and NMPROTON is written as VMOL does.
C
C Natoms, Iz2a(Natoms), Coord(3,Natoms): ZMAT centers (dummies
C included), their MOL_ERDOED atom (0 for dummies) and coordinates.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "flags.com"
#include "machsp.com"

      Integer Iz2a(Natoms)
      Dimension Coord(3,Natoms), Charge(Natm), Cint(50000)
      Character*2 Mtype(Mxatms)
      Character*80 Ecpnam(Mxatms), Cnick(Mxatms)
      Logical Nharm, Seward

      Nharm = (Iflags(62) .Eq. 1)
      Do I = 1, Natm
         Charge(I) = Zatm(I)
      Enddo
      Ierr = 0
      Call Chrges(Charge,Mtype,Natm,Ierr)
      If (Ierr .Ne. 0) Then
         Write(6,"(a)") " @Erd_ecp: error in Chrges"
         Call Errex
      Endif
      Call Ecpdef(Natm,Mtype,Ecpnam,Cnick)
C (Ecppar assigns to its print flag, so pass variables.)

      Ishow = 0
      Lusew = 0
      Seward = .False.
      Call Ecppar(Ishow,Cnick,Charge,Natm,Seward,Lusew)

      Nproton = 0
      Do I = 1, Natm
         Zatm(I) = Charge(I)
         Nproton = Nproton + Idint(Charge(I))
      Enddo
      Call Putrec(20,'JOBARC','NMPROTON',1,Nproton)

      Call Ecp_init
      Call Erd_ecpset(Nshl,Ishatm,Lshl,Nprsh,Nctsh,Alpha,Iofalp,Ccoef,
     &                Iofcc,Nharm)
      Call Getrec(20,'JOBARC','COORD   ',3*Natoms*Iintfp,Coord)
      Call Dzero(Cint,50000)
      Call Ecp_int_driver(Natm,Natoms,Iz2a,Coord,Cint,.False.)

      Return
      End
