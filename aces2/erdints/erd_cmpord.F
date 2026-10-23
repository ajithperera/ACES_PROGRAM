      Subroutine Erd_cmpord(Icmpao,Ncrt,Iz2a,Natoms,Coord,Scr,Imap,
     &                      Ipop,Imem,Iaocen,Iaoang,Iaocn0,Iaoan0,
     &                      Nbsorb,Naoorb)
C
C Computational (VMOL) order of the Cartesian AOs, and the JOBARC
C records describing it. The rest of ACES II assumes this order
C (vscf/setup.F, for instance, rebuilds it from NBASATOM): orbit by
C orbit (joda's COMPPOPV/COMPMEMB), centers of an orbit in REMAP order
C (MAP2ZMAT, as vmol2ja writes it), and within an orbit
C
C    comp = (orbit offset) + (v-1)*npop + member
C
C where v runs over the functions of one atom in VMOL order: angular
C momentum, then Cartesian component, then all contractions of that
C angular momentum. The ZMAT AO order is the ERD Cartesian order
C (atoms in ZMAT order, shells as in MOL_ERDOED, contraction major).
C
C Output: Icmpao(i) is the computational index of ZMAT AO i.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "machsp.com"

      Integer Icmpao(Ncrt), Iz2a(Natoms), Imap(Natoms), Ipop(*),
     &        Imem(Natoms), Iaocen(Ncrt), Iaoang(Ncrt),
     &        Iaocn0(Ncrt), Iaoan0(Ncrt), Nbsorb(*), Naoorb(*)
      Dimension Coord(3,Natoms), Scr(*)
      Integer Ia2z(Mxatms), Ifsh(Mxatms+1)
      Character*4 Cmppgp

C ZMAT centers (dummies included) <-> MOL_ERDOED atoms (real only).

      Call Getrec(20,'JOBARC','ATOMCHRG',Natoms,Imap)
      Nreal = 0
      Do Iz = 1, Natoms
         Iz2a(Iz) = 0
         If (Imap(Iz) .Ne. 0) Then
            Nreal = Nreal + 1
            Iz2a(Iz) = Nreal
            If (Nreal .Le. Mxatms) Ia2z(Nreal) = Iz
         Endif
      Enddo
      If (Nreal .Ne. Natm) Then
         Write(6,"(a,2i5)") " @Erd_cmpord: atom count mismatch ",
     &                      Nreal, Natm
         Call Errex
      Endif
      Call Getrec(20,'JOBARC','COORD   ',3*Natoms*Iintfp,Coord)
      Do Ia = 1, Natm
         Dist = 0.0D0
         Do K = 1, 3
            Dist = Dist + Dabs(Coord(K,Ia2z(Ia)) - Xyzatm(K,Ia))
         Enddo
         If (Dist .Gt. 1.0D-6) Then
            Write(6,"(2a,i4)") " @Erd_cmpord: COORD and MOL_ERDOED ",
     &                        "disagree for atom ", Ia
            Call Errex
         Endif
      Enddo

C Computational centers, exactly as vmol2ja's driver.F does it.

      Call Getcrec(-1,'JOBARC','COMPPTGP',4,Cmppgp)
      Call Getrec(-1,'JOBARC','COMPORDR',1,Iordgp)
      Call Getrec(-1,'JOBARC','COMPNORB',1,Norbit)
      Call Getrec(-1,'JOBARC','COMPPOPV',Norbit,Ipop)
      Call Getrec(-1,'JOBARC','COMPMEMB',Natoms,Imem)
      Call Remap(Cmppgp,Iordgp,Natoms,Norbit,Coord,Scr,Imem,Ipop,
     &           Imap)
      Call Putrec(20,'JOBARC','MAP2ZMAT',Natoms,Imap)

C First shell of every atom.

      Do Ia = 1, Natm+1
         Ifsh(Ia) = Nshl + 1
      Enddo
      Do Is = Nshl, 1, -1
         Ifsh(Ishatm(Is)) = Is
      Enddo
      Ifsh(Natm+1) = Nshl + 1
      Do Ia = Natm, 1, -1
         If (Ifsh(Ia) .Gt. Ifsh(Ia+1)) Ifsh(Ia) = Ifsh(Ia+1)
      Enddo

      Do I = 1, Ncrt
         Icmpao(I) = 0
      Enddo
      Icomp = 0
      Icen  = 0
      Do Io = 1, Norbit
         Npop  = Ipop(Io)
         Norbm0 = Icomp
         Nsoo  = 0
         Do M = 1, Npop
            Icen = Icen + 1
            Iz   = Imap(Icen)
            If (Iz .Le. 0) Then
               Write(6,"(a,i4)") " @Erd_cmpord: REMAP failed, center",
     &                           Icen
               Call Errex
            Endif
            Ia = Iz2a(Iz)
            If (Ia .Eq. 0) Then
               Write(6,"(a)") " @Erd_cmpord: dummy atom in an orbit"
               Call Errex
            Endif

C VMOL order on this atom: l, then component, then contraction.

            Iv = 0
            Lmx = -1
            Do Is = Ifsh(Ia), Ifsh(Ia+1)-1
               Lmx = Max(Lmx,Lshl(Is))
            Enddo
            Do L = 0, Lmx
               Ncmp = (L+1)*(L+2)/2
               Do Ic = 1, Ncmp
                  Do Is = Ifsh(Ia), Ifsh(Ia+1)-1
                     If (Lshl(Is) .Eq. L) Then
                        Do It = 1, Nctsh(Is)
                           Iv = Iv + 1
                           Iao = Iofcrt(Is) + (It-1)*Ncmp + Ic
                           Ic0 = Norbm0 + (Iv-1)*Npop + M
                           Icmpao(Iao) = Ic0
                           Iaocen(Ic0) = Icen
                           Iaoang(Ic0) = L
                           Iaocn0(Iao) = Iz
                           Iaoan0(Iao) = L
                        Enddo
                     Endif
                  Enddo
               Enddo
            Enddo
            Do Is = Ifsh(Ia), Ifsh(Ia+1)-1
               Nsoo = Nsoo + Nfnsh(Is)
            Enddo
            If (M .Eq. 1) Nv1 = Iv
            If (Iv .Ne. Nv1) Then
               Write(6,"(a)") " @Erd_cmpord: orbit members differ"
               Call Errex
            Endif
         Enddo
         Icomp = Icomp + Npop*Nv1
         Naoorb(Io) = Npop*Nv1
         Nbsorb(Io) = Nsoo
      Enddo

      If (Icomp .Ne. Ncrt) Then
         Write(6,"(a,2i6)") " @Erd_cmpord: AO count mismatch ", Icomp,
     &                      Ncrt
         Call Errex
      Endif
      Do I = 1, Ncrt
         If (Icmpao(I) .Eq. 0) Then
            Write(6,"(a,i6)") " @Erd_cmpord: AO not placed ", I
            Call Errex
         Endif
      Enddo

      Call Putrec(20,'JOBARC','CENTERBF',Ncrt,Iaocen)
      Call Putrec(20,'JOBARC','ANGMOMBF',Ncrt,Iaoang)
      Call Putrec(20,'JOBARC','CNTERBF0',Ncrt,Iaocn0)
      Call Putrec(20,'JOBARC','ANMOMBF0',Ncrt,Iaoan0)
      Call Putrec(20,'JOBARC','NBASATOM',Norbit,Nbsorb)
      Call Putrec(20,'JOBARC','NAOBFORB',Norbit,Naoorb)

      Return
      End
