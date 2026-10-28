      Subroutine Erg_ecpgrd(G,Dso,Nsb,Ncrt,Natoms,Iz2a,Coord,Z,P,Scr,
     &                      Gz)
C
C Add the ECP part of the gradient, from ecplib's gradient driver as
C xvdint's Ecp_dgrad uses it. The density it wants is over the
C Cartesian AOs in the ZMAT order, P = Z Dso Z(T) with Z = CMP2ZMAT
C and Dso the total SO density (Erg_dens). Its centers are the ZMAT
C centers (COORD, dummies included); a center carries the shells of
C the atom Iz2a points to (Erg_ecpini loaded the shells by atom), and
C the gradient comes back by ZMAT center. Ecp_dgrad doubles it.
C Cint is ecplib's scratch (Maxmem = 50000 in ecp.par).
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "machsp.com"

      Dimension G(3,Natm), Dso(Nsb,Nsb), Coord(3,Natoms), Z(Ncrt,Nsb),
     &          P(Ncrt,Ncrt), Scr(Ncrt,Nsb), Gz(3,Natoms), Cint(50000)
      Integer Iz2a(Natoms), Igenby(Mxatms)
      Logical Exist

      Call Getrec(20,'JOBARC','CMP2ZMAT',Ncrt*Nsb*Iintfp,Z)
      Call Dgemm('N','N',Ncrt,Nsb,Nsb,1.0D0,Z,Ncrt,Dso,Nsb,0.0D0,
     &           Scr,Ncrt)
      Call Dgemm('N','T',Ncrt,Ncrt,Nsb,1.0D0,Scr,Ncrt,Z,Ncrt,0.0D0,
     &           P,Ncrt)
      Call Getrec(20,'JOBARC','COORD   ',3*Natoms*Iintfp,Coord)
      Do I = 1, Natoms
         Igenby(I) = Iz2a(I)
         If (Igenby(I) .Eq. 0) Igenby(I) = 999
      Enddo

C The driver opens ECPGRD as a new file.

      Inquire(File='ECPGRD',Exist=Exist)
      If (Exist) Then
         Open(Unit=99,File='ECPGRD',Status='Old')
         Close(Unit=99,Status='Delete')
      Endif
      Call Dzero(Gz,3*Natoms)
      Call Dzero(Cint,50000)
      Call Ecp_grdint_driver(Natoms,Igenby,Coord,Cint,P,Ncrt,Gz,
     &                       .True.)

      Do I = 1, Natoms
         Ia = Iz2a(I)
         If (Ia .Gt. 0) Then
            Do K = 1, 3
               G(K,Ia) = G(K,Ia) + 2.0D0*Gz(K,I)
            Enddo
         Endif
      Enddo

      Return
      End
