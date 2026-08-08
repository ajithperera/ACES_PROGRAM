










      Subroutine Pccd_check_1d(Label,Hpq,N,Type,S)

      Implicit Double Precision(A-H,O-Z)
      Character*9 Label 
      Character*2 Type
      Character*1 S
      Dimension Hpq(N*N)

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

      Data Dnull,Ione,Inull /0.0D0,1,0/

      E   = Dnull
      Ioff= Ione
      do Irrep = 1, Nirrep
         If (Type .Eq. "OO") M = Pop(Irrep,1)
         If (Type .Eq. "VV") M = Vrt(Irrep,1)
         Call Pccd_sum(Hpq(Ioff),M,E,S)
         Ioff = Ioff + M*M
      enddo

      write(6,"(a,2(2x,F15.10))") Label, e,e**2 
 
      return
      end 
