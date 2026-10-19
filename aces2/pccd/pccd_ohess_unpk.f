










      Subroutine Pccd_ohess_unpk(Buf,T,List,Is1,Is2,Is3,Is4,Nocc,Nvrt)
C
C Unpack the totally symmetric list List (RHF, rows pq, columns rs,
C each DPD pair ordered by the irrep of its second index, first index
C fastest) into the dense T(p,q,r,s) over the full index spaces
C Is = 1 (occupied, 1..Nocc) or 2 (virtual, 1..Nvrt), orbitals of each
C space numbered by irrep. Buf holds one irrep block.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Buf(*),T(*),N(4),Is(4),Ld(4),Iof(8,2),Nsp(8,2)
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      Is(1) = Is1
      Is(2) = Is2
      Is(3) = Is3
      Is(4) = Is4
      Iob = 0
      Ivb = 0
      Do Irr = 1, Nirrep
         Iof(Irr,1) = Iob
         Iof(Irr,2) = Ivb
         Nsp(Irr,1) = Pop(Irr,1)
         Nsp(Irr,2) = Vrt(Irr,1)
         Iob = Iob + Pop(Irr,1)
         Ivb = Ivb + Vrt(Irr,1)
      Enddo
      Do K = 1, 4
         N(K) = Nocc
         If (Is(K) .Eq. 2) N(K) = Nvrt
      Enddo
      Ld(1) = 1
      Ld(2) = N(1)
      Ld(3) = N(1)*N(2)
      Ld(4) = N(1)*N(2)*N(3)
      Call Dzero(T,N(1)*N(2)*N(3)*N(4))

      Do Ix = 1, Nirrep
         Nrow = Irpdpd(Ix,Isytyp(1,List))
         Ncol = Irpdpd(Ix,Isytyp(2,List))
         If (Nrow*Ncol .Eq. 0) Go to 10
         Call Getlst(Buf,1,Ncol,1,Ix,List)
         Icol = 0
         Do I4r = 1, Nirrep
         I3r = Dirprd(I4r,Ix)
         Do I4 = 1, Nsp(I4r,Is(4))
         Do I3 = 1, Nsp(I3r,Is(3))
            Icol = Icol + 1
            K34 = (Iof(I3r,Is(3))+I3-1)*Ld(3)
     +          + (Iof(I4r,Is(4))+I4-1)*Ld(4)
            Irow = 0
            Do I2r = 1, Nirrep
            I1r = Dirprd(I2r,Ix)
            Do I2 = 1, Nsp(I2r,Is(2))
            Do I1 = 1, Nsp(I1r,Is(1))
               Irow = Irow + 1
               K = K34 + (Iof(I1r,Is(1))+I1)
     +           + (Iof(I2r,Is(2))+I2-1)*Ld(2)
               T(K) = Buf(Irow+(Icol-1)*Nrow)
            Enddo
            Enddo
            Enddo
         Enddo
         Enddo
         Enddo
 10      Continue
      Enddo
      Return
      End
