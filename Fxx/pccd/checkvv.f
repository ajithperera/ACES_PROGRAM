      Subroutine checkvv(H,Nbas,Nvrt)

      Implicit Double Precision (A-H,O-Z)
      Dimension H(Nbas,Nbas)
      Dimension Ioffo(8)
      Dimension Ioffv(8)
      Integer An,Bn
   
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Ioffo(1) = 0
      Ioffv(1) = Nocco(1)

      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo

      Ispin = 1
      e = 0.0D0
      Do Irrep_b = 1, Nirrep
         Irrep_a = Dirprd(Irrep_b,1)
         Do b = 1, Vrt(Irrep_b,Ispin)
            Do a = 1, Vrt(Irrep_a,Ispin)
               An = A + Ioffv(Irrep_a)
               Bn = B + Ioffv(Irrep_b)
               e = e +  (H(an,bn)*H(an,bn))
            Enddo
          Enddo
      Enddo
      Print*, E
      Return
      End 
