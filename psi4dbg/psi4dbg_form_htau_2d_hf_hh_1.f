










      Subroutine Psi4dbg_form_htau_2d_hf_hh_1(Hoo,Dpq,Dhf,Dcc,Work,
     +                                        Maxcor,Nocc,Nvrt,Nbas,
     +                                        List_v)

      Implicit Double Precision(A-H,O-Z)
      Integer D,C,Dn,Cn
      Integer B,A,Bn,An

      Dimension Work(Maxcor)
      Dimension Hoo(Nocc*Nocc)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
      Dimension Dcc(Nbas,Nbas)
      Dimension Ioffo(8),Ioffv(8)

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

      Data One,Onem,Dnull,Ione,Two,Quart,Half/1.0D0,-1.0D0,0.0D0,
     +                                        1,2.0D0,0.250D0,
     +                                        0.50D0/
      Ioffo(1) = 0
      Ioffv(1) = Nocco(1)
      Irrepx   = Ione
      Ispin    = Ione

      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo

C Note that Dpq at entry is Alpha+Beta  

      Do Irrep = 1, Nirrep
         N = Pop(Irrep,1)
         Do I = 1, N
            Dcc(Ioffo(Irrep)+I,Ioffo(Irrep)+I) =  
     +         (Dpq(Ioffo(Irrep)+I,Ioffo(Irrep)+I) - Two)
         Enddo
      Enddo


      call checksum("Htau_oo :",Hoo,Nocc*Nocc)
      Return
      End

 

