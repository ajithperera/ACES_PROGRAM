










      Subroutine Get_mbpt2_nos(Doo,Dvv,Coo,Cvv,Work,Maxcor,Iuhf)
     &                         
      Implicit Double Precision(A-H,O-Z)

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      Dimension Doo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Dvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Work(Maxcor)

      Data Ione,Done,Dnull,Inull,Mone /1,1.0D0,0.0D0,0,-1/
     
      Do Ispin = 1, (Iuhf+1)
         Ioff = (Ispin-1)*Nfmi(1) + Ione
         Joff = (Ispin-1)*Nfea(1) + Ione

         Do Irrep = 1, Nirrep
            Nd1 = pop(Irrep,Ispin)
            Nd2 = vrt(Irrep,Ispin)

            Call Eig(Doo(Ioff),Coo(Ioff),Ione,Nd1,Mone) 
            Call Eig(Dvv(Joff),Cvv(Joff),Ione,Nd2,Mone)

      Write(6,*) 
      If (ispin .eq. 1) then
         Write(6,"(a)") " The occ-occ mbpt(2) alpha natural orbitals"
      Else
         Write(6,"(a)") " The occ-occ mbpt(2) beta natural orbitals"
      Endif 
      call checksum("COO  :",Coo(ioff),nd1*Nd1)
CSSS      Call output(Coo(ioff),1,nd1,1,nd1,nd1,nd1,1)
      If (ispin .eq. 1) then
      Write(6,"(a)") " The vrt-vrt mbpt(2) alpha natural orbitals"
      Else
      Write(6,"(a)") " The vrt-vrt mbpt(2) beta natural orbitals"
      Endif 
CSSS      Call output(Cvv(joff),1,nd2,1,nd2,nd2,nd2,1)
      call checksum("Cvv  :",Cvv(joff),nd2*Nd2)
            Ioff = Ioff + Nd1*Nd1
            Joff = Joff + Nd2*Nd2
         Enddo 
      Enddo


      Return
      End


