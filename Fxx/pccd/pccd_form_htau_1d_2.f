      Subroutine Pccd_form_htau_1d_2(Htau_pq,Htau_qp,Doo,Dvv,Dvo,
     +                               Fds,Work,Maxcor,Nocc,Nvrt,
     +                               Nbas,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)
      Integer A,B
      Logical Sym_packed,Nonhf_ref

      Dimension Hpq(Nbas,Nbas)
      Dimension Doo(Nocc*Nocc)
      Dimension Dvv(Nvrt*Nvrt)
      Dimension Dvo(Nvrt*Nocc)
      Dimension Fds(Nbas*Nbas)
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)

      Common/Nhfref/Nonhf_ref

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

      Data Ione,One,Onem,DNull,Two,Half /1,1.0D0,-1.0D0,0.0D0,2.0D0,
     +                                   0.50D0/

C Extract the occupied and virtual diagonal to two seperate blocks.

      I000 = Ione
      I010 = I000 + Nocc*Nocc
      I020 = I010 + Nvrt*Nvrt
      I030 = I020 + Nvrt*Nocc
      I040 = I030 + Nvrt*Nocc
      I050 = I040 + Nocc*Nocc
      I060 = I050 + Nvrt*Nvrt
      Iend = I060 + Nbas*Nbas
      If (Iend .Gt. Maxcor) Call Insmem("pccd_form_htau_1d_2",Iend,
     +                                   Maxcor)

      Call Dzero(Work(I000),Nocc*Nocc)
      Call Dzero(Work(I010),Nvrt*Nvrt)
      Call Dcopy(Nocc,Fds,Nbas+1,Work(I000),Nocc+1)
      Call Dcopy(Nvrt,Fds(Nocc*Nbas+Nocc+1),Nbas+1,Work(I010),Nvrt+1)

      write(6,*)
      write(6,"(a)") "OO Fock diagonals"
      call output(Work(I000),1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,*)
      write(6,"(a)") "VV Fock diagonals"
      call output(Work(I010),1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      Call Dgemm("T","N",Nocc,Nvrt,Nvrt,One,Dvo,Nvrt,Work(I010),Nvrt,
     +            Dnull,Work(I020),Nocc)

      Call Dgemm("N","N",Nvrt,Nocc,Nocc,One,Dvo,Nvrt,Work(I000),Nocc,
     +            Dnull,Work(I030),Nvrt)

      Write(6,*)
      Write(6,"(a)") "G2(vo)=F(vv)*Dvo contribution"
      call output(Work(I020),1,Nocc,1,Nvrt,Nocc,Nvrt,1)
      Write(6,*)
      Write(6,"(a)") "G2(ov)=F(oo)*D(vo)(t) contribution"
      call output(Work(I030),1,Nvrt,1,Nocc,Nvrt,Nocc,1)

      Call Pccd_proc_ovvo(Work(I020),Work(I030),Work(Iend),Maxcor,
     +                    Nocc,Nvrt)
      
      Call Dgemm("N","N",Nocc,Nocc,Nocc,One,Work(I000),Nocc,Doo,Nocc,
     +            Dnull,Work(I040),Nocc)
      Call Daxpy(Nocc,One,Work(I000),Nocc+1,Work(I040),Nocc+1)

      Call Dgemm("N","N",Nvrt,Nvrt,Nvrt,One,Work(I010),Nvrt,Dvv,Nvrt,
     +            Dnull,Work(I050),Nvrt)
      Write(6,*)
      Write(6,"(a)") "G2(oo)=Dpo(t)*F(oo) contribution"
      call output(Work(I040),1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,*)
      Write(6,"(a)") "G2(vv)=Dvv*F(vv) contribution"
      call output(Work(I050),1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      Sym_packed = .False.
      Call Pccd_frmful(Work(I060),Work(I040),Work(I050),Work(I030),
     +                 Work(I020),Work(Iend),Maxcor,Nocc,Nvrt,Nbas,
     +                 "Ov_like",Sym_paked)
      write(6,*)
      write(6,"(2a)") "F(p,q) =  F(p,r)D(r,q) contributions from the", 
     +                " diagonal of the Fock matrix"
      call output(Work(I060),1,Nbas,1,Nbas,Nbas,Nbas,1)

      Call Dscal(Nbas*Nbas,Onem,Htau_pq,1)
      Call Daxpy(Nbas*Nbas,One,Work(I060),1,Htau_pq,1)

      Do A = Nocc+1, Nbas
         Do I = 1, Nocc
            Index = (A-1)*Nbas + I
            Htau_pq(I,A) = -Htau_pq(I,A)
         Enddo
      Enddo

      Do I = 1, Nocc
         Do A = Nocc+1, Nbas
            Index = (I-1)*Nbas + A
            Htau_pq(A,I) = -Htau_pq(A,I)
         Enddo
      Enddo


      write(6,*)
      write(6,"(2a)")"F(p,q) =  F(p,r)D(r,q) contributions from the", 
     +               " digonal of the Fock matrix added"
      write(6,"(a)") " to the contributions from the non diagonal Fock."
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)

      Return
      End
   
