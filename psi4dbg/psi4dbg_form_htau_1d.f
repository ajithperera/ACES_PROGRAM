










      Subroutine Psi4dbg_form_htau_1d(Htau_pq,Htau_qp,Work,Maxcor,
     +                                Hpq,Dpq,Nocc,Nvrt,Nbas,E)
      
      Implicit Double Precision(A-H,O-Z)

      Dimension Hpq(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)

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

C The expression coded here is H(p,r)D(r,q) - H(r,q)D(p,r). This is
C antisymmetric (from p+q)


      I000 = Ione
      I010 = I000 + Nbas*Nbas
      Iend = I010 + Nbas*Nbas

      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d",Iend,Maxcor)

      Call Dgemm("N","N", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_pq,Nbas)

      Call Dgemm("T","T", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_qp,Nbas)

C Compute the One-Particle contribution to the energy. When normal-order
C this is one particle contribution to the correlation energy.

      write(6,*)
      write(6,"(a)") "F(p,q) =  F(p,r)D(r,q) contribution" 
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)
      write(6,"(a)") "F(q,p) =  F(r,q)D(p,r) contribution" 
      call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)

      E=dnull
      do i = 1, Nbas 
        E = E + Htau_pq(i,i)
      enddo 
      write(6,"(2a,(2x,F15.10))") "One particle contribution to the",
     +                            " MBPT(2) correlation energy :", 
     +                              E


      Return
      End
   
