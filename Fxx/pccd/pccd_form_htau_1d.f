      Subroutine Pccd_form_htau_1d(Htau_pq,Htau_qp,Work,Maxcor,
     +                             Hpq,Dpq,Dhf,Doo,Dvv,Dov,Nocc,
     +                             Nvrt,Nbas,Iuhf,E)
      
      Implicit Double Precision(A-H,O-Z)
      Integer A,B
      Logical Sym_packed 

      Dimension Hpq(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
      Dimension Doo(Nocc,Nocc)
      Dimension Dvv(Nvrt,Nvrt)
      Dimension Dov(Nocc,Nvrt)
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



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

      Data Ione,One,Onem,DNull,Two,Half /1,1.0D0,-1.0D0,0.0D0,2.0D0,
     +                                   0.50D0/

C The expression coded here is H(p,r)D(r,q) - H(r,q)D(p,r). This is
C antisymmetric (from p+q)


      I000 = Ione
      I010 = I000 + Nocc*Nocc
      I020 = I010 + Nvrt*Nvrt
      I030 = I020 + Nvrt*Nocc
      Iend = I030 + Nvrt*Nocc
      Maxcor = Maxcor - Iend 

      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_1",Iend,
     +                                   Maxcor)

      Call Dgemm("T","T", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_qp,Nbas)

C Compute the One-Particle contribution to the energy. When normal-order
C this is one particle contribution to the correlation energy.


      Call Analyze_fock(Work(I000),Work(I010),Work(I030),Nfmi(1),
     +                  Nfea(1),Nt(1),Nonhf)

      Call Pccd_form_fd(Doo,Dvv,Dov,Work(I000),Work(I010),Work(I020),
     +                  Work(I030),Work(Iend),Maxcor,Iuhf,Nonhf)

       Sym_packed = .False. 
       Call Pccd_frmful(Htau_pq,Work(I000),Work(I010),Work(I020),
     +                  Work(I030),Work(Iend),Maxcor,Nocc,Nvrt,
     +                  Nbas,"Ov_like",Sym_packed)

      write(6,*)
      write(6,"(2a)") "F(p,q) =  F(p,r)D(r,q) for HF like Fock matrices"
     +                ," contributions added"
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)

      Return
      End
   
