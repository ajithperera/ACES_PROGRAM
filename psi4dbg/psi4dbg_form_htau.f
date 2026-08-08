










      Subroutine Psi4dbg_form_htau(Htau_qp,Htau_pq,Hoo_qp,Hoo_pq,
     +                             Hvv_qp,Hvv_pq,Hov,Hvo,Work,
     +                             Maxcor,Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Hoo_pq(Nocc*Nocc)
      Dimension Hoo_qp(Nocc*Nocc)
      Dimension Hvv_pq(Nvrt*Nvrt)
      Dimension Hvv_qp(Nvrt*Nvrt)
      Dimension Hvo(Nvrt*Nocc)
      Dimension Hov(Nocc*Nvrt)
      Dimension Work(Maxcor)

      Data Ione,Onem,One,Dnull,Half,Two,Four/1,-1.0D0,1.0D0,
     +                                       0.0D0,0.50D0,
     +                                       2.0D0,4.0D0/

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end


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




      Nmax = Max(Nvrt,Nocc)

      I000 = Ione
      I010 = I000 + Nmax*Nmax
      Iend = I010 + Nbas*Nbas
      If (Iend.Gt.Maxcor) Call Insmem("Pccd_form_htau",Iend,Maxcor)

      Write(6,"(a)") " The Htau_qp (from one-electron Hamiltonian)"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Write(6,"(a)") " The Htau_pq (from one-electron Hamiltonian)"
      Call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Call Dscal(Nocc*Nocc,Half,Hoo_qp,1)
      Call Dscal(Nvrt*Nvrt,Half,Hvv_qp,1)

      Write(6,"(a)") " The Hoo_qp"
      call output(Hoo_qp,1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,"(a)") " The Hvv_qp"
      call output(Hvv_qp,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)

      Maxcor = Maxcor - Iend 
      Call Pccd_frmful(Hoo_qp,Hvv_qp,Work(I010),Work(Iend),Maxcor,
     +                 Nbas,0)
      Call Daxpy(Nbas*Nbas,One,Work(I010),1,Htau_pq,1)

      Write(6,"(a,a)") " The Htau_qp (before OV/VO blocks added",
     +                 " and not anti-symmetrized)"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Write(6,"(a)") "Checking VO and OV blocks"
      Ioff = Ione
      Do Irrep = 1, Nirrep
         N = Pop(Irrep,1)
         M = Vrt(Irrep,1)
         Call output(Hvo(Ioff),1,M,1,N,M,N,1)
         Call output(Hov(Ioff),1,N,1,M,N,M,1)
         Ioff = Ioff + N*M
      Enddo

      Call Pccd_frmful_ov(Hov,Hvo,Htau_pq,Work(Iend),Maxcor,Nbas,
     +                    "ADD ",0)

      Write(6,"(a,a)") " The Htau_qp (after OV/VO blocks added",
     +                 " and not anti-symmetrized)"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)

C Antisymmetrize to generate the orbital rotation gradients matrix. The scalling
C by minus two accomplish two things. The minus is needed to get the correct 
C sign for the VO and OV blocks (VO block must be positive). The two gives the
C orbital rotaion gradient of the Alpha block (and match with half of the value
C obtained numerically). Note that proceeding steps need to be tailored to the
C fact that we are working with alpha block only. 

      Call Transp(Htau_pq,Work(I010),Nbas,Nbas)
      Call Daxpy(Nbas*Nbas,Onem,Work(I010),1,Htau_pq,1)
      Call Dscal(Nbas*Nbas,-Two,Htau_pq,1)

      Write(6,"(a)") " The orbital rotation gradient"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Write(6,"(a)") " Recheck The Hoo_qp"
      call output(Hoo_qp,1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,"(a)") " Recheck The Hvv_qp"
      call output(Hvv_qp,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)

      Call Pccd_frmful(Hoo_qp,Hvv_qp,Work(I010),Work(Iend),Maxcor,
     +                 Nbas,0)

C Scale the diagonals by half (correlation contribution only).

CSSS      Call Dscal(Nbas*Nbas,Half,Work(I010),Nbas+1)
      Call Daxpy(Nbas*Nbas,One,Work(I010),1,Htau_qp,1)

      Call Pccd_frmful_ov(Hov,Hvo,Htau_qp,Work(Iend),Maxcor,Nbas,
     +                    "ADD ",0)
CSSS      Call Dscal(Nbas*Nbas,Two,Htau_qp,1)

      Call Getrec(20, 'JOBARC', 'NUCREP', Iintfp, Zrepl)

      Write(6,"(a)") " The Htau_qp"
      Call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
      e=dnull
      do i = 1, Nbas 
         e = e + Htau_qp(i,i)*Two
      enddo
      write(6,*)
      write(6,"(a,1x,F15.10)") "The electronic energy :",e
    
      write(6,*)
      e = e 
      write(6,"(a,1x,F15.10)") "The correlation energy:", e
      write(6,*)
      Return
      End


