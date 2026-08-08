      Subroutine Pccd_form_htau_1d_1(Htau_pq,Htau_qp,Work,
     +                               Maxcor,Hpq,Dpq,Fds,
     +                               Nocc,Nvrt,Nbas,Iuhf,E)
      
      Implicit Double Precision(A-H,O-Z)
      Integer A,B
      Logical Nonhf_ref
      Logical Nonhf_terms
      Logical Nonhf

      Dimension Hpq(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Fds(Nbas,Nbas)
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
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

      Common/Nhfref/Nonhf_ref 

      Data Ione,One,Onem,DNull,Two,Half /1,1.0D0,-1.0D0,0.0D0,2.0D0,
     +                                   0.50D0/

C The expression coded here is H(p,r)D(r,q) - H(r,q)D(p,r). This is
C antisymmetric (from p+q)

      Call Dgemm("N","N", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_qp,Nbas)
      Call Dscal(Nbas,Dnull,Htau_qp,Nbas+1)
      Call Dgemm("N","N", Nbas,Nbas,Nbas,One,Dpq,Nbas,Hpq,Nbas,
     +            One,Htau_qp,Nbas)

C Compute the One-Particle contribution to the energy. When normal-order
C this is one particle contribution to the correlation energy.

      write(6,*)
      write(6,"(a)") "F(q,p) =  F(r,q)D(p,r) contribution" 
      call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)

      E=dnull
      do i = 1, Nbas 
        E = E + Htau_qp(i,i)
      enddo 
      write(6,"(2a,(2x,F15.10))") "One particle contribution to the",
     +                            " energy :", E 

      I000 = Ione
      I010 = I000 + Nocc*Nocc
      I020 = I010 + Nvrt*Nvrt
      I030 = I020 + Nvrt*Nocc
      Iend = I030 + Nvrt*Nocc
      Maxcor = Maxcor - Iend

      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_1",
     +                                   Iend,Maxcor)

      Call Analyze_fock(Work(I000),Work(I010),Work(I030),Nfmi(1),
     +                  Nfea(1),Nt(1),Nonhf_terms)

     
      Nonhf_ref = Nonhf_terms 
      Nonhf     = Nonhf_ref

      I000 = Ione
      I010 = I000 + Nbas*Nbas
      Iend = I010 + Nbas*Nbas
      Maxcor = Maxcor - Iend

      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_1",
     +                                   Iend,Maxcor)
C Compute the contributions from the diagonals of the Fock matrix for Htau_pq.
C If there are non-diagonal elements then the nonhf get set and they are
C computed below and get added to Htau_pq. The Htau_qp is simply to keep 
C track of the trace the reference contributions are not added. 

      Call Dzero(Work(I000),Nbas*Nbas)
      Call Dzero(Fds,Nbas*Nbas)
      Call Dzero(Htau_pq,Nbas*Nbas)
      Call Dcopy(Nbas,Hpq,Nbas+1,Fds,Nbas+1)
      Call Dscal(Nbas*Nbas,Onem,Fds,1)
      Call Dscal(Nbas,Dnull,Hpq,Nbas+1)

      write(6,*)
      write(6,"(a)") "Fock matrix sans diagonals" 
      call output(Hpq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,"(a)") "Fock diagonals" 
      call output(Fds,1,Nbas,1,Nbas,Nbas,Nbas,1)

      Call Dgemm("N","N", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Work(I000),Nbas)

      write(6,*)
      write(6,"(2a)")"F(p,q)= F(p,r)D(r,q) p.ne.r contribution",
     +               " non-diagonal fock sans diagonal"
      Write(6,"(a)") "and not augemented with bare Fock elements"
      call output(Work(I000),1,Nbas,1,Nbas,Nbas,Nbas,1)
      Do A = Nocc+1, Nbas
         Do I = 1, Nocc
            Index = (A-1)*Nbas + I 
            Htau_pq(I,A) = - Work(I000-1+Index)
         Enddo
      Enddo 

      Do I = 1, Nocc 
         Do A = Nocc+1, Nbas
            Index = (I-1)*Nbas + A
            Htau_pq(A,I) =  Work(I000-1+Index) + Hpq(A,I)
         Enddo
      Enddo 

      Do I = 1, Nocc 
         Do J = 1, Nocc
            Index = (I-1)*Nbas + J 
            Htau_pq(I,J) =  Work(I000-1+Index) + 
     +                      Hpq(I,J)
         Enddo
      Enddo 

      Kndex = 0
      Do A = Nocc+1, Nbas
         Jindex = 0
         Kindex = Kindex + 1
         Do B = Nocc+1, Nbas
            Jindex = Jindex + 1
            Index = Nocc*Nbas + (Kindex-1)*Nbas + Nocc + Jindex 
            Htau_pq(A,B) =  Work(I000-1+Index) 
         Enddo
      Enddo 

      Call Pccd_proc_htau(Htau_pq,Work(Iend),Maxcor,Nbas,Nocc,Nvrt)

      write(6,*)
      write(6,"(2a)") "F(p,q)=F(p,r)D(r,q) contributions for ",
     +                " non-diagonal Fock matrices sans diagonals",
     +                "and bare Fock elements added"
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)

      Return
      End
   
