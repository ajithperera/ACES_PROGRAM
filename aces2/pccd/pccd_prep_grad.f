













































































































































































































      Subroutine Pccd_prep_grad(grad,grad_oo,grad_vv,grad_vo,
     +                          grad_ov,work,Imemleft,Lenoo,
     +                          Lenvv,Lenvo,Nocc,Nvrt,Nbfns)

      Implicit Double Precision (A-H, O-Z)

      Dimension Grad(Nbfns,Nbfns)
      Dimension Grad_oo(Nocc,Nocc),Grad_vv(Nvrt,Nvrt)
      Dimension Grad_vo(Nvrt,Nocc),Grad_ov(Nocc,Nvrt)
      Dimension Work(Imemleft)
      Dimension Ioffo(8),Ioffv(8)
      Integer Oocount,Vvcount,Vocount



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

      Irrepx  = 1
      Oocount = 0
      Vvcount = 0
      Vocount = 0
      oVcount = 0
      Length  = Nbfns*Nbfns


      Call Zero(Grad_oo,Nocc*Nocc)
      Call Zero(Grad_vv,Nvrt*Nvrt)
      Call Zero(Grad_vo,Nvrt*Nocc)
      Call Zero(Grad_ov,Nocc*Nvrt)

      Do I = 1, Nocc
         Do J = 1, Nocc
            Grad_oo(J,I) = Grad(J,I)
         Enddo
      Enddo 

      Ioff = 0
      Do I = Nocc+1, Nbfns 
         Ioff = Ioff + 1
         Joff = 0
         Do J = Nocc+1, Nbfns 
            Joff = Joff + 1
            Grad_vv(Joff,Ioff) = Grad(J,I)
         Enddo
      Enddo 

      Ioff = 0
      Do I = 1, Nocc
         Ioff = Ioff + 1
         Joff = 0
         Do J = Nocc+1, Nbfns
            Joff = Joff + 1
            Grad_vo(Joff,Ioff) = Grad(J,I)
         Enddo
      Enddo 

      Ioff = 0
      Do I = 1, Nocc 
         Ioff = Ioff + 1
         Joff = 0
         Do J = Nocc+1, Nbfns 
            Joff = Joff + 1
            Grad_ov(Ioff,Joff) = Grad(I,J)
         Enddo
      Enddo 
      Write(6,*)
      Write(6,"(a)") "The OO-MO gradient matrices"
      call output(Grad_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,"(a)") "The VV-MO gradient matrices"
      call output(Grad_vv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      Write(6,"(a)") "The VO-MO gradient matrices"
      call output(Grad_vo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
      Write(6,"(a)") "The OV-MO gradient matrices"
      call output(Grad_ov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)
      Return
      End
