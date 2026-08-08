













































































































































































































      Subroutine Pccd_frmblocks(grad,grad_oo,grad_vv,grad_vo,
     +                          grad_ov,work,Imemleft,Lenoo,
     +                          Lenvv,Lenvo,Nocc,Nvrt,Nbas)

      Implicit Double Precision (A-H, O-Z)

      Dimension Grad(Nbas,Nbas)
      Dimension Grad_oo(Lenoo),Grad_vv(Lenvv)
      Dimension Grad_vo(Lenvo),Grad_ov(Lenvo)
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
      Ovcount = 0
      Length  = Nbas*Nbas

      write(6,*)
      write(6,"(a,a)") "@-pccd_frmblocks: The incomming full matrix"
      call output(Grad,1,Nbas,1,Nbas,Nbas,Nbas,1)

C The occupied/occupied block

      Do Iocc = 1, Nocc
         Do Jocc = 1, Nocc
            Oocount = Oocount + 1
            Grad_oo(Oocount) = Grad(Jocc,Iocc)
         Enddo
      Enddo

C The virtual/virtual block

      Do Ivrt = 1, Nvrt
         Do Jvrt = 1, Nvrt
            Vvcount = Vvcount + 1
            Grad_vv(Vvcount) = Grad(Jvrt+Nocc,Ivrt+Nocc)
         Enddo
      Enddo

C The virtual-occupied (Vo)block

      Do Iocc = 1, Nocc
         Do Jvrt = 1, Nvrt
            Vocount = Vocount + 1
            Grad_vo(Vocount) = Grad(Jvrt+Nocc,Iocc)
         Enddo
      Enddo

C The occupied-virtual (ov) block

      Do Ivrt = 1, Nvrt
         Do Jocc = 1, Nocc
             Ovcount = Ovcount + 1
             Grad_ov(Ovcount) = Grad(Jocc,Ivrt+Nocc)
         Enddo
      Enddo


      Write(6,*)
      Write(6,"(a)") "The OO-MO block of the matrix"
      call output(Grad_oo,1,nocc,1,nocc,nocc,nocc,1)
      Write(6,"(a)") "The VV-MO block of the matrix"
      call output(Grad_vv,1,nvrt,1,nvrt,nvrt,nvrt,1)
      Write(6,"(a)") "The VO-MO block of the matrix"
      call output(Grad_vo,1,nvrt,1,nocc,nvrt,nocc,1)
      Write(6,"(a)") "The OV-MO block of the matrix"
      call output(Grad_ov,1,nocc,1,nvrt,nocc,nvrt,1)
      Return
      End
