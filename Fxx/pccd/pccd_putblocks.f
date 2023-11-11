



























































































































































































      Subroutine Pccd_Putblocks(Grad,Grad_oo,Grad_vv,Grad_vo,
     +                          Grad_ov,work,Imemleft,Lenoo,
     +                          Lenvv,Lenvo,Nocc,Nvrt,Nbas,
     +                          String,Sym_packed)

      Implicit Double Precision (A-H, O-Z)
      Character*7 String 
      Logical Symmetry,Sym_packed 

      Dimension Grad(Nbas,Nbas)
      Dimension Grad_oo(Lenoo),Grad_vv(Lenvv)
      Dimension Grad_vo(Lenvo),Grad_ov(Lenvo)
      Dimension Work(Imemleft)
      Dimension Ioffo(8),Ioffv(8)
      Integer Oocount,Vvcount,Vocount,Ovcount



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
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

      Common/Symm/Symmetry

      Data Ione /1/

      Irrepx  = 1
      Oocount = 0
      Vvcount = 0
      Vocount = 0
      Ovcount = 0
      
      Call Getrec(0,"JOBARC","SPNINDEX",Len,Ispin)
      If (Len .Gt .0) Then 
         Call Getrec(20,"JOBARC","SPNINDEX",Ione,Ispin)
      Else
         Ispin = Ione
      Endif 

      Ioff_oo = 1 + (Ispin-1) * Irpdpd(Irrepx,21)
      Ioff_vv = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
      Ioff_vo = 1 + (Ispin-1) * Irpdpd(Irrepx,9)

      Ndim_oo = Irpdpd(Irrepx,20+Ispin)
      Ndim_vv = Irpdpd(Irrepx,18+Ispin)
      Ndim_vo = Irpdpd(Irrepx,8+Ispin)
      Ndim_ov = Ndim_vo

C Calculate offsets occ and vrts.

      Ioffo(1) = 0
      Ioffv(1) = Nocco(Ispin)

      Do Irrep = 1, Nirrep-1
         Ioffo(Irrep+1) = Ioffo(Irrep) + Pop(Irrep,Ispin)
         Ioffv(Irrep+1) = Ioffv(Irrep) + Vrt(Irrep,Ispin)
      Enddo

      write(6,*)
      write(6,"(a,i2)") "@-pccd_putblocks, The Ispin:", Ispin
      write(6,"(a,1x,8i2)") "ioffo:",(Ioffo(i),i=1,Nirrep)
      write(6,"(a,1x,8i2)") "ioffv:",(Ioffv(i),i=1,Nirrep)

      Write(6,*)
      Write(6,"(a)") "Printing from Pccd_putblocks"
      Write(6,"(a)") "The OO-MO gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_oo(i),i=1,Lenoo)
      Write(6,"(a)") "The VV-MO  gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vv(i),i=1,Lenvv)
      Write(6,"(a)") "The VO-MO  gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vo(i),i=1,Lenvo)
      Write(6,"(a)") "The OV-MO  gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_ov(i),i=1,Lenvo)

      If (Sym_packed) Then
          Ioffoo = Ione
          Ioffvv = Ione
          Ioffvo = Ione

          Nocc2 = Nocc*Nocc
          Nvrt2 = Nvrt*Nvrt
          Nvovo = Nvrt*Nocc
          Call Dzero(Grad,Nbas*Nbas)

          I000 = Ione
          I010 = I000 + Nocc2 
          Iend = I010 + Nocc2 
          Call Putblk(Grad,Grad_oo(Ioffoo),"OO",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,20+Ispin),
     +                Nocc2,Ispin)
          I000 = Ione
          I010 = I000 + Nvrt2 
          Iend = I010 + Nvrt2
          Call Putblk(Grad,Grad_vv(Ioffvv),"VV",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,18+Ispin),
     +                Nvrt2,Ispin)
          I000 = Ione
          I010 = I000 + Nvovo
          Iend = I010 + Nvovo
          Call Putblk(Grad,Grad_vo(Ioffvo),"VO",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,8+Ispin),
     +                Nvovo,Ispin)

          If (String .Eq. "Special") Then

          Do Irrepr = 1, Nirrep
             Irrepl = Dirprd(Irrepr,Irrepx)
             Do Iocc = 1, Pop(Irrepr,Ispin)
                Do Ivrt = 1, Vrt(Irrepl,Ispin)
                  Ovcount =  Ovcount + 1 
                  Grad(iocc+ioffo(irrepr),ivrt+ioffv(irrepl)) = 
     +                                    Grad_ov(Ovcount)
                Enddo
             Enddo
          Enddo

          Else

          I000 = Ione
          I010 = I000 + Nvovo
          Iend = I010 + Nvovo
          Call Putblk(Grad,Grad_ov(Ioffvo),"OV",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,8+Ispin),
     +                Nvovo,Ispin)
          Endif 

          Ioffo  = Ioffo  + Irpdpd(Irrepx,20+Ispin)
          Ioffv  = Ioffv  + Irpdpd(Irrepx,18+Ispin)
          Ioffvo = Ioffvo + Irpdpd(Irrepx,8+Ispin)
             
      Else

C The occupied/occupied block

          Call Dzero(Grad,Nbas*Nbas)
          Do Iocc = 1, Nocc
             Do Jocc = 1, Nocc
                 Oocount = Oocount +1
                 Grad(Jocc,Iocc) = Grad_oo(Oocount)
             Enddo
          Enddo 

C The virtual/virtual block

          Do Ivrt = 1, Nvrt
             Do Jvrt = 1,Nvrt
                Vvcount = Vvcount +1
                Grad(Jvrt+Nocc,Ivrt+Nocc) = Grad_vv(Vvcount)
            Enddo
          Enddo 

C The virtual-occupied (Vo)block

          Do Iocc = 1, Nocc
             Do Jvrt = 1, Nvrt
                Vocount = Vocount + 1
                Grad(Jvrt+Nocc,Iocc) = Grad_vo(Vocount)
             Enddo
          Enddo

C The occupied-virtual (ov) block

          If (String .Eq. "Special") Then

             Do Jocc = 1, Nocc
                Do Ivrt = 1, Nvrt
                   Ovcount = Ovcount + 1
                   Grad(Jocc,Ivrt+Nocc) = Grad_ov(Ovcount)
                Enddo
             Enddo

          Else

             Do Ivrt = 1, Nvrt
                Do Jocc = 1, Nocc
                   Ovcount = Ovcount + 1
                   Grad(Jocc,Ivrt+Nocc) = Grad_ov(Ovcount)
                Enddo
             Enddo 

          Endif 
      ENdif 


      Return
      End
