













































































































































































































      Subroutine Tdee_prep_qdpole_ints(Qdp_ints,Qdp_ints_full,Scr1,
     +                              Scr2,Qdp_ints_oo,Qdp_ints_vv,
     +                              Qdp_ints_vo,Work,Imemleft,Lenoo,
     +                              Lenvv,Lenvo,Nbfns,Naobfns,Iuhf,
     +                              Ipert,Label,Irrepx)

      Implicit Double Precision (A-H, O-Z)

      Dimension Scr1(Naobfns*Nbfns),Scr2(Naobfns*Nbfns)
      Dimension Qdp_ints(Naobfns*(Naobfns+1)/2)
      Dimension Qdp_ints_full(Nbfns,Nbfns)
      Dimension Qdp_ints_oo(Lenoo),Qdp_ints_vv(Lenvv)
      Dimension Qdp_ints_vo(Lenvo)
      Dimension Work(Imemleft)
      Dimension Ioffo(8),Ioffv(8)
      Character*8 Label 
      Integer Oocount,Vvcount,Vocount

      Logical Vprop, Abacus 



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

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

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

      Vprop   = .True.
      Abacus  = .False.
      Length  = (Naobfns+1)*Naobfns/2
      Oocount = 0
      Vvcount = 0
      Vocount = 0

      Do Ispin = 1, Iuhf+1

         Ioff_oo = 1 + (Ispin-1) * Irpdpd(Irrepx,21)
         Ioff_vv = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
         Ioff_vo = 1 + (Ispin-1) * Irpdpd(Irrepx,9)

         Ndim_oo = Irpdpd(Irrepx,20+Ispin)
         Ndim_vv = Irpdpd(Irrepx,18+Ispin)
         Ndim_vo = Irpdpd(Irrepx,8+Ispin)

C Calculate offsets occ and vrts.

         Ioffo(1) = 0
         Ioffv(1) = Nocco(Ispin)
      
         Do Irrep = 1, Nirrep-1
            Ioffo(Irrep+1) = Ioffo(Irrep) + Pop(Irrep,Ispin)
            Ioffv(Irrep+1) = Ioffv(Irrep) + Vrt(Irrep,Ispin)
         Enddo

         If (Vprop) Then 
            Call Getrec(20,"JOBARC",Label,Length*Iintfp,
     +                 Qdp_ints)
            Call Expnd2(Qdp_ints,Qdp_ints_full,Naobfns)
            Call Ao2mo2(Qdp_ints_full,Qdp_ints_full,Scr1,Scr2,
     +                    Nbfns,Naobfns,Ispin)
         Else if (Abacus) Then

            Call Getlst(Qdp_ints,Ipert,1,1,Irrepx, 198)
            Call Expnd2(Qdp_ints,Qdp_ints_full,Naobfns)
            Call Ao2mo2(Qdp_ints_full,Qdp_ints_full,Scr1,Scr2,
     +                    Nbfns,Naobfns,Ispin)
         Endif
      write(6,"(a)") "The MO basis quadrapole integrals"
      call output(Qdp_ints_full,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

C The occupied/occupied block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jocc = 1, Pop(Irrepl,Ispin)
                  Oocount = Oocount +1
                  Qdp_ints_oo(Oocount) = Qdp_ints_full
     +                                           (Jocc+Ioffo(Irrepl),
     +                                            Iocc+Ioffo(Irrepr))
               Enddo
            Enddo
         Enddo 

C The virtual/virtual block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Ivrt = 1, Vrt(Irrepr,Ispin)
               Do Jvrt = 1,Vrt(Irrepl,Ispin)
                  Vvcount = Vvcount +1
                  Qdp_ints_vv(Vvcount) = Qdp_ints_full
     +                                           (Jvrt+Ioffv(Irrepl),
     +                                            Ivrt+Ioffv(Irrepr))
               Enddo
            Enddo
         Enddo 

C The virtual/occupied block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jvrt = 1, Vrt(Irrepl,Ispin)
                  Vocount = Vocount +1
                  Qdp_ints_vo(Vocount) = Qdp_ints_full
     +                                           (Jvrt+IoffV(Irrepl),
     +                                            Iocc+Ioffo(Irrepr))
               Enddo
            Enddo
         Enddo 

C Save the MO basis,sym-packed dipole integral blocks 

         Call Putlst(Qdp_ints_oo(Ioff_oo),1,1,1,Ispin,370)
         Call Putlst(Qdp_ints_vv(Ioff_vv),1,1,1,Ispin,371)
         Call Putlst(Qdp_ints_vo(Ioff_vo),1,1,1,Ispin,372)

      Write(6,*)
      Write(6,"(a)") "The OO-MO quadrapole matrices"
      write(6,"(6(1x,F12.6))") (Qdp_ints_oo(Ioff_oo-1+i),
     +                          i=1,Ndim_oo)
      Write(6,"(a)") "The VV-MO quadrapole matrices"
      write(6,"(6(1x,F12.6))") (Qdp_ints_vv(Ioff_vv-1+i),
     +                          i=1,Ndim_vv)
      Write(6,"(a)") "The VO-MO quadrapole matrices"
      write(6,"(6(1x,F12.6))") (Qdp_ints_vo(Ioff_vo-1+i),
     +                          i=1,Ndim_vo)

      Enddo
         
      Return
      End
