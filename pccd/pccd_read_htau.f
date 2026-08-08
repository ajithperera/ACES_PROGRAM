










      Subroutine Pccd_form_htau(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Logical Sym_packed 
      Logical Symmetry
      Logical OO_constr,Opt_orbs 
      Character*5 Spin(2)

      Dimension Work(Maxcor)

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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
      Common /Symm/Symmetry 
      Common /OO_info/Opt_orbs,OO_constr

      Data Spin /"Alpha", "Beta "/
      Data Ione,Onem,One,Dnull,Half,Two/1,-1.0D0,1.0D0,0.0D0,0.50D0,
     +                                  2.0D0/

      Irrepx = Ione
      Nbas = Nocco(1) + Nvrto(1)

      Nocca   = Nocco(1)
      Nvrta   = Nvrto(1)
      Noccb   = Nocco(2)
      Nvrtb   = Nvrto(2)

      If (Symmetry) Then
         Ndim_ooa = Irpdpd(Irrepx,21)
         Ndim_vva = Irpdpd(Irrepx,19)
         Ndim_voa = Irpdpd(Irrepx,9)
         Ndim_oob = Irpdpd(Irrepx,22)
         Ndim_vvb = Irpdpd(Irrepx,19)
         Ndim_vob = Irpdpd(Irrepx,10)
      Else
         Ndim_ooa   = Nocca*Nocca
         Ndim_vva   = Nvrta*Nvrta
         Ndim_voa   = Nocca*Nvrta
         Ndim_oob   = Noccb*Noccb
         Ndim_vvb   = Nvrtb*Nvrtb
         Ndim_vob   = Noccb*Nvrtb
      Endif

      Ndim_oo = Ndim_ooa+Iuhf*Ndim_oob
      Ndim_vv = Ndim_vva+Iuhf*Ndim_vvb
      Ndim_vo = Ndim_voa+Iuhf*Ndim_vob
      Length  = Nbas*Nbas 

      I000 = Ione
      I010 = I000 + Ndim_oo
      I020 = I010 + Ndim_vv
      I030 = I020 + Ndim_vo
      I040 = I030 + Ndim_vo
      I050 = I040 + Length+Iuhf*Length
      Iend = I050 
      Maxcor = Maxcor - Iend 
      If (Iend.Gt.Maxcor) Call Insmem("Pccd_form_uhtau",Iend,Maxcor)

      Call Getrec(20,"JOBARC","ROTGRDOO",Ndim_oo*Iintfp,Work(I000))
      Call Getrec(20,"JOBARC","ROTGRDVV",Ndim_vv*Iintfp,Work(I010))
      Call Getrec(20,"JOBARC","ROTGRDVO",Ndim_vo*Iintfp,Work(I020))
      Call Getrec(20,"JOBARC","ROTGRDOV",Ndim_vo*Iintfp,Work(I030))
      Call Dscal(Ndim_vo,Onem,Work(I030),1)

      write(6,*)
      call checksum("RTGRDOOA:",Work(I000),Ndim_ooa)
      call checksum("RTGRDOOB:",Work(I000+Ndim_ooa),Ndim_oob)
      call checksum("RTGRDVVA:",Work(I010),Ndim_vva)
      call checksum("RTGRDVVB:",Work(I010+Ndim_vva),Ndim_vvb)
      call checksum("RTGRDVOB:",Work(I020),Ndim_voa)
      call checksum("RTGRDVOB:",Work(I020+Ndim_voa),Ndim_vob)
      call checksum("RTGRD0VB:",Work(I030),Ndim_voa)
      call checksum("RTGRDVVB:",Work(I030+Ndim_voa),Ndim_vob)
C
      If (Symmetry) Sym_packed = .True.
      Ioff    = I040
      Ioff_oo = I000
      Ioff_vv = I010
      Ioff_vo = I020
      Ioff_ov = I030

      Do Ispin = 1, Iuhf + 1 
         If (Ispin .Eq. 1) Then
            Nocc = Nocca
            Nvrt = Nvrta
         Else if (Ispin .Eq. 2) Then
            Nocc = Noccb
            Nvrt = Nvrtb
         Endif 
         Call Putrec(20,"JOABRC","SPNINDEX",Ione,Ispin)

         Call Pccd_frmful(Work(Ioff),Work(Ioff_oo),Work(Ioff_vv),
     +                    Work(Ioff_vo),Work(Ioff_ov),Work(Iend),
     +                    Maxcor,Nocc,Nvrt,Nbas,"Special",Sym_packed)
         Write(6,*)
         Write(6,"(a,5a)") "The orbital rotation gradients for spin: ",
     &                      Spin(Ispin)
         Call output(Work(Ioff),1,Nbas,1,Nbas,Nbas,Nbas,1)
         Call DScal(Nbas,Dnull,Work(Ioff),Nbas+1)
         Call Pccd_antisymmetrize(Work(Ioff),Work(Iend),Maxcor,Nbas,
     &                           Nocc,Nvrt)
         Call Pccd_htau_trnspose(Work(Ioff),Work(Iend),Nbas,Nocc,Nvrt)
 
        If (OO_constr) Call pccd_insist_no_oovv(Work(Ioff),Nocc,Nvrt,
     &                                          Nbas)

         Write(6,"(2a,5a)") " The anti-symm. orbital rotation gradient",
     +                     " for spin: ",Spin(Ispin)
         Call output(Work(Ioff),1,Nbas,1,Nbas,Nbas,Nbas,1)
         If (Ispin .Eq. 1) Call Putrec(20,"JOBARC","OBRTGRDA",
     +                                 Length*Iintfp,Work(Ioff))
         If (Ispin .Eq. 2) Call Putrec(20,"JOBARC","OBRTGRDB",
     +                                 Length*Iintfp,Work(Ioff))

         Ioff    = Ioff    + Nbas*Nbas
         Ioff_oo = Ioff_oo + Ndim_ooa
         Ioff_vv = Ioff_vv + Ndim_vva
         Ioff_vo = Ioff_vo + Ndim_voa
         Ioff_ov = Ioff_ov + Ndim_voa
      Enddo 

      Return
      End


