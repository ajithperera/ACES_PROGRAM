










      Subroutine Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                             IHess_calc,Scale)
       
      Implicit Double Precision(A-H,O-Z)

      Logical Nonhf_ref,Symmetry
      Dimension Work(Maxcor)
      Dimension Nbfirr(8)

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

      Common/Nhfref/Nonhf_ref
      Common/Symm/Symmetry 

      Data Ione,Itwo,Dnull,Two /1,2,0.0D0,2.0D0/

      Nbas = Nocco(1) + Nvrto(1)

      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)

      If (Idrop .Ne. 0) Then
         Call Getrec(20,"JOBARC","NUMDROPA",1,Ndrop)
         Nbas = (Nbas - Ndrop) 
      Endif 

      Itrlen = Nbas*(Nbas+1)/Itwo 
      Isqlen = Nbas*Nbas 
      Nocc   = Nocco(1)
      Nvrt   = Nvrto(1)

      Noccsqf = Nocc*Nocc
      Nvrtsqf = Nvrt*Nvrt
      Nvrtocf = Nocc*Nvrt

      If (Symmetry) Then
         Irrepx = Ione
         Noccsqf = Irpdpd(Irrepx,21)
         Nvrtsqf = Irpdpd(Irrepx,19)
         Nvrtocf = Irpdpd(Irrepx,9)
      Endif 
   
      Ndim_oo = Nocc*Nocc
      Ndim_vv = Nvrt*Nvrt
      Ndim_vo = Nvrt*Nocc

      I000 = Ione
      I010 = I000 + Isqlen
      I020 = I010 + Isqlen
      I030 = I020 + Isqlen
      I040 = I030 + Isqlen
      I050 = I040 + Isqlen

      I060 = I050 + Ndim_oo
      I070 = I060 + Ndim_vv
      I080 = I070 + Isqlen
      I090 = I080 + Isqlen 
      Iend = I090 + Ndim_vo
      Maxcor = Maxcor - Iend 
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_htau_pq",Iend,Maxcor)

      Call Getrec(20,"JOBARC",'FOCKA   ',Isqlen*IIntfp,Work(I020))
      Call Ao2mo2(Work(I020),Work(I020),Work(I070),Work(I080),Nbas,
     +            Nbas,1)

      Write(6,"(a)") "The MO basis FOCK matrix"
      call output(Work(I020),1,nbas,1,nbas,nbas,nbas,1)
C 
      Call Getrec(20,"JOBARC","DENSOO  ",Noccsqf*IIntfp,Work(I050))
      Call Getrec(20,"JOBARC","DENSVV  ",Nvrtsqf*IIntfp,Work(I060))
      Write(6,"(a,1x,l)") "Is nonhf set?",nonhf_ref
      If (Nonhf_ref) Call Getrec(20,"JOBARC","DENSVO  ",
     +                           Nvrtocf*IIntfp,Work(I090))

      Call Pccd_form_dpq(Work(I030),Work(I040),Work(I050),Work(I060),
     +                   Work(I090),Work(Iend),Maxcor,Nocc,Nvrt,Nbas,
     +                   Iuhf,Nonhf_ref)

C Built H_tau(pq) = h(r,p)*D(q,r). Note that D(o,v)=D(v,o) = 0 (i.e.occ-vrt or 
C vrt-occ block of H_tau is zero). 

      Call Pccd_form_htau_1D_1(Work(I000),Work(I010),Work(Iend),
     +                         Maxcor,Work(I020),Work(I030),
     +                         Work(I040),Nocc,Nvrt,Nbas,Iuhf)
      If (Symmetry) Then
         Call Pccd_symexp(Work(I050),Work(I060),Work(I090),
     +                    Work(Iend),Maxcor,Nocc,Nvrt,Nbas)
      Endif 

      Call Pccd_form_htau_1D_2(Work(I000),Work(I010),Work(I050),
     +                         Work(I060),Work(I090),Work(I040),
     +                         Work(Iend),Maxcor,Nocc,Nvrt,Nbas,
     +                         Iuhf)


      I050 = I040 + Max(Noccsqf,Nvrtocf)
      I060 = I050 + Max(Noccsqf,Nvrtocf)
      I070 = I060 + Nvrtsqf
      I080 = I070 + Nvrtsqf
      I090 = I080 + Nvrtocf
      I100 = I090 + Nvrtocf
      I110 = I100 + Isqlen
      Iend = I110 + Isqlen
      Maxcor = Maxcor - Iend 
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_htau_pq",Iend,Maxcor)

      Call Pccd_form_htau_oo(Work(I050),Work(I060),Work(I090),
     +                       Work(I040),Work(I050),Work(Iend),
     +                       Maxcor,Nocc,Nvrt,Noccsqf,
     +                       Nvrtsqf,Nvrtocf,Iuhf,Nonhf_ref)

      Call Pccd_form_htau_ov(Work(I050),Work(I060),Work(I090),
     +                       Work(I040),Work(I050),Work(Iend),
     +                       Maxcor,Nocc,Nvrt,Noccsqf,
     +                       Nvrtsqf,Nvrtocf,Iuhf,Nonhf_ref)

      Call Pccd_form_htau_2D(Work(I000),Work(I010),Work(I040),
     +                       Work(I050),Work(I060),Work(I070),
     +                       Work(I080),Work(I090),Work(I030),
     +                       Work(I100),Work(I110),Work(Iend),
     +                       Nocco(1),Nvrto(1),Maxcor,Nbas,E)

      Call Putrec(20,"JOBARC","OBRTGRDA",Isqlen*IIntfp,Work(I000))

      Write(6,*)
      Call Pccd_sortgam(Work(Iend),Maxcor,Iuhf)
      Call Pccd_gfock_debug(Work(Iend),Maxcor,Iuhf)
      If (IHess_calc .EQ. 0) Return

C This is part of the Hessian that is computed from H*{p+q}. This
C is done here since all the pieces that is needed for this (Htau_pq 
C D_pq, f_pq) is available. 

      I100 = Iend
      I110 = I100 + Noccsq
      I120 = I110 + Noccsq
      I130 = I120 + Noccsq
      I140 = I130 + Nvrtsq
      I150 = I140 + Nvrtsq
      I160 = I150 + Nvrtsq
      I170 = I160 + Nvrtoc
      I180 = I170 + Nvrtoc
      I190 = I180 + Nvrtoc
      I200 = I190 + Nvrtoc
      I210 = I200 + Nvrtoc
      Iend = I210 + Nvrtoc
      Maxcor = Maxcor - Iend 
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_htau_pq",Iend,Maxcor)

      Call Pccd_form_htau_1D_pqrs(Work(I000),Work(I010),
     +                            Work(I020),Work(I030),
     +                            Work(I040),Work(I050),
     +                            Work(I060),Work(I070),
     +                            Work(I080),Work(I090),
     +
     +                            Work(I100),Work(I010),
     +                            Work(I120),Work(I130),
     +                            Work(I140),Work(I150),
     +                            Work(I160),Work(I170),
     +                            Work(I180),Work(I190),
     +                            Work(I200),Work(I210),
     +                            Work(Iend),Nocco(1),
     +                            Nvrto(1),Maxcor,Nbas)

      Return 
      End 
    
     
