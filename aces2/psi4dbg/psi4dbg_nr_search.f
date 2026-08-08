













































































































































































































      Subroutine Psi4dbg_nr_search(Grd,Grd_oo,Grd_vv,Grd_vo,
     +                             Grd_ov,Grd_stat,Lenoo,Lenvv,Lenvo,
     +                             Nbas,Nocc,Nvrt,Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      
      Dimension Grd_stat(6)
      Dimension Grd(Nbas,Nbas)
      Dimension Grd_oo(Lenoo)
      Dimension Grd_vv(Lenvv)
      Dimension Grd_vo(Lenvo)
      Dimension Grd_ov(Lenvo)
      Dimension Work(Maxcor)



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
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end

      Data Ione,Izero,One,Half,Dnull,Two/1,0,1.0D0,0.5D0,0.0D0,
     +                                   2.0D0/
      Data Onem/-1.0D0/

      I000 = Ione
      I010 = I000 + Nbas
      I020 = I010 + Lenvo
      I030 = I020 + Lenvo
      Iend = I030 + Nbas*Nbas
      Memleft = Maxcor - Iend 

      Call Getrec(20,"JOBARC","SCFEVALA",Nbas*Iintfp,Work(I000))

      Indi   = Izero
      Inda   = Izero 
      Indv   = Izero 
      Indo   = Izero 
      Indv0  = Izero 
      Indo0  = Izero 
      Ijunk  = Izero
      Irrepx = Ione 

      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Nocci = Pop(Irrepr,1)
         Nvrti = Vrt(Irrepl,1)

         Do indxi = 1, Nocci
            Indi = Indi + 1
            Indv = Indv0
            Do indxa = 1, Nvrti
               Inda = Inda + 1
               Indv = Indv + 1
               Work(I010+Inda-1) = (Work(I000+Indv-1+Nocc) - 
     +                             Work(I000+Indi-1))
            Enddo
         Enddo
         Indv0 = Indv0 + Nvrti
      Enddo 

      Nsize = Inda
      Ioff = I010
      Joff = I020 
      Do Irrepr= 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Nocci = Pop(Irrepr,1)
         Nvrti = Vrt(Irrepl,1)
         Call Transp(Work(Ioff),Work(Joff),Nocci,Nvrti)
         Ioff = Ioff + Nocci*Nvrti
         Joff = Joff + Nocci*Nvrti
      Enddo 

      If (Lenvo .Ne. Nsize) Then
         Write(6,"(2a)") " There is an internal inconsistency in OV/VO",
     +                   " length."
         Call Errex
      Endif

      Write(6,"(a)") "The VO/OC orbital rotation hessians"
      Ioff = I010
      Joff = I020 
      Do Irrepr= 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Nocci = Pop(Irrepr,1)
         Nvrti = Vrt(Irrepl,1)
         Call output(work(ioff),1,nvrti,1,nocci,nvrti,nocci,1)
         Call output(work(joff),1,nocci,1,nvrti,nocci,nvrti,1)
         Ioff = Ioff + Nocci*Nvrti
         Joff = Joff + Nocci*Nvrti
      Enddo 

      Call Psi4dbg_vstat(Grd,Grd_stat,Nsize)
      Write(6,*)
      Write(6,"(2a)")  " The statistics of the orbital rotation",
     +                 " gradient matrix"
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"
      Write(6,"(5x,a,5xa,5xa)") "Minimum grad.", "Maximum grad.", 
     +                           "RMS grad"
      Write(6,*)
      Write(6,"(3(5x,E12.6))") Grd_stat(3), Grd_stat(4), Grd_stat(5)
      Write(6,*)
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"
      Write(6,*)

      If (Grd_stat(5) .Lt. Conv_tol) Then
          Write(6,"(2a)") " The micro-iterations of orbital rotation",
     +                    " gradients reached convergence" 
          Call Putrec("20","JOBARC","MICRO_CV",Iintfp,Ione)
      Endif 

C Note that I am using 1/(ea-ei) instead of 1/2(ea-ei) since I work with only
C the alpha (or beta block).  The main loop that rotates orbitals goes over 
C only alpha block only (in the RHF context). See pccd_rotg.F.
      
      Do I = 1, Nsize
         D = -One/Work(I010-1+I)
         Grd_vo(I) = Grd_vo(I)*D
         D = -One/Work(I020-1+I)
         Grd_ov(I) = Grd_ov(I)*D
      Enddo

      Call Dzero(Grd_oo,Lenoo)
      Call Dzero(Grd_vv,Lenvv)
      Call Pccd_frmful(Grd_oo,Grd_vv,Grd,Work(Iend),Memleft,Nbas,0)
      Call Pccd_frmful_ov(Grd_ov,Grd_vo,Grd,Work(Iend),Memleft,
     +                    Nbas,"COPY",0)
      Write(6,"(a)") "Scalled gradients (G=-H^-1g)" 
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)

C#ifdef _NOSKIP
C Form  G-G(^t). This appears double dipping but help convergence tremendoulsy. 

      Call Transp(Grd,Work(I030),Nbas,Nbas)
      Call Daxpy(Nbas*Nbas,Onem,Work(I030),1,Grd,1)

      write(6,*)
      Write(6,"(a)") "Antisymmetrized Scalled gradients (G-G(^t))" 
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Call Dgemm("N","T",Nbas,Nbas,Nbas,Half,Grd,Nbas,Grd,Nbas,
     +            Dnull,Work(I030),Nbas)

      Do I = 1, Nbas
         Grd(I,I) = Grd(I,I) + One
      Enddo

      Call Daxpy(Nbas*Nbas,One,Work(I030),1,Grd,1)

      Write(6,"(a)") "U=(1+K+1/2K*K)" 
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Call Pccd_gramschmidt(Grd,Nbas,Nbas)
C#endif 


      Write(6,"(a)") "Unitary check of Kappa" 
      Call Dgemm("N","T",Nbas,Nbas,Nbas,One,Grd,Nbas,Grd,Nbas,
     +           Dnull,Work(I030),Nbas)
      Call output(Work(i030),1,Nbas,1,Nbas,Nbas,Nbas,1)
   
      Return
      End

