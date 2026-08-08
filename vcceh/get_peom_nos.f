










      Subroutine Get_peom_nos(Work,Maxcor,Iwork,Imaxcor,Iuhf,Npert,
     +                        IPert,Numpert,Irrepx,Iside,Build_dens,
     +                        Save_vecs) 

      Implicit Double Precision(A-H,O-Z) 

      Dimension Work(Maxcor)
      Dimension IWork(Imaxcor)
      Dimension Npert(8)
      Logical Realfreq
      Logical Build_dens
      Logical SS,SD,DS,DD
      Logical CCSD,Parteom,Nodavid
      Logical Save_vecs
      Dimension Lenvv(2),Lenoo(2)



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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      Common/Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common/Drvhbar/SS,SD,DS,DD
      Common/Eominfo/CCSD,PARTEOM,NODAVID

      Data Ione,Dnull,Inull,Icen /1,0.0D0,0,100/
  

C The variables that keeps active space info is defined as having a
C maximum of 100 perturbation per irrep (see active_space.com). 
C Let make sure that we stay in that limit.

      If (Numpert .Gt. Icen) Then
          Write(6,"(a,I3,a,I3,a)") " The number of perturbations ",
     +                 Numpert, "in irrep ", Irrepx, " is greater" 
          Write(6,"(a,I3)") " than the maximum allowed ", Icen
          Call Errex
      Endif 

C The Maxexp is maximum number of iterations in linear equation
C solver.

      Freq     = 0.0D0
      Realfreq = .True.
      DD       = .False.
      Parteom  = .True.
      Nodavid  = .True.
      
      N1aa = Irpdpd(Irrepx,9)
      ndim = N1aa
      If (Iuhf .Ne. 0) Ndim = Ndim + Irpdpd(Irrepx,10)
      Kmax = Min(N1aa,Maxexp-1)
      
      I000 = Ione 
      I010 = I000 + Kmax*Ndim
      I020 = I010 + Kmax*Ndim
      I030 = I020 + Numpert*Ndim
      I040 = I030 + Numpert*Ndim 
      I050 = I040 + Ndim 
      Iend = I050 + Numpert
     
      If (Iend .Ge. Maxcor) Call Insmem("get_peom_nos",Iend,Maxcor)

      Memleft = Maxcor - Iend
      Imax = Iend + (2*Nsizec+Nsizec/4)

      If (Imax .Ge. Memleft) Call Insmem("get_peom_nos",Iend,Memleft)

      Call Blcklineq(Iuhf,Irrepx,Iside,Numpert,Kmax,N1aa,Ndim,
     +               Work(I030),Work(I020),Work(I000),Work(I010),
     +               Work(I040),Work(I050),Work(Iend),Memleft,
     +               Freq,Realfreq,Save_vecs)

      If (.Not. Build_dens .And. .Not. Save_vecs) Then
         Parteom = .False.
         Nodavid = .False.
         DD      = .True.
         Call Dcopy(Nsizec,Work(Iend),1,Work(I000),1)
         Return
      Endif 

      Do Ipert = 1, Numpert 

         Ioffset = Ipert
 
         Call setup_4nos(Work(I000),Maxcor,Irrepx,Ioffset,Iside,Iuhf)

         Do Ispin = 1, Iuhf+1
            Lenvv(Ispin) = Irpdpd(1,18+Ispin)
            Lenoo(Ispin) = Irpdpd(1,20+Ispin)
         Enddo

         I000 = Ione 
         I010 = I000 + Lenoo(1)+Iuhf*Lenoo(2)
         I020 = I010 + Lenvv(1)+Iuhf*Lenvv(2)
         I030 = I020 + Lenoo(1)+Iuhf*Lenoo(2)
         Iend = I030 + Lenvv(1)+Iuhf*Lenvv(2)

         If (Iend .Ge. Maxcor) Call Insmem("get_peom_nos",Iend,Maxcor)
         Maxcor = Maxcor - Iend
 
         Call Get_perturb_rspns_dens(Work(I000),Work(I010),Work(Iend),
     +                               Maxcor,Irrepx,Iuhf,Lenvv,Lenoo)
      
         Call Get_perturb_nos(Work(I000),Work(I010),Work(I020),
     +                        Work(I030),Work(Iend),Maxcor,Iwork,
     +                        Imaxcor,Lenvv,Lenoo,Ioffset,Iuhf)

C Save the OO and VV natural orbitals on disk. 

         Call Putlst(Work(I020),Ioffset,1,1,Irrepx,380)
         Call Putlst(Work(I030),Ioffset,1,1,Irrepx,381)

         If (Iuhf .Ne. 0) Then

            Call Putlst(Work(I020+Lenoo(1)),Ioffset,1,1,Irrepx,
     +                       382)
            Call Putlst(Work(I030+Lenvv(1)),Ioffset,1,1,Irrepx,
     +                   383)
         Endif 
      Enddo 

      Parteom = .False.
      Nodavid = .False.
      DD      = .True.

      Return
      End 
  
       
