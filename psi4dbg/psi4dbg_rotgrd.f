













































































































































































































      Subroutine Psi4dbg_rotgrd(Work,Maxcor,Nbas,Nocc,Nvrt)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Tol,X,Fndlrgab
      Dimension Grad_stat(6)
      Dimension Work(Maxcor)
      Logical Gradcnv 



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

      Data Ione, Izero /1,0/

      Write(6,*) "-----------------pCCD_rotgrd-------------------"
      Irrepx   = Ione
      Length   = Nbas*Nbas
      Lenoo    = Irpdpd(Irrepx,21) 
      Lenvv    = Irpdpd(Irrepx,19) 
      Lenvo    = Irpdpd(Irrepx,9)  

      Tol=10.0**-(Iflags((76)))

      I000 = Ione
      I010 = I000 + Lenoo
      I020 = I010 + Lenvv
      I030 = I020 + Lenvo
      I040 = I030 + Lenvo
      Iend = I040 + Length
      Memleft = Maxcor - Iend 

      Call Getrec(20,"JOBARC","ORBRTGRD",Length*Iintfp,Work(I040))
      write(6,"(a,i2)") "The orbital rotation gradients",Nvrt
      call output(Work(I040),1,Nbas,1,Nbas,Nbas,Nbas,1)
      Call Pccd_prep_grad(Work(I040),Work(I000),Work(I010),
     +                    Work(I020),Work(I030),Work(Iend),
     +                    Memleft,Lenoo,Lenvv,Lenvo,Nocc,Nvrt,
     +                    Nbas)

       Call Psi4dbg_rotg(Work(I040),Work(I000),Work(I010),
     +                   Work(I020),Work(I030),Grad_stat,Work(Iend),
     +                   Memleft,Lenoo,Lenvv,Lenvo,Nocc,Nvrt,
     +                   Nbas,0)
C Grad_stat(1) : Larget absolute
C Grad_stat(2) : Smallest absolute
C Grad_stat(3) : Largest value
C Grad_stat(4) : Smallest value
C Grad_stat(5) : RMS gradient
C Grad_stat(6) : Dynamic range (abs. min-abs. max)

      Write(6,*)
      Write(6,90)
90    FORMAT(T3,'Determination of pCCD/CCD/LCCD optimum orbitals')
      Write(6,300)
300   Format(T3,"Orbital rotation gradients are used to determine",
     +          " new orbitals.")
      Y  = Grad_stat(1)
      X  = Grad_stat(5)
      Write(6,100) Y,X,Tol
100   FORMAT(T3,"Largest and RMS gradients and tolerance are: ",
     +           E12.6,1x,E12.6,1x,E12.6,'.')

      Call Getrec(-20,'JOBARC','ORBOPITR',Ione,Ncycle)
      If (X .Lt. Tol) Then
         Call Putrec(20,'JOBARC','GRADTEST',Ione,Ione)
         Write(6,200) Ncycle
200      Format(T3,"Pccd/CCD/LCCD gradient rotation converged in ", i3,
     +             " optimization cycles.")
         GradCNV = .TRUE.
      Else
         Call Putrec(20,'JOBARC',"GRADTEST",IONE,0)
         Gradcnv = .FALSE.
      Endif

      RETURN
      END
