










      Subroutine Pccd_nodiis(Grd,Work,Maxcor,Ispin,Nbas,Ncycle)

      Implicit Double Precision(A-H,O-Z)
    
      Dimension Work(Maxcor)
      Dimension Grd(Nbas,Nbas)


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
      Common /Symm/Symmetry

      Data Ione,Onem,half,Dnull,One/1,-1.0D0,0.50D0,0.0D0,1.0D0/

C Form  K-K(^t). This appears double dipping but help convergence
C tremendoulsy.
      Nbas2 = Nbas*Nbas
      I000 = Ione
      I010 = I000 + Nbas2
      I020 = I010 + Nbas2
      I030 = I020 + Nbas2
      I040 = I030 + Nbas2
      I050 = I040 + Nbas2
      
      If (Iend .Gt. Maxcor) Call Insmem("pccd_nodiis",Iend,Maxcor)
      Maxcor = Maxcor-Iend
    
CSSS      Call Transp(Grd,Work(I000),Nbas,Nbas)
CSSS      Call Daxpy(Nbas*Nbas,Onem,Work(I000),1,Grd,1)

      write(6,*)
      Write(6,"(a)") "Antisymmetrized Scalled gradients (K-K(^t))"
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)

C 1/2 K*K(^t)
      Call Dgemm("N","T",Nbas,Nbas,Nbas,Half,Grd,Nbas,Grd,Nbas,
     +            Dnull,Work(I000),Nbas)

      Do I = 1, Nbas
         Grd(I,I) = Grd(I,I) + One
      Enddo

      Call Daxpy(Nbas*Nbas,One,Work(I000),1,Grd,1)

      Write(6,"(a)") "U=(1+G+1/2G*K)"
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)

      Call Pccd_gramschmidt(Grd,Nbas,Nbas)

      Write(6,"(a)") "Unitary check of Kappa"
      Call Dgemm("N","T",Nbas,Nbas,Nbas,One,Grd,Nbas,Grd,Nbas,
     +           Dnull,Work(I000),Nbas)
      Call output(Work(i000),1,Nbas,1,Nbas,Nbas,Nbas,1)

      If (Ispin .Eq. 1) Call Getrec(20,"JOBARC","SCFEVECA",Nbas2*Iintfp,
     +                              Work(I010))
      If (Ispin .Eq. 2) Call Getrec(20,"JOBARC","SCFEVECB",Nbas2*Iintfp,
     +                              Work(I010))


      If (Ispin .Eq. 1) Call Getrec(20,"JOBARC","EVECOAOA",Nbas2*Iintfp,
     +                              Work(I010))
      If (Ispin .Eq. 2) Call Getrec(20,"JOBARC","EVECOAOB",Nbas2*Iintfp,
     +                              Work(I010))

      If (Ispin .Eq. 1) Call Putrec(20,"JOBARC","EVCOAOXA",Nbas2*Iintfp,
     +                              Work(I010))
      If (Ispin .Eq. 2) Call Putrec(20,"JOBARC","EVCOAOXB",Nbas2*Iintfp,
     +                              Work(I010))

      Call Mo2oao(Grd,Work(I000),Work(I010),Work(I020),Nbas,Ispin)
      call output(Work(i000),1,Nbas,1,Nbas,Nbas,Nbas,1)

      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(I000),Nbas,
     +            Work(I010),Nbas,Dnull,Work(I020),Nbas)

      If (Ispin .Eq. 1) Call Putrec(20,"JOBARC","SCFEVCA0",Nbas2*Iintfp,
     +                              Work(I020))
      If (Ispin .Eq. 2) Call Putrec(20,"JOBARC","SCFEVCB0",Nbas2*Iintfp,
     +                              Work(I020))
      Write(6,"(a)") "Orthogonalized new rotated Vectos in OAO basis"
      call output(Work(i020),1,Nbas,1,Nbas,Nbas,Nbas,1)
      call Dgemm("T","N",Nbas,Nbas,Nbas,One,Work(I020),Nbas,Work(I020),
     +            Nbas,Dnull,Work(I030),Nbas)
      Ones = Ddot(Nbas*Nbas,Work(I030),1,Work(I030),1)
      If (Int(Ones) .Eq. Nbas) Write(6,"(2a)") "Normalization test",
     +                                         " passes" 
CSSS      Call Pccd_rotate(Grd,Work(Iend),Maxcor,Ispin,Nbas,Ncycle)

      Return
      End
