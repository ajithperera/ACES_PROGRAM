













































































































































































































      Subroutine Prep_scfmos(Work,Maxcor,Coo,Cvv,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Logical Dropmo



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
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Dimension Work(Maxcor)
      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
 
      Data Ione,Inull /1,0/

      Call Getrec(20,"JOBARC","NAOBASFN",Ione,Naobfns)
      Call Getrec(20,"JOBARC","NBASTOT ",Ione,Nbfns_tot)
      Call Getrec(0,"JOBARC","NUMDROPA",Length,Ndrop)
  
      Dropmo = (Length .Gt. Inull)
         
      If (Dropmo) Then 
         Call Getrec(20,"JOBARC","NUMDROPA",Ione,Ndrop)
         If (Iuhf .Ne. 0) Call Getrec(20,"JOBARC","NUMDROPB",Ione,
     +                                Ndrop)
      Endif 

      Nbfns = Nbfns_tot - Ndrop 
      If (Nbfns .Ne. Nocco(1)+nvrto(2)) Then
          Write(6,"(2a)") " Internal inconsistency in number of basis",
     +                    " functions."
          Write(6,"(a,I3,a,I3)") " Read: ",Nbfns, " Internally: ", 
     +                            Nocco(1)+nvrto(1)
      Endif 
 
      Nbfns2 = Nbfns_tot*Nbfns
      I000 = Ione
      I010 = 1000 + (Iuhf+1)*Nbfns2
      I020 = 1010 + (Iuhf+1)*Nbfns
      Iend = I020 + Nbfns2
       
      Call Getrec(20,"JOBARC","SCFEVECA",Nbfns2,Work(I000))
      If (Iuhf .Ne.0) Call Getrec(20,"JOBARC","SCFEVECB",NBfns2,
     +                            Work(I000+Nbfns2))

      Call Getrec(20,"JOBARC","SCFEVALA",Nbfns,Work(I010))
      If (Iuhf .Ne.0) Call Getrec(20,"JOBARC","SCFEVALB",NBfns,
     +                            Work(I010+Nbfns))

      Write(6,"(a)") " The alpha SCF eigenvalues"
      Write(6,"(5(1x,F15.5))") (Work(I010+i-1),i=1,Nbfns)
      Write(6,"(a)") " The alpha SCF vectors"
      Call output(Work(I000),1,Nbfns_tot,1,Nbfns,Nbfns_tot,Nbfns,1)
      If (Iuhf .Ne. 0) Then
      Write(6,"(a)") " The beta SCF eigenvalues"
        Write(6,"(5(1x,F15.5))") (Work(I010+Nbfns+i-1),i=1,Nbfns)
         Write(6,"(a)") " The beta SCF vectors"
         Call output(Work(I000+Nbfns2),1,Nbfns_tot,1,Nbfns,Nbfns_tot,
     +               Nbfns,1)
      Endif 


      Return
      End
